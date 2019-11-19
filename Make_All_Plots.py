params = {"figure.figsize": (12,9),
          "font.size": 20,
          "font.weight": "normal",
          "xtick.major.size": 9,
          "xtick.minor.size": 4,
          "ytick.major.size": 9,
          "ytick.minor.size": 4,
          "xtick.major.width": 4,
          "xtick.minor.width": 3,
          "ytick.major.width": 4,
          "ytick.minor.width": 3,
          "xtick.major.pad": 8,
          "xtick.minor.pad": 8,
          "ytick.major.pad": 8,
          "ytick.minor.pad": 8,
          "lines.linewidth": 3,
          "lines.markersize": 10,
          "axes.linewidth": 4,
          "legend.loc": "best",
          "text.usetex": False,    
          "xtick.labelsize" : 20,
          "ytick.labelsize" : 20,
          }

import matplotlib
matplotlib.rcParams.update(params)
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.offsetbox import AnchoredText

import numpy as np

import astropy as ap
import astropy.coordinates as coord
import astropy.units as u
from astropy.time import Time


### MKT J170456.2-482100 paper plots ###
# Use this script to make Figures 2, 4, 5, 6, 8, 10, 11 and 12
# from https://doi.org/10.1093/mnras/stz3027

### Required files:
# TYC_optical_semesters.npy
# TYC_LS_periods.npy
# TYC_LS_periodErrors.npy
# TYC_optical_binned.npy


def phase_up(ref_mjd, mjds, period):
    '''
    Phase up and fold a light curve
    
    Folds the MJDs of a light curve to a
    certain period (in days)
    and phases the MJDs to a specific MJD.
    
    Args:
    ref_mjd (float): the MJD that you are using as
                     your reference. Your light curve will
                     be phased to this MJD
    mjds (array): an array of floats, the MJDs of
                  your light curve
    period (float): the period (in days) that you are
                    folding to

    Returns:
    An array corresponding to the phase of each
    MJD after it has been folded and phased to
    ref_mjd.
    '''
    return (((mjds - ref_mjd)/period)-
            ((mjds - ref_mjd)//period))


def weighted_mean(values, uncertainties):
    '''
    Find the weighted mean of an array of values.
    
    Args:
    values (array): the array of values that you would
                    like to find the weighted mean of
    uncertainties (array): the array of uncertainties
                           corresponding to the values.
                           These will be used as the weights
    Returns:
    mean (float): the weighted mean of the values
    error (float): the uncertainty on the weighted mean
    '''
    # Find the weight of each value using the uncertainties
    weights = 1./(uncertainties**2)
    
    # Find the weighted mean of the values
    mean = (np.sum(values*weights))/(np.sum(weights))
    # Find the uncertainty on the weighted mean
    error = 1./(np.sqrt(np.sum(weights)))
    
    return mean, error


# Save all of the Figures
save_figures = True

# Set up the colour-scheme that will be used in all of the plots
season_colours = ['#4363d8', '#911eb4',
                '#f032e6', '#e6194b',
                '#f58231', '#ffe119']
season_colours *= 6
# Set up the markers that you'll use throughout the plots
season_markers = ['o', 's', '^', 'd']
season_markers *= 6

# Load in the optical observations that have been split into
# into semesters, plus the results of the Lomb-Scargle analysis
all_optical_info = np.load('TYC_optical_semesters.npy')
optical_mjd_chunks = all_optical_info[0]
optical_mag_chunks = all_optical_info[1]
optical_mag_e_chunks = all_optical_info[2]
optical_name_chunks = all_optical_info[3]

mean_all_V = np.mean(optical_mag_chunks[-1])

periods = np.load('TYC_LS_periods.npy')
period_errors = np.load('TYC_LS_periodErrors.npy')

# Sort the MJDs into chronological order
first_mjds = []
for i, time in enumerate(optical_mjd_chunks[:-1]):
    first_mjds.append(int(np.min(time)))
first_mjds = np.array(first_mjds)
time_order = np.argsort(first_mjds)

### Make Figure 4 of the article ###

# Do you want to see what the period is
# for each semester? (False means no)
print_periods = False

# Set the frame colour and data point size
framecol= 'Black'
pointsize = 3

fig, ax = plt.subplots(1, 1, figsize=(8, 8))

for i, j in enumerate(time_order):
    if print_periods:
        # Print the information for each semester in order
        print((u'{0} \tperiod: {1:.4f} '
                '\u00B1 {2:.4f} days').format(optical_name_chunks[j],
                                                periods[j],
                                                period_errors[j]))
    # Make sure that each instrument has the
    # same marker style every time
    if optical_name_chunks[j]=='ASAS':
        m = 'o'
    elif optical_name_chunks[j]=='KELT':
        m = 'd'
    elif optical_name_chunks[j]=='ASAS-SN V':
        m = 'v'
    else:
        m = 'p'

    # Plot the period and uncertainty for the semester
    ax.errorbar(periods[j],
                int(np.min(optical_mjd_chunks[j])),
                xerr=period_errors[j],
                c=season_colours[j],
                marker=m,
                markeredgecolor='Grey')

    if periods[j] < 16:
        # One KELT season has more power in an harmonic,
        # so plot the real period here too.
        ax.errorbar(periods[j]*2,
                    int(np.min(optical_mjd_chunks[j])),
                    xerr=period_errors[j], fmt=m,
                    c=framecol, markeredgecolor='Grey')

# Plot the best period (from the Lomb-Scargle of
# all the V-band observations put together) as
# a grey dashed line
ax.plot(np.ones(10)*periods[-1],
        np.linspace(np.min(optical_mjd_chunks[-1])-500,
                    np.max(optical_mjd_chunks[-1])+500,
                    10), '--', c='Grey', zorder=0)

# Set all the ticks and spines to the frame colour
ax.tick_params(color=framecol, labelcolor=framecol)
for spine in ax.spines.values():
    spine.set_edgecolor(framecol)

# Label the left y-axis with the first MJD
# in each semester. I miss some semesters
# here so that the MJDs don't overlap
first_mjds = np.array([51933, 52440,
                        52664, 53036,
                        53404, 53759,
                        54136, 54501,
                        54869, 56427,
                        56705, 57084,
                        57457, 57777,
                        58150,
                        58508])
plt.yticks(first_mjds, first_mjds)
ax.set_xlabel('Period (days)', fontsize=18, color=framecol)

# Label the right y-axis with the date corresponding
# to the first MJD in each semester that is shown
# on the left y-axis
ymds = []
for date in Time(first_mjds.astype(float), format='mjd').isot:
    ymd = date[:date.index('T')]
    ymds.append(ymd)
ax2 = ax.twinx()
plt.yticks(first_mjds, ymds)

# Set the y-lim on both sides so everything
# lines up nicely
ax.set_ylim(np.min(optical_mjd_chunks[-1])-200,
            np.max(optical_mjd_chunks[-1])+300)
ax2.set_ylim(np.min(optical_mjd_chunks[-1])-200,
            np.max(optical_mjd_chunks[-1])+300)

ax.set_ylabel('MJD', fontsize=18, color=framecol)

# Make a legend that only shows the shape of
# the markers for each different instrument
asas_dot = mlines.Line2D([], [], color='None', marker='o',
                            markeredgecolor='DarkGrey',
                            markersize=pointsize+5, markeredgewidth=3,
                            label='ASAS')
kelt_dot = mlines.Line2D([], [], color='None', marker='d',
                            markeredgecolor='DarkGrey',
                            markersize=pointsize+5, markeredgewidth=3,
                            label='KELT')
asassnV_dot = mlines.Line2D([], [], color='None', marker='v',
                            markeredgecolor='DarkGrey',
                            markersize=pointsize+5, markeredgewidth=3,
                            label='ASAS-SN V band')
asassng_dot = mlines.Line2D([], [], color='None', marker='p',
                            markeredgecolor='DarkGrey',
                            markersize=pointsize+5, markeredgewidth=3,
                            label='ASAS-SN g band')
# Put the legend outside of the plot
leg0 = ax.legend(handles=[asas_dot, kelt_dot, asassnV_dot, asassng_dot],
                    fontsize=15, frameon=True, loc='lower left', ncol=2,
                    bbox_to_anchor= (0.0, 1.01), borderaxespad=0)

# Flip the y-axes
ax.invert_yaxis()
ax2.invert_yaxis()

plt.tight_layout()

if save_figures:
    # Save the figure how you like, but you do need
    # bbox_inches='tight' so that the legend isn't cut
    plt.savefig(('Optical_LombScargleResults_'
                    'FrameColour{}.eps').format(framecol),
                transparent=False,
                bbox_inches='tight')
plt.close()

### Make Figure 5 of the article ###

# The MJD that you're going to phase everything up to:
zero = 53571.094339999836
# The best period (from performing a Lomb-Scargle
# analysis on all of the V-band observations together)
best_period = periods[-1]
print('Period: {:.6f} days'.format(best_period))

# Get the binned optical data
binned_optical_data = np.load('TYC_optical_binned.npy')
# The MJDs of the data pints
binned_times = binned_optical_data[0]
# The magnitudes
binned_mags = binned_optical_data[1]
# The errors on the magnitudes
binned_errs = binned_optical_data[2]

# Make the plot

# Set the point colour
pointcolor = 'Black'
# Set the point size
pointsize = 10

# Set the number of subplots so that you have one
# per semester
fig, ax = plt.subplots(len(binned_times)//2, 2,
                        figsize=(14, 20))

for i, times in enumerate(binned_times):
    # Make sure that you keep the symbols consistent
    # for each instrument
    if optical_name_chunks[i]=='ASAS':
        m = 'o'
    elif optical_name_chunks[i]=='KELT':
        m = 'd'
    elif optical_name_chunks[i]=='ASAS-SN V':
        m = 'v'
    else:
        m = 'p'
    # This is so that the sub plots are in
    # order.
    if i < len(binned_times)//2:
        c = 0
        r = i
    else:
        c = 1
        r = i - (len(binned_times)//2)

    # Make a label that has the MJD
    # and ISOT date of the first data point
    # in the semester
    min_mjd = np.int(np.min(times))
    min_date = Time(min_mjd, format='mjd')
    min_date = (min_date.isot).split('T')[0]        
    info_text = AnchoredText((u'First MJD: {0}\n{1}').format(min_mjd, min_date),
                                loc='lower right', prop=dict(fontsize=14))
    # Be consistent with your colours for each season
    if len(optical_mag_chunks[i])>35:
        col = season_colours[i]
    else:
        # Some seasons are included that didn't have
        # enough data points to perform a Lomb-Scargle.
        # Colour these seasons grey
        col = 'Grey'

    btimes = np.array(binned_times[i])
    bmags = np.array(binned_mags[i])
    berrs = np.array(binned_errs[i])

    # Plot two folded periods by phasing the data
    # to the specified MJD and folding
    # to the best period
    ax[r, c].errorbar(phase_up(zero, btimes, best_period),
                        bmags,
                        yerr=berrs,
                        c=col, fmt=m,
                        markersize=pointsize, markeredgecolor='DarkGrey')
    ax[r, c].errorbar(phase_up(zero, btimes, best_period)+1,
                        bmags,
                        yerr=berrs,
                        c=col, fmt=m,
                        markersize=pointsize, markeredgecolor='DarkGrey')


    if np.int(np.min(times)) == 58150:
        ax[r, c].set_ylim(10.88, 11.06)

    # Since we're working with magnitudes, flip
    # the y-axis
    ax[r, c].invert_yaxis()

    ax[r, c].set_xlim(0, 2)
    # Add the label
    ax[r, c].add_artist(info_text)

    # Unless the plot is one of the bottom two plots,
    # you don't want x-tick labels, and you wan the
    # x-ticks to be inside the plot
    if ((r, c) != (len(binned_times)//2, 0) and
        (r, c) != (len(binned_times)//2 - 1, 0) and
        (r, c) != (len(binned_times)//2 - 1, 1)):
        ax[r, c].axes.get_xaxis().set_ticklabels([])
        ax[r, c].tick_params(axis='x', direction='in')
    # Give the plots on the left-hand column y-labels
    if c == 0:
        ax[r, c].set_ylabel('Apparent'+'\n'+'Magnitude', fontsize=18)

# For the bottom two plots, keep the x-tick-labels and
# label the x-axis
ax[len(binned_times)//2 - 1, 0].set_xlabel('Phase', fontsize=20)
ax[len(binned_times)//2 - 1, 1].set_xlabel('Phase', fontsize=20)

fig.tight_layout()
fig.subplots_adjust(hspace=0)

if save_figures:
    plt.savefig('Optical_FoldedSemesters.eps')
plt.close()

### Make Figure 11 from the paper ###

fig, ax = plt.subplots(figsize=(15, 8))

pointsize = 10
for i, time in enumerate(optical_mjd_chunks[:-1]):
    # Make sure that you keep the symbols consistent
    # for each instrument
    if optical_name_chunks[i]=='ASAS':
        m = 'o'
    elif optical_name_chunks[i]=='KELT':
        m = 'd'
    elif optical_name_chunks[i]=='ASAS-SN V':
        m = 'v'
    else:
        m = 'p'

    # Be consistent with your colours for each season
    if len(optical_mag_chunks[i])>35:
        col = season_colours[i]
    else:
        # Some seasons are included that didn't have
        # enough data points to perform a Lomb-Scargle.
        # Colour these seasons grey
        col = 'Grey'

    # Calculate the weighted mean and uncertainties on the weighted mean
    # for each semester
    wmean, werror = weighted_mean(optical_mag_chunks[i],
                                    optical_mag_e_chunks[i])
    (_, caps, _) = ax.errorbar((((np.max(time)-np.min(time))/2.)+
                                np.min(time)),
                                wmean,
                                yerr=np.std(optical_mag_chunks[i]),
                                c=col, fmt=m,
                                markeredgecolor='DarkGrey',
                                capsize=20)
    # Change the size of the lines on top of
    # the "errorbars"
    for cap in caps:
        cap.set_markeredgewidth(3)

ax.set_xlim(51800, 58850)

# Make some custom x-ticks so you can put
# the months on the top of the plot
ax2 = ax.twiny()
ax2.set_xticks([51969, 53065, 54160,
                55256, 56352, 57448, 58543])
ax2.set_xticklabels(['Mar 2001', 'Mar 2004',
                        'Mar 2007', 'Mar 2010',
                        'Mar 2013', 'Mar 2016',
                        'Mar 2019'])
ax2.set_xlim(51800, 58850)

# Make a custom legend
asas_dot = mlines.Line2D([], [], color='None', marker='o',
                            markeredgecolor='DarkGrey',
                            markersize=pointsize+2, markeredgewidth=3,
                            label='ASAS')
kelt_dot = mlines.Line2D([], [], color='None', marker='d',
                            markeredgecolor='DarkGrey',
                            markersize=pointsize+2, markeredgewidth=3,
                            label='KELT')
asassnV_dot = mlines.Line2D([], [], color='None', marker='v',
                            markeredgecolor='DarkGrey',
                            markersize=pointsize+2, markeredgewidth=3,
                            label='ASAS-SN V band')
asassng_dot = mlines.Line2D([], [], color='None', marker='p',
                            markeredgecolor='DarkGrey',
                            markersize=pointsize+2, markeredgewidth=3,
                            label='ASAS-SN g band')
leg0 = ax.legend(handles=[asas_dot, kelt_dot, asassnV_dot, asassng_dot],
                    fontsize=26, frameon=True, loc='lower left', ncol=2)
# Because we're working with magnitudes, invert
# the y-axis
ax.invert_yaxis()

ax.set_ylabel('Apparent magnitude', fontsize=30)
ax.set_xlabel('MJD', fontsize=30)

plt.tight_layout()

if save_figures:
    plt.savefig('Optical_StandDevs.eps')
plt.close()

### Make Figure 8 from the paper ###

# Get the radial velocities
SALT_values = np.load('SALT_radial_velocities.npy')
LCO_values = np.load('LCO_radial_velocities.npy')

SALT_rv_jds = SALT_values[:, 0]
SALT_rv_times = Time(SALT_rv_jds, format='jd')
SALT_rv_mjds = SALT_rv_times.mjd

SALT_rvs = SALT_values[:, 1]
SALT_rv_errs = SALT_values[:, 2]

LCO_rv_jds = LCO_values[:, 0]
LCO_rv_times = Time(LCO_rv_jds, format='jd')
LCO_rv_mjds = LCO_rv_times.mjd

LCO_rvs = LCO_values[:, 1]
LCO_rv_errs = LCO_values[:, 2]

# Plot
# Set the pointsize
pointsize = 10

# We want to plot just the most recent g-band
# ASAS-SN semester.
i = -2
times = optical_mjd_chunks[i]
mags = np.copy(optical_mag_chunks[i])
mag_e = np.copy(optical_mag_e_chunks[i])

fig, ax1 = plt.subplots(1, 1, figsize=(12, 8))

# Make a a second y-axes to plot
# the RVs and optical magnitudes
# together
ax2 = ax1.twinx()

# On the first axis, plot the folded and phased
# optical magnitudes (two full periods)
ax1.errorbar(phase_up(zero, times, best_period),
                mags, yerr=mag_e,
                fmt='o', c='DarkGrey')
ax1.errorbar(phase_up(zero, times, best_period)+1,
                mags, yerr=mag_e,
                fmt='o', c='DarkGrey')

# Put the SALT and LCO radial velocities together
rv_mjds = np.concatenate((SALT_rv_mjds, LCO_rv_mjds))
rv_vals = np.concatenate((SALT_rvs, LCO_rvs))
rv_mjds = np.linspace(np.min(rv_mjds), np.max(rv_mjds), 256)
# Make a by-eye fit of the radial velocities
# with a sine curve
rv_ys = 43*np.sin(((2*np.pi)/best_period)*(rv_mjds+1.4)) + 9
rv_xs = phase_up(zero, rv_mjds, best_period)
# Put this in order so you can plot it
# as a line
rv_order = np.argsort(rv_xs)
# Plot the by-eye sine curve fit
ax2.plot(rv_xs[rv_order], rv_ys[rv_order],
            '--', c='Blue', zorder=0)
ax2.plot(rv_xs[rv_order]+1, rv_ys[rv_order],
            '--', c='Blue', zorder=0)
# Plot two periods of the SALT and LCO RVs folded to the same
# phase and period as the optical magnitudes
ax2.errorbar(phase_up(zero, SALT_rv_mjds, best_period),
                SALT_rvs, yerr=SALT_rv_errs,
                c='Black', fmt='o', label='SALT RVs')
ax2.errorbar(phase_up(zero, SALT_rv_mjds, best_period)+1,
                SALT_rvs, yerr=SALT_rv_errs,
                c='Black', fmt='o')
ax2.errorbar(phase_up(zero, LCO_rv_mjds, best_period),
                LCO_rvs, yerr=LCO_rv_errs,
                c='Black', fmt='v', label='LCO RVs')
ax2.errorbar(phase_up(zero, LCO_rv_mjds, best_period)+1,
                LCO_rvs, yerr=LCO_rv_errs,
                c='Black', fmt='v')

# Set the y limits so that the plots
# look nice
ax1.set_ylim(11.3, 11.58)
ax2.set_ylim(-60, 60)

# Add the labels and x limit. Make the colour
# of the x-labels match the colours of the
# plots
ax1.set_ylabel('Apparent g magnitude', fontsize=28, color='Grey')
ax2.set_ylabel('Radial velocity (km/s)', fontsize=28)
ax1.set_xlabel('Phase', fontsize=28)
ax1.set_xlim(0, 2)
ax2.set_xlim(0, 2)

# Make a legend
SALT_dot = mlines.Line2D([], [], markerfacecolor='Black', marker='o',
                            markersize=pointsize+4,
                            label='SALT RVs', color='w')
LCO_dot = mlines.Line2D([], [], markerfacecolor='Black', marker='v',
                            markersize=pointsize+4,
                            label='LCO RVs', color='w')
asassn_dot = mlines.Line2D([], [], markerfacecolor='Grey', marker='o',
                            markersize=pointsize+4,
                            label='ASAS-SN g-band'+'\n'+'photometry', color='w')
leg1 = ax2.legend(handles=[SALT_dot, LCO_dot, asassn_dot],
                    fontsize=16, frameon=True, loc='lower left')
# As the first axis has the magnitudes,
# flip the y-axis
ax1.invert_yaxis()

plt.tight_layout()

if save_figures:
    plt.savefig('RV_curve.eps', transparent=False)
plt.close()

swift_isotime = np.array(['2019-04-18T20:45:17', '2019-04-18T22:18:47',
                            '2019-05-05T03:30:48', '2019-05-05T05:11:06',
                            '2019-05-05T09:41:44'])
swift_ap_time = Time(swift_isotime, scale='utc')
swift_mjds = swift_ap_time.mjd

### The MeerKAT observations ###
# We need to the scale the MeerKAT flux densities as the source is
# near the FWHM of the MeerKAT primary beam. This is also where we
# mark the measurements that are non-detections.

meerkat_values = np.load('TYC_MeerKAT_fluxes.npy')
# The times in MJD
thunderkat_times = Time(meerkat_values[:, 0], format='mjd', scale='utc')
thunderkat_mjds = thunderkat_times.mjd
# the peak flux in Jy
thunderkat_flux = meerkat_values[:, 1]
# the uncertainty on the peak flux in Jy
thunderkat_flux_err = meerkat_values[:, 2]
# The RMS noise local to the source
tkt_local_rms = np.load('TYC_MeerKAT_local_RMS.npy')

# The scale value calculated in the paper
# as the source is at the FWHM of the
# primary beam of MeerKAT
tkt_scale = 2.0017550340402122
tkt_scale_err = 0.20242296072260868

# For some points the flux is an extreme outlier.
# This is because of the noise in the image and
# is a artefact. Set the value of any of these points
# to the detection threshold.
outlier = np.where(thunderkat_flux_err>6e-4)[0]
thunderkat_flux_err[outlier] = tkt_local_rms[outlier]
thunderkat_flux[outlier] = 2.*tkt_local_rms[outlier]

# Treat non-detections (below 2 times the local RMS)
# in the same way as outliers
below_noise = np.where((thunderkat_flux/(2*tkt_local_rms))<1)[0]
thunderkat_flux[below_noise] = 2*tkt_local_rms[below_noise]
thunderkat_flux_err[below_noise] = tkt_local_rms[below_noise]

# Plot the outliers and non-detections
# as upper limits
thunderkat_upperlims = np.zeros(len(thunderkat_flux))
thunderkat_upperlims[below_noise] = 1
thunderkat_upperlims[outlier] = 1
# Turn this into a boolean array for
# plotting
uls = thunderkat_upperlims.astype(bool)

# Scale the flux density of the source as it
# is close to the FWHM of the primary beam
# of MeerKAT
tkt_scaled_flux = thunderkat_flux*tkt_scale
tkt_scaled_err = np.abs(tkt_scaled_flux*
                        np.sqrt((thunderkat_flux_err/
                                    thunderkat_flux)**2+
                                (tkt_scale_err/tkt_scale)**2))
scaled_values = np.hstack((np.expand_dims(thunderkat_mjds, axis=1),
                            np.expand_dims(tkt_scaled_flux, axis=1),
                            np.expand_dims(tkt_scaled_err, axis=1)))
np.save('TYC_MeerKAT_ScaledFlux', scaled_values)

### Make Figure 2 ###

fig, ax = plt.subplots(1, 1, figsize=(14, 6))

# Plot the detections as black points
ax.errorbar(thunderkat_mjds[~uls], (tkt_scaled_flux[~uls])*1e3,
            yerr=(np.abs(tkt_scaled_err[~uls]))*1e3, 
            fmt='o', uplims=thunderkat_upperlims[~uls],
            c='Black', label='Detections')
# Plot the upper limits as grey points
ax.errorbar(thunderkat_mjds[uls], (tkt_scaled_flux[uls])*1e3,
            yerr=(np.abs(tkt_scaled_err[uls]))*1e3, 
            fmt='o', uplims=thunderkat_upperlims[uls],
            c='Grey', label='Non-detections')

ax.set_ylabel('Flux density (mJy)', fontsize=20)
ax.set_xlabel('MJD', fontsize=20)

# Make a copy of the axes so that you can
# plot the dates along the top (as well
# as the MJDs along the bottom)
ax2 = ax.twiny()
ax2.xaxis.tick_top()
ax2.set_xticks([58088, 58150, 58209,
                58270, 58331, 58392,
                58453, 58515, 58574,
                58635, 58696, 58757])
ax2.set_xticklabels(['Dec 17','Feb 18','Apr 18',
                        'Jun 18','Aug 18', 'Oct 18',
                        'Dec 18', 'Feb 19', 'Apr 19',
                        'Jun 19', 'Aug 19', 'Oct 19'])

ax.set_ylim(0, 0.85)
ax2.set_ylim(0, 0.85)
ax.set_xlim(58060, 58685)
ax2.set_xlim(58060, 58685)

# Make a nice legend
detection_dot = mlines.Line2D([], [], color='Black', marker='o',
                            markeredgecolor='Black', linestyle='None',
                            markersize=8, markeredgewidth=3,
                            label='Detections')
nondetection_dot = mlines.Line2D([], [], color='Grey', marker='o',
                            markeredgecolor='Grey', linestyle='None',
                            markersize=8, markeredgewidth=3,
                            label='Non-detections')
leg1 = ax.legend(handles=[detection_dot, nondetection_dot],
                    fontsize=18, frameon=True, loc='upper left')

plt.tight_layout()

if save_figures:
    plt.savefig('Radio_LightCurve.eps',
                transparent=False)
plt.close()

### Make Figure 10 ###

pointcol = 'Black'
framecol = 'Black'
fontsize = 22
pointsize = 10

fig, ax = plt.subplots(2, 1, figsize=(19, 11))
j = -2
p = periods[-1]
print(p, optical_name_chunks[j])
zero = optical_mjd_chunks[-2][0]
mjd = optical_mjd_chunks[j]
mag = optical_mag_chunks[j]

ax[0].errorbar(thunderkat_mjds[~uls], (tkt_scaled_flux[~uls])*1e3,
            yerr=(np.abs(tkt_scaled_err[~uls]))*1e3, 
            fmt='o', uplims=thunderkat_upperlims[~uls],
            c='Black', label='Detections')
ax[0].errorbar(thunderkat_mjds[uls], (tkt_scaled_flux[uls])*1e3,
            yerr=(np.abs(tkt_scaled_err[uls]))*1e3, 
            fmt='o', uplims=thunderkat_upperlims[uls],
            c='Grey', label='Non-detections')
j = 16
ax[1].errorbar(optical_mjd_chunks[j], optical_mag_chunks[j],
                yerr=optical_mag_e_chunks[j],
                fmt='v', c=season_colours[j],
                markeredgecolor='DarkGrey')
j = 18
ax[1].errorbar(optical_mjd_chunks[j], optical_mag_chunks[j],
                yerr=optical_mag_e_chunks[j],
                fmt='p', c=season_colours[j],
                markeredgecolor='DarkGrey')
j = 19
ax[1].errorbar(optical_mjd_chunks[j], optical_mag_chunks[j],
                yerr=optical_mag_e_chunks[j],
                fmt='p', c=season_colours[j],
                markeredgecolor='DarkGrey')


for mjd in SALT_rv_mjds:
    SALT_leg = ax[0].plot(np.ones(10)*mjd, np.linspace(-1, 1, 10), '--',
                c='DarkGrey', zorder=0)
    ax[1].plot(np.ones(10)*mjd, np.linspace(10, 12, 10), '--',
                c='DarkGrey', zorder=0)
for mjd in LCO_rv_mjds:
    LCO_leg = ax[0].plot(np.ones(10)*mjd, np.linspace(-1, 1, 10), ':',
                c='Blue', zorder=0, alpha=0.5, label='Epochs of LCO spectra')
    ax[1].plot(np.ones(10)*mjd, np.linspace(10, 12, 10), ':',
                c='Blue', zorder=0, alpha=0.5)

for mjd in swift_mjds:
    swift_leg = ax[0].plot(np.ones(10)*mjd, np.linspace(-1, 1, 10), '-.',
                c='Red', zorder=0, alpha=0.3, label='Epochs of Swift observations')
    ax[1].plot(np.ones(10)*mjd, np.linspace(10, 12, 10), '-.',
                c='Red', zorder=0, alpha=0.3)

parkes_date = Time('2019-01-30T00:00:00.0')
ax[0].plot(np.ones(10)*parkes_date.mjd, np.linspace(-1, 1, 10),
            linestyle=(0, (3, 1, 1, 1, 1, 1)), c='Yellow', zorder=0)
ax[1].plot(np.ones(10)*parkes_date.mjd, np.linspace(10, 12, 10),
            linestyle=(0, (3, 1, 1, 1, 1, 1)), c='Yellow', zorder=0)

parkes_line = mlines.Line2D([], [], color='Yellow', linestyle=(0, (3, 1, 1, 1, 1, 1)),
                            label='Epoch of Parkes\nobservation')
salt_line = mlines.Line2D([], [], color='DarkGrey', linestyle='--',
                            label='Epochs of SALT\nspectra')
lco_line = mlines.Line2D([], [], color='Blue', linestyle=':',
                            label='Epochs of LCOGT\nspectra')
swift_line = mlines.Line2D([], [], color='Red', linestyle='-.',
                            label='Epochs of Swift\nobservations')
TKT_dot = mlines.Line2D([], [], markerfacecolor='Black', marker='o',
                            markersize=pointsize+2,
                            label='ThunderKAT\ndetections', color='w')

TKT_non_dot = mlines.Line2D([], [], markerfacecolor='Grey', marker='o',
                            markersize=pointsize+2,
                            label='ThunderKAT\nnon-detections', color='w')
leg0 = ax[0].legend(handles=[TKT_dot, TKT_non_dot, parkes_line, salt_line, lco_line, swift_line],
                    fontsize=17, frameon=True, loc='upper left')

asassnV_dot = mlines.Line2D([], [], color='None', marker='v',
                            markeredgecolor='DarkGrey',
                            markersize=pointsize, markeredgewidth=3,
                            label='ASAS-SN V band')
asassng_dot = mlines.Line2D([], [], color='None', marker='p',
                            markeredgecolor='DarkGrey',
                            markersize=pointsize, markeredgewidth=3,
                            label='ASAS-SN g band')
leg1 = ax[1].legend(handles=[asassnV_dot, asassng_dot],
                    fontsize=17, frameon=True, loc='lower left')

ax[1].set_xlabel('MJD', fontsize=fontsize, color=framecol)
ax[0].set_ylabel('Flux density (mJy)', fontsize=fontsize, color=framecol)
ax[1].set_ylabel('Apparent magnitude', fontsize=fontsize, color=framecol)

for i in range(2):
    ax[i].set_xlim(58278, 58705)
    ax[i].tick_params(axis='both', colors=framecol)
    ax[i].spines['bottom'].set_color(framecol)
    ax[i].spines['top'].set_color(framecol)
    ax[i].spines['left'].set_color(framecol)
    ax[i].spines['right'].set_color(framecol)

ax[0].set_ylim(-0.02, 0.9)
ax[1].set_ylim(10.75, 11.71)

ax[0].xaxis.tick_top()
ax[0].set_xticks([58331, 58392, 58453,
                    58515, 58574, 58635])
ax[0].set_xticklabels(['Aug 2018', 'Oct 2018',
                        'Dec 2018', 'Feb 2019',
                        'Apr 2019', 'Jun 2019'])

ax[1].invert_yaxis()
plt.tight_layout()
plt.subplots_adjust(hspace=0)

if save_figures:
    plt.savefig('Optical_Radio_LightCurve.eps',
                transparent=False)
plt.close()

### The Gudel Plot ###

# Also known as the Gudel-Benz plot, this is where the X-ray flux of 
# active stars is plotted against the radio flux. I've included the 
# links to the sources of these values in this notebook. The original
# Gudel-Benz plots were shown in Guedel & Benz 1993 and Benz & Guedel
# 1994. This plot (Figure 12 in the paper) was quite difficult to
# reproduce, and the plot shown here and in the paper does not include
# all the points and differet types of object as in the original versions

jy_to_erg = 1e-23 # erg s-1 cm-2 Hz-1
pc_to_cm = 3.086e+18

dMes_radio = np.array(([13.20, 12.74, 13.87, 13.70,
                        13.96, 13.06, 12.86, 12.19,
                        13.34, 13.59, 12.65, 13.59,
                        13.76, 12.15, 13.09, 14.27,
                        14.01, 12.28, 12.63, 13.70,
                        12.57, 13.53]))
dMes_radio_lims = np.array([0, 0, 0, 0,
                            0, 0, 0, 1,
                            0, 0, 1, 0,
                            0, 0, 0, 0,
                            0, 1, 0, 0,
                            1, 0])

dMes_xray = np.array([27.63, 28.18, 28.95, 29.03,
                        29.61, 28.47, 28.78, 27.00,
                        29.21, 28.93, 26.92, 28.64,
                        29.37, 27.75, 29.02, 29.36,
                        29.38, 27.39, 27.73, 29.26,
                        27.12, 28.84])


# Flux at 4.9 GHz in Jy
Flux_4 = np.array([1957, 87, 266, 254,
                    449, 209, 250, 505,
                    np.nan, 135, 1220, np.nan,
                    185, 219, 216, np.nan,
                    86, 2015, 1200, np.nan,
                    np.nan, np.nan, 365, 253,
                    np.nan, 371, 324])*1e-6
Flux_4_err = np.array([40, 38, 53, 31,
                        50, 36, 50, 89,
                        np.nan, 36, 33, np.nan,
                        39, 39, 44, np.nan,
                        36, 100, 600, np.nan,
                        np.nan, np.nan, 33, 33,
                        np.nan, 27, 27])*1e-6

# Flux at 8.5 GHz in Jy
Flux_8 = np.array([1650, 311, 307, 163,
                    268, 323, np.nan, 222,
                    np.nan, 115, 475, np.nan,
                    np.nan, np.nan, 179, 165,
                    87, 1985, 1000, np.nan,
                    82, 82, 360, 245,
                    np.nan, 459, 226])*1e-6
Flux_8_err = np.array([30, 32, 38, 25,
                        43, 28, np.nan, 54,
                        np.nan, 29, 32, np.nan,
                        np.nan, np.nan, 41, 33,
                        29, 100, 400, np.nan,
                        27, 27, 26, 26,
                        np.nan, 21, 21])*1e-6
# distances in pc
distances = np.array([2.7, 4.8, 14.7, 14.2,
                        14.5, 6.0, 4.9, 2.4,
                        2.4, 12.1, 6.2, 6.2,
                        10.5, 10.5, 15.6, 2.9,
                        10.9, 8.8, 8.8, 3.9,
                        4.0, 4.0, 8.3, 8.3,
                        3.6, 6.4, 6.4])

dMes_xray = 10.**np.array([27.63, 28.18, 28.95, 29.03,
                            29.61, 28.47, 28.78, np.nan, 
                            np.nan,
                            29.21, 28.93, np.nan, np.nan, 28.64,
                            29.37, 27.75, 29.02, 29.36,
                            29.38, np.nan, np.nan, 27.73, np.nan, 29.26,
                            np.nan, np.nan, 28.84])

flux_4_weights = 1./(Flux_4_err**2.)
flux_8_weights = 1./(Flux_8_err**2.)

fw4  =  Flux_4 * flux_4_weights
fw8  =  Flux_8 * flux_8_weights

weighted_mean_flux = (fw4+fw8)/(flux_4_weights+flux_8_weights)
weighted_mean_flux_err = np.sqrt(1./(flux_4_weights+flux_8_weights))

weighted_mean_flux[np.isnan(weighted_mean_flux)] = Flux_4[np.isnan(weighted_mean_flux)]
weighted_mean_flux_err[np.isnan(weighted_mean_flux_err)] = Flux_4_err[np.isnan(weighted_mean_flux_err)]

weighted_mean_flux[np.isnan(weighted_mean_flux)] = Flux_8[np.isnan(weighted_mean_flux)]
weighted_mean_flux_err[np.isnan(weighted_mean_flux_err)] = Flux_8_err[np.isnan(weighted_mean_flux_err)]

weighted_mean_erg = (weighted_mean_flux * jy_to_erg) * 4*np.pi*((distances*pc_to_cm)**2)
weighted_mean_erg_err = (weighted_mean_flux_err * jy_to_erg) * 4*np.pi*((distances*pc_to_cm)**2)

# 8.5 GHz flux density in Jy
radioflux_1992AA = np.array([3.02, 0.08, 0.08,
                                5.22, 0.28, 0.29,
                                0.08, 0.12, 0.12,
                                0.21, 0.10, 0.26,
                                2.5])*1e-3
radioflux_err_1992AA = np.array([0.034, 0.027, 0.025,
                                    0.070, 0.044, 0.033,
                                    0.027, 0.040, 0.038,
                                    0.026, 0.033, 0.020,
                                    0.5])*1e-3
radioflux_ul_1992AA = np.array([0, 1, 1,
                                0, 0, 0,
                                1, 0, 1,
                                0, 1, 0,
                                0])

# X-rays in erg s-1
xrays_1992AA = 10.**np.array([29.62, 28.3, 27.27,
                                29.8, 29.6, 29.1,
                                27.6, 29.4, 28.08,
                                30.1, 27.2, 29.7,
                                30.38])
# distances in pc
distances_1992AA = np.array([11.4, 3.3, 4.9,
                                16.7, 13.9, 13.0,
                                4.5, 16.4, 5.1,
                                24, 3.4, 19.4,
                                24])

# radio luminosity in erg s-1 cm-2
radioerg_1992AA = radioflux_1992AA * jy_to_erg * (4.*np.pi*(distances_1992AA*pc_to_cm)**2)
radioerg_err_1992AA = radioflux_err_1992AA * jy_to_erg * (4.*np.pi*(distances_1992AA*pc_to_cm)**2)

xray_1989ApJS = 10.**np.array([31.32, np.nan, 30.21, 31.11, np.nan,
                            31.69, np.nan, 30.38, 31.25, np.nan,
                            31.01, 28.91, 31.25, 30.56, np.nan,
                            29.95, 30.67, 30.55, np.nan, np.nan,
                            np.nan, np.nan, np.nan, 31.28, 29.16,
                            np.nan, np.nan, 30.77, 30.96, 30.72,
                            np.nan, 31.31, 30.80, 30.90, 31.22,
                            np.nan, 31.43, 29.35, 30.48, 30.84,
                            28.77, 32.03, 31.69, 30.89, np.nan,
                            31.38, 30.55, 29.82, 31.29, 30.95,
                            31.32, np.nan, 30.95, 29.87, 32.21,
                            np.nan, 31.08, 31.10, 30.76, 30.71,
                            30.90, np.nan, 30.82, np.nan, 30.94,
                            30.49, 29.94, np.nan, np.nan, np.nan,
                            np.nan, np.nan, 31.17, 31.15, np.nan,
                            np.nan, np.nan, 30.66, 30.35, np.nan,
                            np.nan, 31.06, 30.48, 32.16, 32.01,
                            31.56, np.nan, 31.69, 31.44, 30.81,
                            np.nan, 31.10, np.nan, 30.61, np.nan,
                            31.45, 30.61, np.nan, np.nan, np.nan,
                            np.nan, np.nan, np.nan, np.nan, 30.23,
                            31.32, np.nan, 31.99, 30.05, 28.67,
                            np.nan, np.nan, 31.23, np.nan,
                            29.50, 29.79, 30.55, 31.15, 29.99, 30.07, 29.14, 29.02])

xray_1989ApJS_ul = np.array([1, np.nan, 0, 1, np.nan,
                                1, np.nan, 0, 0, np.nan,
                                1, 0, 0, 0, np.nan,
                                1, 0, 0, np.nan, np.nan,
                                np.nan, np.nan, np.nan, 0, 0,
                                np.nan, np.nan, 0, 0, 0,
                                np.nan, 1, 0, 0, 0,
                                np.nan, 0, 0, 0, 0,
                                0, 0, 0, 0, np.nan,
                                0, 0, 0, 0, 0,
                                0, np.nan, 0, 0, 0,
                                np.nan, 0, 0, 0, 0,
                                0, np.nan, 0, np.nan, 0,
                                0, 0, np.nan, np.nan, np.nan,
                                np.nan, np.nan, 0, 0, np.nan,
                                np.nan, np.nan, 0, 0, np.nan,
                                np.nan, 0, 0, 0, 0,
                                0, np.nan, 0, 0, 0,
                                np.nan, 0, np.nan, 0, np.nan,
                                0, 0, np.nan, np.nan, np.nan,
                                np.nan, np.nan, np.nan, np.nan, 0,
                                0, np.nan, 0, 1, 0,
                                np.nan, np.nan, 0, np.nan, 0,
                                0, 0, 0, 0, 0, 0, 0])



radio_1989ApJS = np.array([0.3, np.mean(np.array([2.72, 3.23])), 1.0, 0.40, 
                            0.36, 0.16, 0.23, 4.97, 0.19, 0.18, 0.21, 0.25,
                            0.89, 0.98, 0.89, 0.17,
                            0.77, np.mean(np.array([8.4, 15.0, 19.6])),
                            1.55, 0.4, 1.56, np.mean(np.array([2.91, 4.89, 0.8])),
                            4.3, np.mean(np.array([2.39, 3.0, 8.0, 11.0, 16.5])),
                            0.15, np.mean(np.array([2.30, 0.67])), 0.37,
                            np.mean(np.array([5.2, 7.8, 17.2, 40.4, 152])),
                            np.mean(np.array([5.0, 7.5, 11.8, 12.0, 180, 21.9])),
                            0.23, np.mean(np.array([0.51, 0.59])), 0.78, 0.70,
                            1.04, 1.08, 0.28, np.mean(np.array([6.5, 37])),
                            0.4, 1.6, 3.6, 0.45, 0.40, np.mean(np.array([2.96, 7.2, 7.0])),
                            0.4, 0.5, np.mean(np.array([0.76, 2.8])), 0.17, 0.15,
                            np.mean(np.array([1.59, 1.15])), 0.4,
                            np.mean(np.array([10.1, 13.7, 20.7, 50.5, 145, 195])),
                            0.79, np.mean(np.array([2.3, 6.4, 9.2, 18.4, 26.7])),
                            0.18, 0.17, 0.24, 3.0, 0.22, 2.24, np.mean(np.array([0.47, 0.56])),
                            0.86, np.mean(np.array([4.67, 0.34, 3.77])), 0.56,
                            np.mean(np.array([1.34, 0.51])), np.mean(np.array([0.26, 0.61])),
                            np.mean(np.array([2.0, 0.8])), 3.65, 0.18, 0.22, 5.0, 0.8, 0.19,
                            np.mean(np.array([2.35, 2.91])), 0.20, 0.40, 0.40, 0.21, 0.24,
                            0.26, 0.17, 0.16, np.mean(np.array([0.62, 0.67, 3.03, 19])),
                            np.mean(np.array([0.61, 2.87, 23.7, 0.84, 2.0])),
                            np.mean(np.array([0.50, 3.59])), np.mean(np.array([1.2, 0.64])),
                            0.27, 0.6, 4.57, 8.65, np.mean(np.array([0.20, 2.7])),
                            np.mean(np.array([2.1, 0.20])), 1.49, np.mean(np.array([17.9, 28.7, 0.87])),
                            0.40, 0.15, 0.14, 0.42, 0.40,
                            0.40, 0.40, 0.80, 0.80, 7.17, 0.32, 0.68, np.mean(np.array([2.5, 18.2, 5.2])),
                            0.7, np.mean(np.array([1.7, 0.7, 1.76, 0.61])), 0.17, 0.22, 0.21, 0.20, 2.02,
                            0.23, 0.16, 0.15, np.mean(np.array([0.18, 0.25])), 0.40, 0.22, 0.25,
                            0.17, 0.28])*1e-3
radio_1989ApJS_ul = np.array([0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1,
                                0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
                                0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1,
                                0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1,
                                1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1,
                                1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 
                                1, 1, 1, 0, 0, 1, 1, 1, 1])

distances_1989ApJS = np.array([100, 65, 50, 85, 88,
                                105, 90, 45, 650, 140,
                                84, 60, 165, 42, 90,
                                80, 30, 21, 63, 110,
                                36, 40, 59, 50, 18,
                                46, 35, 46, 33, 205,
                                35, 250, 55, 22, 250,
                                45, 100, 8, 95, 180,
                                15, 630, 145, 255, 100,
                                205, 290, 50, 460, 350,
                                50, 60, 26, 55, 290,
                                60, 76, 310, 190, 130,
                                300, 310, 325, 360, 220,
                                140, 38, 200, 120, 263,
                                210, 28, 70, 105, 95,
                                90, 110, 110, 54, 240,
                                105, 56, 25, 200, 455, 
                                330, 600, 155, 100, 80,
                                190, 90, 185, 290, 400,
                                310, 75, 340, 135, 450,
                                30, 250, 360, 500, 63,
                                67, 320, 200, 40, 70,
                                195, 250, 160, 42, 155,
                                48, 15, 300, 50, 95,
                                30, 33])

# radio luminosity in erg s-1 cm-2
radioerg_1989ApJS = radio_1989ApJS * jy_to_erg * (4.*np.pi*(distances_1989ApJS*pc_to_cm)**2)

xray_err_1989ApJS = np.zeros(len(xray_1989ApJS_ul))
xray_err_1989ApJS[np.where(xray_1989ApJS_ul==1)[0]] = 0.2*xray_1989ApJS[np.where(xray_1989ApJS_ul==1)[0]]

radio_err_1989ApJS = np.zeros(len(radio_1989ApJS_ul))
radio_err_1989ApJS[np.where(radio_1989ApJS_ul==1)[0]] = 0.2*radio_1989ApJS[np.where(radio_1989ApJS_ul==1)[0]]
radioerg_err_1989ApJS = radio_err_1989ApJS * jy_to_erg * (4.*np.pi*(distances_1989ApJS*pc_to_cm)**2)

meerkat_values = np.load('TYC_MeerKAT_ScaledFlux.npy')
meerkat_flux = meerkat_values[:, 1]
meerkat_flux_err = meerkat_values[:, 2]

meerkat_weights = 1./(meerkat_flux_err**2.)
meerkat_mean_err = np.sqrt(1./np.sum(meerkat_weights))
meerkat_mean = np.sum(meerkat_flux*meerkat_weights)/np.sum(meerkat_weights)

meerkat_erg = meerkat_mean * jy_to_erg * (4.*np.pi*(555*pc_to_cm)**2)
meerkat_erg_err = meerkat_mean_err * jy_to_erg * (4.*np.pi*(555*pc_to_cm)**2)

# These are the X-ra values from Swift
meerkat_xray = 1.72e-12 # erg s-1 cm-2
meerkat_xray = meerkat_xray * (4.*np.pi*(555*pc_to_cm)**2)

meerkat_values = np.load('TYC_MeerKAT_ScaledFlux.npy')
meerkat_flux = meerkat_values[:, 1]
meerkat_flux_err = meerkat_values[:, 2]

meerkat_weights = 1./(meerkat_flux_err**2.)
meerkat_mean_err = np.sqrt(1./np.sum(meerkat_weights))
meerkat_mean = np.sum(meerkat_flux*meerkat_weights)/np.sum(meerkat_weights)

meerkat_erg = meerkat_mean * jy_to_erg * (4.*np.pi*(555*pc_to_cm)**2)
meerkat_erg_err = meerkat_mean_err * jy_to_erg * (4.*np.pi*(555*pc_to_cm)**2)

# These are the X-ra values from Swift
meerkat_xray = 1.72e-12 # erg s-1 cm-2
meerkat_xray = meerkat_xray * (4.*np.pi*(555*pc_to_cm)**2)

### Plot the Spectral Energy Distribution ###
# This SED was made and analysed by Iain McDonald. This is just the plotting
# method for the data points and model. Iain used the software described in
# McDonald et al. (2012, 2017) to fit BT-SETTL model atmospheres to the star.
# The best fit model is plotted here. To find the sources of the data points
# used, please refer to the paper

sed_datapoints = np.genfromtxt('TYC_SED_datapoints.txt',
                                skip_header=15)[:, 3:-1]
sed_model_points = np.genfromtxt('TYC_SED_model.txt')

model_wl = sed_model_points[:, 0] # in angstrom
model_flux = (sed_model_points[:, 1] * model_wl**2.)/1.0e18
model_wl = (model_wl/10.)/1.000275346 # in nm and vacuum wls

fluxes= []
wls = []
for i in range(len(model_flux)//500+1):
    a = i*500
    b = (i+1)*500
    wls.append(np.nanmean(model_wl[a:b]))
    fluxes.append(np.nanmean(model_flux[a:b]))    
wls = np.array(wls)
fluxes = np.array(fluxes)

fig, ax = plt.subplots(1, 1, figsize=(8, 5))

ax.errorbar(sed_datapoints[:, 0], sed_datapoints[:, 3], fmt='o', c='Black')
ax.plot(wls, fluxes, 'Grey', alpha=0.5)
ax.errorbar(2600/10., 3.50, xerr=(500./10.)/2., c='Black', fmt='-')

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim(200, 10000)
ax.set_ylim(0.9, 1500)

ax.set_ylabel((r'Reddening-corrected flux'+'\n'+
                r'$\left(\mathrm{mJy}\right)}=0.51\right)$'), fontsize=20)
ax.set_xlabel(r'Wavelength $\left(\mathrm{nm}\right)$', fontsize=20)

params = (r'$T_{\mathrm{eff}}=5519^{+405}_{-222}\,\mathrm{K}$'+'\n'+
            r'$E \left( B-V \right) =0.51^{+0.08}_{-0.10}\,\mathrm{mag}$'+'\n'+
            r'$L=49.6^{+9.5}_{-8.7}\,\mathrm{L_{\odot}}$')

info_text = AnchoredText(params,
                            loc='lower right',
                            prop=dict(fontsize=20),
                            frameon=False)
ax.add_artist(info_text)

fig.tight_layout()

if save_figures:
    fig.savefig('TYC_SED.eps')
plt.close()