import numpy as np
import glob
from astropy.time import Time
import Optical_Processing_Functions as opf


### Process the optical observations of TYC 8332-2529-1

# This script processes the optical observations and write
# the data as numpy files

# The file TYC_optical_semesters.npy has the name,
# MJDs, magnitudes, and errors on the magnitudes of the
# ASAS, KELT, and ASAS-SN observations.

# The file TYC_optical_binned.npy has the same
# semesters, but binned such that there is only one MJD,
# magntitude, and error per day (for each day that the
# source is observed).

# The file TYC_optical_binned_noOutliers.npy is
# the same as TYC_optical_binned.npy, but the
# outlier magnitudes have been removed.

if __name__ in '__main__':
    # First load in the KELT data and convert from
    # HJD to MJD
    east_tfa = np.genfromtxt('KELT_S36_lc_027057_V01_east_tfa.dat')
    west_tfa = np.genfromtxt('KELT_S36_lc_027056_V01_west_tfa.dat')

    east_tfa_time = east_tfa[:, 0]-2400000.5
    east_tfa_V = east_tfa[:, 1] - 2.07
    east_tfa_V_err = east_tfa[:, 2]

    west_tfa_time = west_tfa[:, 0]-2400000.5
    west_tfa_V = west_tfa[:, 1] - 2.07
    west_tfa_V_err = west_tfa[:, 2]

    # Now upload the ASAS data and process it to
    # remove the low-quality epochs and outliers
    print('Processing the ASAS observations')
    asasfilename = 'ASAS_data.tsv'

    (ASAS_times,
     ASAS_mags,
     ASAS_mag_errs) = opf.process_ASAS(asasfilename,
                                       good_ratings=True,
                                       cut_outliers=True)

    # Now upload the ASAS-SN data and split into
    # g-band and V-band
    asassnfilename = 'ASASSN.csv'
    (asassn_mjd_V,
     asassn_mag_V,
     asassn_mag_V_err,
     asassn_mjd_g,
     asassn_mag_g,
     asassn_mag_g_err) = opf.process_ASASSN(asassnfilename)

    # Put all of the optical observations together
    optical_names = ['ASAS', 'KELTs East', 'KELTs West',
                     'ASAS-SN V', 'ASAS-SN g']
    optical_mjds = [ASAS_times, east_tfa_time, west_tfa_time,
                    asassn_mjd_V, asassn_mjd_g]
    optical_mags = [ASAS_mags, east_tfa_V, west_tfa_V,
                    asassn_mag_V, asassn_mag_g]
    optical_mag_errs = [ASAS_mag_errs, east_tfa_V_err, west_tfa_V_err,
                        asassn_mag_V_err, asassn_mag_g_err]

    # Split the observations into semesters
    (optical_mjd_chunks,
     optical_mag_chunks,
     optical_mag_e_chunks,
     optical_name_chunks) = opf.split_semesters(optical_names, optical_mjds,
                                                optical_mags, optical_mag_errs,
                                                sufficient_points=20)
    # Save these as a numpy file
    np.save('TYC_optical_semesters', np.array([optical_mjd_chunks,
                                                optical_mag_chunks,
                                                optical_mag_e_chunks,
                                                optical_name_chunks]))

    # Bin the data so that there is only one data point
    # per day that the source in observed
    (binned_times,
     binned_mags,
     binned_errs,
     no_times,
     no_mags,
     no_errs) = opf.bin_per_day(optical_mjd_chunks,
                                optical_mag_chunks,
                                optical_mag_e_chunks,
                                optical_name_chunks)
    # Save these as numpy files
    np.save('TYC_optical_binned', np.array([binned_times,
                                            binned_mags,
                                            binned_errs]))
    np.save('TYC_optical_binned_noOutliers', np.array([no_times,
                                                       no_mags,
                                                       no_errs]))

    ############################################################
    # Print out some information for each data set so that
    # you have an idea of what you're working with
    print('\n------------------------------------'
          '\nInformation for each data set')
    total_days = 0
    num_points = 0
    for i, name in enumerate(optical_names):
        print('{0} Start MJD: {1:.2f} End MJD: {2:.2f}'.format(name,
                                                               np.min(optical_mjds[i]),
                                                               np.max(optical_mjds[i])))
        print('Number of data points: {}'.format(len(optical_mjds[i])))
        print('Start ISOT: ', Time(np.min(optical_mjds[i]), format='mjd').isot)
        print('End ISOT', Time(np.max(optical_mjds[i]), format='mjd').isot)
        days = (Time(np.max(optical_mjds[i]),
                     format='mjd')-
                Time(np.min(optical_mjds[i]),
                     format='mjd'))
        print('Total number of days: {0}'.format(int((days.sec)/60./60./24.)))
        print()

        num_points += len(optical_mjds[i])
        total_days += (days.sec)/60./60./24.

    print('Total number of days of observations: {}'.format(int(total_days)))
    print('Total number of data points: {}'.format(num_points))
    print('Total span of observations: {:.2f} years'.format((Time(np.max(optical_mjds[-1]),
                                                              format='mjd')-
                                                         Time(np.min(optical_mjds[0]),
                                                              format='mjd')).sec/60./60./24./365.25))
