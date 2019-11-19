import numpy as np
from astropy.stats import LombScargle
from scipy.interpolate import UnivariateSpline
import Optical_LombScargle_Functions as lsf


def run_LS(optical_mjd_chunks,
           optical_mag_chunks,
           optical_mag_e_chunks,
           optical_name_chunks):
    '''
    Run a Lomb-Scargle on each semester of optical obs.
    
    Args:
    optical_mjd_chunks (list): list of arrays, where each
                               array is an array of MJDs
                               for a semester of observations
    optical_mag_chunks (list): list of arrays, where each
                               array is an array of magnitudes
                               for a semester of observations
    optical_mag_e_chunks (list): list of arrays, where each
                                 array is an array of errors on 
                                 the magnitudes
                                 for a semester of observations
    optical_name_chunks (list): list of arrays, where each
                               array is an array of magnitudes
                               for a semester of observations
    Returns:
    ls_freqs (list)
    ls_pows (list)
    The frequencies and powers from the Lomb-Scargle
    analyses.
    '''
    
    ls_freqs = []
    ls_pows = []

    for i in range(len(optical_name_chunks)):
        if len(optical_mjd_chunks[i])>35:
            (ls_frequency,
             ls_power) = LombScargle(optical_mjd_chunks[i],
                                     optical_mag_chunks[i]).autopower(minimum_frequency=0,
                                                                      maximum_frequency=4)
            ls_freqs.append(ls_frequency)
            ls_pows.append(ls_power)
        else:
            ls_freqs.append([np.nan])
            ls_pows.append([np.nan])

    return ls_freqs, ls_pows


def get_periods(ls_freqs, ls_pows,
                widths_right, widths_left,
                optical_name_chunks):
    '''
    Find the period and the error on the period.
    
    Use the results of the Lomb-Scargle analysis
    to find the period and the error on the period.
    
    Args:
    ls_freqs (list): a list of arrays, where each array
                     gives the frequencies from the LS
                     analysis of a semester of optical
                     observations
    ls_pows (list): the list of arrays of the powers
                    corresponding to the ls_freqs from
                    the LS analysis
    widths_right (list): the approximate positions of
                         the edge of the peak (right-
                         hand side)
    widths_left (list): the approximate positions of
                        the edge of the peak (left-
                        hand side)
    optical_name_chunks (list): list of arrays, where each
                                array is an array of magnitudes
                                for a semester of observations
    Returns:
    The period in days and the error on the period
    '''
    centres = []
    HWHMs = []
    maxpos_periods = []

    for i, freqs in enumerate(ls_freqs):
        if len(freqs)>1:
            # Cut the frequencies and powers so that the arrays
            # aren't huge
            freqs_for_max = freqs[np.where((1/freqs)<100)[0][0]:
                                  np.where((1/freqs)>5)[0][-1]]
            pows_for_max = ls_pows[i][np.where((1/freqs)<100)[0][0]:
                                      np.where((1/freqs)>5)[0][-1]]

            # Get the value and index of the maximum
            # power, and hence the period of the max power
            maxpos = np.nanargmax(pows_for_max)
            maxval = np.nanmax(pows_for_max)
            max_period = 1/freqs_for_max[maxpos]

            # Cut around the max power position to fit for the FWHM
            periods_temp = np.flip(1/freqs_for_max[maxpos-widths_left[i]:
                                                   maxpos+widths_right[i]],
                                   axis=0)
            powers_temp = np.flip((pows_for_max[maxpos-widths_left[i]:
                                                maxpos+widths_right[i]]-
                                   maxval/2.),
                                  axis=0)
            # Find the roots of the FWHM
            spline = UnivariateSpline(periods_temp, powers_temp)
            root1, root2 = spline.roots()
            # The FWHM, and hence, HWHM
            FWHM = np.abs(root1 - root2)
            HWHM = FWHM/2.
            # The middle of the FWHM gives the period value
            centre = root1+HWHM

            period_lower = max_period - HWHM
            period_upper = max_period + HWHM

            print((u'{0} \tperiod: {1:.6f} '
                   '\u00B1 {2:.6f} days \n').format(optical_name_chunks[i],
                                                    centre,
                                                    period_upper-max_period))
            # Put everything into the lists
            centres.append(centre)
            HWHMs.append(HWHM)
            maxpos_periods.append(max_period)
        else:
            # If the LS didn't work, put in nans
            centres.append(np.nan)
            HWHMs.append(np.nan)
            maxpos_periods.append(np.nan)
            print((u'{0} \tperiod: {1:.6f} '
                   '\u00B1 {2:.6f} days \n').format(optical_name_chunks[i],
                                                    np.nan,
                                                    np.nan))
    return centres, HWHMs
