import numpy as np
from astropy.stats import LombScargle
from scipy.interpolate import UnivariateSpline
import Optical_LombScargle_Functions as lsf


###########################################################
### Perform a Lomb-Scargle analysis of optical observations of TYC 8332-2529-1

# This python script performs Lomb-Scargle (LS) analyses
# of each semester of optical observations
# of the source TYC 8332-2529-1 and writes two numpy files.

# The file TYC_LS_periods.npy has a numpy array of
# the periods (in days) from each semester in it.

# The file TYC_LS_periodErrors.npy has a numpy array
# of the errors on the periods (in days) from each semester in it.

if __name__ in '__main__':
    # First load in the optical data that was made using
    # Process_Optical_Data.ipynb (or Process_Optical_Data.py).
    all_optical_info = np.load('TYC_optical_semesters.npy')

    optical_mjd_chunks = all_optical_info[0]
    optical_mag_chunks = all_optical_info[1]
    optical_mag_e_chunks = all_optical_info[2]
    optical_name_chunks = all_optical_info[3]
    mean_all_V = np.mean(optical_mag_chunks[-1])

    # Now perform the Lomb-Scargle analysis on each
    # semester (and all of the V-band observations).
    print('\n---------------------------------------------'
          '\nThe period and uncertainty for each semester:')
    (ls_freqs, ls_pows) = lsf.run_LS(optical_mjd_chunks,
                                     optical_mag_chunks,
                                     optical_mag_e_chunks,
                                     optical_name_chunks)

    # Next find the peak of the power spectrum to get the
    # period and the uncertainties (in units of days), and
    # save these values in numpy files.

    # The approximate positions of the base
    # of the LS peaks for finding the centre
    # of the peak and the FWHM
    widths_right = [7, 5, 5, 6, 6, 6, 5, 7, 8,
                    7, 6, 4, 5, 6, 6, 6, 5, 5,
                    5, 4, 5, 10]
    widths_left  = [7, 4, 5, 5, 6, 5, 5, 5, 5,
                    7, 4, 4, 5, 6, 6, 6, 4, 5,
                    5, 4, 5, 10]

    periods, period_errors = lsf.get_periods(ls_freqs, ls_pows,
                                             widths_right, widths_left,
                                             optical_name_chunks)

    # Save the periods and the uncertainties
    # on the periods
    np.save('TYC_LS_periods', np.array(periods))
    np.save('TYC_LS_periodErrors', np.array(period_errors))
