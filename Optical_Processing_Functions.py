import numpy as np
import glob
from astropy.time import Time


def phase_up(ref_mjd, mjds, period):
    '''
    Fold and phase up MJDs to a certain period and MJD

    Takes an array of MJDs and folds it to a period,
    phasing to a reference MJD.

    Args:
    ref_mjds (float): the reference MJD to phase the
                      array of MJDs to
    mjds (array/float): an MJD or array of MJDs to be
                        folded and phased
    period (float): the period to fold to in units
                    of days
    '''
    return (((mjds - ref_mjd)/period)-
            ((mjds - ref_mjd)//period))


def time_chunks(time_array, max_time_gap=30):
    '''
    Cuts a time array when there's a gap bigger than max_time_gap

    Args:
    time_array (array): a numpy array of MJDs
    kwargs:
    max_time_gap (int): the gap between dates where a cut will
                        be made
    '''
    return (np.where((time_array[1:] -
                      time_array[:-1])>30)[0])


def compress_KELT(original_array, kelt_1, kelt_2, kelt_3, names_list):
    '''
    Pops the new, compressed KELT data back into the optical
    data array

    A very specific function, only useful pretty much 
    this one time...
    '''
    names_list = np.array(names_list)

    KELT_start = np.where(names_list=='KELTs East')[0][0]
    KELT_end = np.where(names_list=='KELTs West')[0][-1]
    new_array = original_array[:KELT_start]
    new_array.append(kelt_1)
    new_array.append(kelt_2)
    new_array.append(kelt_3)
    new_array += original_array[KELT_end+1:]
    
    return new_array


def process_ASAS(asasfilename, good_ratings=True, cut_outliers=True):
    '''
    Take the ASAS data file and get the magnitudes and times.
    
    ASAS data has rating A, B, C, D where A data is considered
    the highest quality, and D is considered the lowest quality.
    It also sometimes has issues where there are extremely
    far outliers in magnitude, likely due to instrument
    issues. This function takes the ASAS data and converts it into
    nice, user-friendly numpy arrays.
    
    Args:
    asasfilename (str): the full path and filename of the ASAS
                        tsv file
    kwargs*
    good_rating (bool): Choose whether you would like to use all
                        of the ASAS data regardless of rating, or
                        only the "good" data. If True, only
                        A and B rated data is kept. Otherwise, all
                        A, B, C, and D data is kept.
                        Default: True
    cut_outliers (bool): Choose whether you would like to keep all
                         the data, including extreme outliers. If True,
                         extreme outliers are cut. If False, all data
                         is kept.
                         Default: True
    Returns:
    times (array): numpy array of floats that correspond to the MJD of
                   the magnitudes
    data_points (array): numpy array of floats corresponding to the
                         magnitudes recorded at each MJD
    errs (array): 
    '''
    
    
    asasfile = np.genfromtxt(asasfilename, dtype=str, comments='#')
    
    t = asasfile[:, 0].astype(float)
    asasfile = asasfile[np.argsort(t)]
    
    # The individual ASAS observations are rated
    # with a letter rating for how good the
    # quality of the observation is
    As = np.where(asasfile[:, 11]=='A')[0]
    Bs = np.where(asasfile[:, 11]=='B')[0]
    Cs = np.where(asasfile[:, 11]=='C')[0]
    Ds = np.where(asasfile[:, 11]=='D')[0]

    print('Before processing:')
    print('Length times: ', len(t), '\n')

    if good_ratings:
        # Find which observations are rated A and B
        # i.e. the good quality observations,
        # and remove the rest of the observations
        good_rating = np.sort(np.concatenate((As, Bs)))
        good_obs = asasfile[good_rating]
    else:
        # Use all of the data, regardless
        # of rating
        good_obs = asasfile

    # Remove the letter grading column (column 11) so that
    # the array is now an array of only numbers
    good_obs = np.delete(good_obs, 11, 1)
    good_obs = good_obs.astype(float)

    # For each epoch there are 5 magnitudes
    # for the object. So get those and
    # take the mean (or the measurements
    # and the uncertainties)
    mags = np.mean(good_obs[:, 1:6], axis=1)
    errs = np.mean(good_obs[:, 6:11], axis=1)
    # Get the time of each epoch (in HJD)
    times = good_obs[:, 0]+2450000-2400000.5

    times_original = np.copy(times)
    mags_original = np.copy(mags)
    errs_original = np.copy(errs)
    
    if cut_outliers:
        upper_outliers = np.where(mags<15)[0]
        data_points = mags[upper_outliers]
        times = times[upper_outliers]
        errs = errs[upper_outliers]

        lower_outliers = np.where(data_points>5)[0]
        data_points = data_points[lower_outliers]
        times = times[lower_outliers]
        errs = errs[lower_outliers]
    else:
        data_points = mags

    print('After processing:')
    print('Length times: ', len(times))
    print('Length data points: ', len(data_points))
    print('Length uncertainties: ', len(errs))
    print('Faintest magnitude: {}'.format(np.max(data_points)))
    print('Brightest magnitude: {}'.format(np.min(data_points)))
    
    return times, data_points, errs


def process_ASASSN(asassnfilename):
    '''
    Process the ASAS-SN data into the right format.
    
    Takes the ASAS-SN file in ASAS-SN csv format and
    splits the bands. Also removes extreme outliers
    caused by bad observations.
    
    Args:
    asassnfilename (str): the full path and name of the
                          ASAS-SN csv file (in the format
                          as downloaded from
                          https://asas-sn.osu.edu/)

    Returns:
    Arrays: (asassn_mjd_V, asassn_mag_V, asassn_mag_V_err,
             asassn_mjd_g, asassn_mag_g, asassn_mag_g_err)
            These are the arrays for the V- and g-band MJDs,
            magnitudes, and errors on the magnitudes.
    '''
    asassn_latest = np.genfromtxt(asassnfilename, delimiter=',', dtype=str, skip_header=1)

    # Remove outliers/bad data
    nans = np.where(asassn_latest[:, 5].astype(float)>20)[0]
    asassn = np.delete(asassn_latest, nans, 0)

    # Split the data into V and g bands
    asassn_V = asassn[np.where(asassn[:, -1]=='V')[0]]
    asassn_g = asassn[np.where(asassn[:, -1]=='g')[0]]

    # Convert HJD (Heliocentric Julian Date) to MJD
    asassn_hjd_V = asassn_V[:, 0].astype(float)
    asassn_hjd_g = asassn_g[:, 0].astype(float)
    asassn_mjd_V = asassn_hjd_V-2400000.5
    asassn_mjd_g = asassn_hjd_g-2400000.5

    # Separate out the V mag and V mag error
    asassn_mag_V = asassn_V[:, 5].astype(float)
    asassn_mag_V_err = asassn_V[:, 6].astype(float)

    # Separate out the g mag and g mag error
    asassn_mag_g = asassn_g[:, 5].astype(float)
    asassn_mag_g_err = asassn_g[:, 6].astype(float)
    
    return (asassn_mjd_V, asassn_mag_V, asassn_mag_V_err,
            asassn_mjd_g, asassn_mag_g, asassn_mag_g_err)


def split_semesters(optical_names, optical_mjds,
                    optical_mags, optical_mag_errs,
                    sufficient_points=20):
    '''
    Split the optical observations into semesters.
    
    A semester is the part of the year when an optical
    source is a night-time source (so is visible to
    optical telescopes). Here, we split the optical
    data into semesters, while still maintaining the
    band and telescope information.
    
    Args:
    optical_names (list): a list of strings, where each
                          string is the name given to the
                          observation
    optical_mjds (list): list of arrays, where each array is
                         the array of MJDs corresponding to
                         the name in optical_names
    optical_mags (list): list of arrays, where each array is
                         the array of magnitudes corresponding to
                         the name in optical_names
    optical_mag_errs (list): list of arrays, where each array is
                             the array of errors for the magnitudes
                             corresponding to
                             the name in optical_names
    Returns:
    (optical_mjd_chunks_new,
    optical_mag_chunks_new,
    optical_mag_e_chunks_new,
    optical_name_chunks_new)
    Each of these is a list of arrays. Where each array corresponds to
    the data for a semester for the name in optical_name_chunks_new.
    '''    
    optical_mjd_chunks = []
    optical_mag_chunks = []
    optical_mag_e_chunks = []
    optical_name_chunks = []

    # Split the data into semesters based on
    # the time gap between observations
    for i, obs_mjds in enumerate(optical_mjds):
        tcs = time_chunks(obs_mjds)
        tcs = tcs+1
        tcs = np.insert(tcs, 0, 0)

        for c, chunk in enumerate(tcs[:-1]):
            mjd_chunk = obs_mjds[chunk:tcs[c+1]]

            if len(mjd_chunk)>sufficient_points:
                optical_mjd_chunks.append(mjd_chunk)
                optical_mag_chunks.append(optical_mags[i][chunk:tcs[c+1]])
                optical_mag_e_chunks.append(optical_mag_errs[i][chunk:tcs[c+1]])
                optical_name_chunks.append(optical_names[i])
        if len(obs_mjds[tcs[c+1]:])>sufficient_points:
            optical_mjd_chunks.append(obs_mjds[tcs[c+1]:])
            optical_mag_chunks.append(optical_mags[i][tcs[c+1]:])
            optical_mag_e_chunks.append(optical_mag_errs[i][tcs[c+1]:])
            optical_name_chunks.append(optical_names[i])

    # Create a giant array for all of the V-band magnitudes together
    all_V_mjds = np.concatenate((optical_mjds[0],
                                 optical_mjds[1],
                                 optical_mjds[2],
                                 optical_mjds[3]))
    all_V_mags = np.concatenate((optical_mags[0],
                                 optical_mags[1],
                                 optical_mags[2],
                                 optical_mags[3]))
    all_V_errs = np.concatenate((optical_mag_errs[0],
                                 optical_mag_errs[1],
                                 optical_mag_errs[2],
                                 optical_mag_errs[3]))

    # Add the V-band array to the ther array
    optical_mjd_chunks.append(all_V_mjds)
    optical_mag_chunks.append(all_V_mags)
    optical_mag_e_chunks.append(all_V_errs)
    optical_name_chunks.append('All V data')

    ### Fix up the KELT data so it's one data set
    ### instead of 2.
    # KELT data comes in "East" and "West" arrays for
    # the same night of observations. Here we take the
    # two separate East and West arrays for each semester
    # and combine them.
    kelt_1_mjd = np.concatenate((optical_mjd_chunks[9], optical_mjd_chunks[12]))
    kelt_1_errs = np.concatenate((optical_mag_e_chunks[9], optical_mag_e_chunks[12]))
    kelt_1_mags = np.concatenate((optical_mag_chunks[9],
                                  (optical_mag_chunks[12]-
                                   (-np.mean(optical_mag_chunks[9])+
                                    np.mean(optical_mag_chunks[12])))))
    kelt_2_mjd = np.concatenate((optical_mjd_chunks[10], optical_mjd_chunks[13]))
    kelt_2_errs = np.concatenate((optical_mag_e_chunks[10], optical_mag_e_chunks[13]))
    kelt_2_mags = np.concatenate((optical_mag_chunks[10],
                                  (optical_mag_chunks[13]-
                                   (-np.mean(optical_mag_chunks[10])+
                                    np.mean(optical_mag_chunks[13])))))
    kelt_3_mjd = np.concatenate((optical_mjd_chunks[11], optical_mjd_chunks[14]))
    kelt_3_errs = np.concatenate((optical_mag_e_chunks[11], optical_mag_e_chunks[14]))
    kelt_3_mags = np.concatenate((optical_mag_chunks[11],
                                  (optical_mag_chunks[14]-
                                   (-np.mean(optical_mag_chunks[11])+
                                    np.mean(optical_mag_chunks[14])))))
    
    # Put the new KELT arrays back into the final arrays
    optical_mjd_chunks_new = compress_KELT(optical_mjd_chunks,
                                           kelt_1_mjd, kelt_2_mjd, kelt_3_mjd,
                                           optical_name_chunks)
    optical_mag_chunks_new = compress_KELT(optical_mag_chunks,
                                           kelt_1_mags, kelt_2_mags, kelt_3_mags,
                                           optical_name_chunks)
    optical_mag_e_chunks_new = compress_KELT(optical_mag_e_chunks,
                                             kelt_1_errs, kelt_2_errs, kelt_3_errs,
                                             optical_name_chunks)
    optical_name_chunks_new = compress_KELT(optical_name_chunks,
                                            'KELT', 'KELT', 'KELT',
                                            optical_name_chunks)
    
    return (optical_mjd_chunks_new,
            optical_mag_chunks_new,
            optical_mag_e_chunks_new,
            optical_name_chunks_new)


def bin_per_day(optical_mjd_chunks,
                optical_mag_chunks,
                optical_mag_e_chunks,
                optical_name_chunks):
    '''
    For nicer plotting, bin the optical so that there
    is only one data point per day.
    
    This function takes the mean magnitude and propogated
    error for each day of optical observations. This is
    just so that the plots are clearer and less noisy.
    This also produces the same binned semesters,
    but with extreme outliers removed.
    
    Args:
    optical_mjd_chunks (list): the list of semesters
                               as produced by the
                               split_semesters function
    optical_mag_chunks (list): the list of magnitudes
                               for each semester
                               as produced by the
                               split_semesters function
    optical_mag_e_chunks (list): the list of magnitude errors
                                 for each semester
                                 as produced by the
                                 split_semesters function
    optical_name_chunks (list): the list of names
                                for each semester
                                as produced by the
                                split_semesters function
    Returns:
    (binned_times,
    binned_mags,
    binned_errs,
    no_times,
    no_mags,
    no_errs)
    The binned times, magnitudes, and errors for each semester.
    The same binned times, magnitudes, and errors, but with
    large outliers removed for clarity.
    '''
    binned_times = []
    binned_mags = []
    binned_errs = []

    # NO means outliers are not included
    no_times = []
    no_mags = []
    no_errs = []

    chunks_to_plot = optical_mjd_chunks[:-1]
    for i, times in enumerate(chunks_to_plot):
        mags = np.copy(optical_mag_chunks[i])
        errs = np.copy(optical_mag_e_chunks[i])

        outliers_1 = np.where(mags<(np.mean(mags)-3*np.std(mags)))[0]
        outliers_2 = np.where(mags>(np.mean(mags)+3*np.std(mags)))[0]
        outliers = np.concatenate((outliers_1, outliers_2))

        mags[outliers] = np.nan
        errs[outliers] = np.nan

        no_times.append(times)
        no_mags.append(mags)
        no_errs.append(errs)

        days_of_obs = np.unique(times.astype(int))

        b_times = []
        b_mags = []
        b_errs = []
        for day in days_of_obs:
            same_day = times[np.where(times.astype(int)==day)[0]]
            same_day_mags = mags[np.where(times.astype(int)==day)[0]]
            same_day_errs = errs[np.where(times.astype(int)==day)[0]]

            day_weights = 1/(same_day_errs**2)
            day_mean_mag = np.sum(same_day_mags*day_weights)/np.sum(day_weights)
            day_mean_err = np.sqrt(1/np.sum(day_weights))

            b_times.append(day)
            b_mags.append(day_mean_mag)
            b_errs.append(day_mean_err)

        binned_times.append(b_times)
        binned_mags.append(b_mags)
        binned_errs.append(b_errs)

    return (binned_times,
            binned_mags,
            binned_errs,
            no_times,
            no_mags,
            no_errs)
