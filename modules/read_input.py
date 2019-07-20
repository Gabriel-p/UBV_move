
import numpy as np
from astropy.io import ascii


def rem_bad_stars(ids, x, y, V, eV, BV, eBV, UB, eUB):
    '''
    Remove stars from all lists that have too large magnitude or color
    values (or their errors) which indicates a bad photometry.
    '''
    # Set photometric range for accepted stars.
    min_lim, max_lim = -50., 50.

    # Store indexes of stars that should be removed.
    lists_arr = list(zip(*[V, eV, BV, eBV, UB, eUB]))
    del_indexes = [i for i, t in enumerate(lists_arr) if
                   any(e > max_lim for e in t) or any(e < min_lim for e in t)]

    # Remove stars from id list first since this are strings.
    id_clean = np.delete(np.array(ids), del_indexes)
    # Remove stars from the coordinates lists.
    x_clean, y_clean = np.delete(np.array([x, y]), del_indexes, axis=1)
    # Remove stars from the rest of the lists.
    V_clean = np.delete(np.array(V), del_indexes)
    eV_clean = np.delete(np.array(eV), del_indexes)
    BV_clean = np.delete(np.array(BV), del_indexes)
    eBV_clean = np.delete(np.array(eBV), del_indexes)
    UB_clean = np.delete(np.array(UB), del_indexes)
    eUB_clean = np.delete(np.array(eUB), del_indexes)

    return id_clean, x_clean, y_clean, V_clean, eV_clean, BV_clean, eBV_clean,\
        UB_clean, eUB_clean


def main(col_numbers, data_file):
    '''
    Read the file that stores the photometric data for all stars.
    '''
    # Loads the data in 'data_file'.
    data = ascii.read(
        data_file,
        fill_values=[
            ('', '0'), ('--', '0'), ('INDEF', '0'), ('nan', '0')]).filled(99.9)
    # Separate into the appropriate types of data given by 'col_numbers'.
    all_data = []
    for cn in col_numbers:
        all_data.append(data.columns[cn].tolist())
    id_star, x_star, y_star, m_obs, e_m, bv_obsrv, e_bv, ub_obsrv, e_ub =\
        all_data
    print("N (all stars in file) = {}".format(len(id_star)))

    # If any magnitude or color value (or their errors) is too large, discard
    # that star.
    ids, x, y, V, eV, BV, eBV, UB, eUB = rem_bad_stars(
        id_star, x_star, y_star, m_obs, e_m, bv_obsrv, e_bv, ub_obsrv, e_ub)

    print("N (stars w/ proper data) = {}".format(len(ids)))

    return ids, x, y, V, eV, BV, eBV, UB, eUB
