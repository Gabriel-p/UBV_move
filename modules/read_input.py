
import numpy as np


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


def main(data_file):
    '''
    Read the file that stores the photometric data for all stars.
    '''
    # Loads the data in 'data_file' as a list of N lists where N is the number
    # of columns. Each of the N lists contains all the data for the column.
    data = np.genfromtxt(data_file, dtype=float, filling_values=99.999,
                         unpack=True)
    id_star, x_star, y_star, m_obs, e_m, bv_obsrv, e_bv, ub_obsrv, e_ub =\
        data[0], data[1], data[2], data[3], data[4], data[5], \
        data[6], data[7], data[8]
    print("N (all stars in file) = {}".format(len(id_star)))

    # If any magnitude or color value (or their errors) is too large, discard
    # that star.
    ids, x, y, V, eV, BV, eBV, UB, eUB = rem_bad_stars(
        id_star, x_star, y_star, m_obs, e_m, bv_obsrv, e_bv, ub_obsrv, e_ub)

    print("N (stars w/ proper data) = {}".format(len(ids)))

    return ids, x, y, V, eV, BV, eBV, UB, eUB
