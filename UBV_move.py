
import numpy as np
import scipy.spatial.distance as sp
from os import listdir
from os.path import isfile, join
from modules import ext_dist_histo
from modules import prob_membs
from modules import output_file
from modules import make_plots


def get_files():
    '''
    Store the paths and names of all the input clusters stored in the
    input folder.
    '''
    cl_files = [join('input/', f) for f in listdir('input/') if
                isfile(join('input/', f))]

    return cl_files


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


def read_input(data_file):
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


def params_input():
    """
    Read data from file.
    """
    with open('params_input.dat', "r") as f_dat:
        # Iterate through each line in the file.
        for l, line in enumerate(f_dat):
            if not line.startswith("#") and line.strip() != '':
                reader = line.split()
                if reader[0] == 'EM':
                    extin_max = float(reader[1])
                if reader[0] == 'SG':
                    ebv_sig, dm_sig = float(reader[1]), float(reader[2])
                if reader[0] == 'EF':
                    extin_fix = float(reader[1])
                if reader[0] == 'DM':
                    dm_fix = float(reader[1])
                if reader[0] == 'PL':
                    plot_v_seg = True if reader[1] is 'y' else False
                    plot_tcd_chart = True if reader[2] is 'y' else False
                    plot_dens_map = True if reader[3] is 'y' else False

    return extin_max, ebv_sig, dm_sig, extin_fix, dm_fix, plot_v_seg,\
        plot_tcd_chart, plot_dens_map


def intrsc_values(bv_obsrv, ub_obsrv, e_bv):
    '''
    Takes *observed* color and magnitude lists and returns corrected/intrinsic
    lists according to a given extinction law and value.
    '''
    # For UBVI system.
    #
    # E(U-B) = 0.72*E(B-V) +0.05*E(B-V)**2
    #
    # E(B-V) = (B-V) - (B-V)o
    # E(U-B) = (U-B) - (U-B)o
    #
    bv_intrsc = np.array(bv_obsrv) - e_bv
    e_ub = 0.72 * e_bv + 0.05 * e_bv ** 2
    ub_intrsc = np.array(ub_obsrv) - e_ub

    return bv_intrsc, ub_intrsc


def track_distance(zams_inter, bv_intrsc, ub_intrsc):
    '''
    Function that takes as input two arrays of (U-B)o and (B-V)o colors for
    a group of stars and returns for each one the point in the ZAMS closest
    to it.
    '''
    # Create color-color list.
    bv_ub = [bv_intrsc, ub_intrsc]

    # Get distances of *every* star to *every* point in this
    # ZAMS. The list is made of N sub-lists where N is the
    # number of cluster stars and M items in each sub-list where M is
    # the number of points in the ZAMS/track/isochrone.
    # dist_m = [[star_1], [star_2], ..., [star_N]]
    # [star_i] = [dist_1, dist_2, ..., dist_M]
    dist_m = np.array(sp.cdist(list(zip(*bv_ub)),
                               list(zip(*zams_inter[:2])), 'euclidean'))

    # Identify the position of the closest point in track for each star and
    # save the distance to it.
    min_dists = dist_m.min(axis=1)

    # Save index of closest point in ZAMS to each star along with the
    # distance value.
    min_dist_indxs = []
    for indx, dist in enumerate(min_dists):
        # Check if distance star-ZAMS_point is less than this value. This
        # is considered an 'intersection' between the extinction vector and
        # the ZAMS.
        if dist <= 0.01:
            min_dist_indxs.append(dist_m[indx].tolist().index(dist))
        else:
            # If no point in ZAMS was located near this star shifted this
            # value of E(B-V), return '999999'.
            min_dist_indxs.append(999999)

    return min_dist_indxs


def get_track():
    '''
    Read the file with the ZAMS to be used.
    '''
    ZAMS_file = 'modules/zams_SK'
    bv_o, ub_o, M_abs, sp_type = [], [], [], []
    with open(ZAMS_file, mode="r") as z_f:
        for line in z_f:
            li = line.strip()
            # Jump comments.
            if not li.startswith("#"):
                reader = li.split()
                bv_o.append(float(reader[2]))
                ub_o.append(float(reader[1]))
                M_abs.append(float(reader[3]))
                sp_type.append(str(reader[0]))

    return bv_o, ub_o, M_abs, sp_type


def main():
    '''
    Take a cluster data file with photometric (U-B) and (B-V) colors and move
    each star over a given ZAMS thus obtaining its parameters: extinction,
    distance, absolute magnitude and spectral type.
    '''
    extin_max, ebv_sig, dm_sig, extin_fix, dm_fix, plot_v_seg,\
        plot_tcd_chart, plot_dens_map = params_input()

    cl_files = get_files()
    for clust_file in cl_files:

        clust_name = clust_file.split('/')[1].split('.')[0]
        print("\nProcessing: {}".format(clust_name))

        # Define range for E(B-V) value.
        extin_range = np.arange(0, extin_max, 0.01)

        # Get interpolated ZAMS points.
        bv_o, ub_o, M_abs, sp_type = get_track()
        zams = [bv_o, ub_o, M_abs]

        # Interpolate extra colors and absolute magnitudes.
        N = 1000
        bv_col, ub_col = np.linspace(0, 1, len(bv_o)), np.linspace(0, 1, N)
        # One-dimensional linear interpolation.
        bv_col_i, ub_col_i, abs_mag_i = (np.interp(ub_col, bv_col,
                                                   zams[i]) for i in range(3))
        # Store zams' interpolated values.
        zams_inter = np.asarray([bv_col_i, ub_col_i, abs_mag_i])
        print('ZAMS interpolated.')

        # Get input data from file.
        id_star, x_star, y_star, m_obs, e_m, bv_obsrv, e_bv, ub_obsrv, e_ub = \
            read_input(clust_file)

        # List that holds indexes for points in the interpolated ZAMS.
        zams_indxs = [[] for _ in range(len(id_star))]
        zams_old = [[9999.] * len(id_star), [9999.] * len(id_star)]
        # Holds the extinction values assigned to each star.
        extin_list = [[] for _ in range(len(id_star))]

        print('Obtaining extinction solutions.')
        # Loop through all extinction values in range.
        for extin in extin_range:

            # Get intrinsic values for both colors with this extinction value.
            bv_intrsc, ub_intrsc = intrsc_values(bv_obsrv, ub_obsrv, extin)

            # For each star in its intrinsic position find the point in the
            # ZAMS closest to it and store that point's index in the ZAMS and
            # the distance star-point. The point must fall inside a small
            # rectangle around the star.
            min_indxs = track_distance(zams_inter, bv_intrsc, ub_intrsc)

            # For each star, check that an intersection point with the ZAMS was
            # found.
            for indx, n_indx in enumerate(min_indxs):
                if n_indx != 999999:
                    # If the star was positioned over the ZAMS, see that this
                    # new position is significantly different from the previous
                    # one. If it is, store this extinction value as a new
                    # solution for the star. Otherwise the new value is
                    # discarded as an interpolating point in the ZAMS very
                    # close to the previous one.
                    dist = np.sqrt(
                        (zams_old[0][indx] - zams_inter[0][n_indx]) ** 2 +
                        (zams_old[1][indx] - zams_inter[1][n_indx]) ** 2)
                    if dist > 0.1:
                        # Append new value.
                        zams_indxs[indx].append(n_indx)
                        # Update value.
                        zams_old[0][indx] = zams_inter[0][n_indx]
                        zams_old[1][indx] = zams_inter[1][n_indx]
                        # Store extinction value solution.
                        extin_list[indx].append(round(extin, 2))

            print('E(B-V) = %0.2f processed' % extin)

        print('Extinction range processed.')

        # Lists for plotting.
        x_uniq, y_uniq, m_obs_uniq, bv_obs_uniq, ub_obs_uniq, bv_int_uniq,\
            ub_int_uniq = [], [], [], [[], [], [], []], [[], [], [], []], [],\
            []
        # Generate final lists for writing to file.
        M_abs_final = [[] for _ in range(len(id_star))]
        sp_type_final = [[] for _ in range(len(id_star))]
        dist = [[] for _ in range(len(id_star))]
        ext_dist_all = [[], []]

        for indx, star_indxs in enumerate(zams_indxs):

            if star_indxs and len(star_indxs) <= 4:
                for indx2, ind in enumerate(star_indxs):
                    M_abs_final[indx].append(round(zams_inter[2][ind], 3))
                    # Calculate distance.
                    A_v = 3.1 * extin_list[indx][indx2]
                    dist_mod = m_obs[indx] - zams_inter[2][ind]
                    dist[indx].append(
                        round((10 ** (0.2 * (dist_mod + 5 - A_v))) / 1000., 3))
                    # Get spectral type.
                    sp_in = min(list(range(len(M_abs))),
                                key=lambda i: abs(
                                    M_abs[i] - zams_inter[2][ind]))
                    sp_type_final[indx].append(sp_type[sp_in])
                    # Store extin and dist values.
                    ext_dist_all[0].append(extin_list[indx][indx2])
                    ext_dist_all[1].append(
                        round((10 ** (0.2 * (dist_mod + 5 - A_v))) / 1000., 3))

            # Identify stars with a unique solution for plotting.
            if len(star_indxs) == 1:
                # x_uniq.append(x_star[indx])
                # y_uniq.append(y_star[indx])
                # m_obs_uniq.append(m_obs[indx])
                bv_obs_uniq[0].append(bv_obsrv[indx])
                bv_obs_uniq[1].append(e_bv[indx])
                ub_obs_uniq[0].append(ub_obsrv[indx])
                ub_obs_uniq[1].append(e_ub[indx])
                bv_intrsc, ub_intrsc = intrsc_values(
                    bv_obsrv[indx], ub_obsrv[indx], extin_list[indx][0])
                # Corrected values.
                bv_int_uniq.append(bv_intrsc)
                ub_int_uniq.append(ub_intrsc)

            # Fill empty solutions with '--'.
            if len(star_indxs) < 4:
                for _ in range(4 - len(star_indxs)):
                    M_abs_final[indx].append('--')
                    sp_type_final[indx].append('--')
                    dist[indx].append('--')
                    extin_list[indx].append('--')

            if len(star_indxs) > 4:
                print(
                    'Star with too many solutions (>4):\n', int(id_star[indx]),
                    x_star[indx], y_star[indx], extin_list[indx])
                for i in range(4):
                    M_abs_final[indx].append('--')
                    sp_type_final[indx].append('--')
                    dist[indx].append('--')
                    extin_list[indx][i] = 'ERR'

        print('Distances, spectral types and absolute magnitudes obtained.')
        print("N (stars w/ unique solutions) = {}".format(len(bv_int_uniq)))

        d_max, e_max, hist, xedges, yedges = ext_dist_histo.main(ext_dist_all)

        E_BV, dist_mod, x_prob, y_prob, m_prob, bv_prob, ub_prob =\
            prob_membs.main(
                ebv_sig, dm_sig, extin_fix, dm_fix, x_star, y_star, m_obs,
                bv_obsrv, ub_obsrv, extin_list, dist, hist, xedges, yedges)

        output_file.main(
            clust_name, id_star, x_star, y_star, m_obs, e_m, bv_obsrv, e_bv,
            ub_obsrv, e_ub, extin_list, M_abs_final, dist, sp_type_final)
        make_plots.main(
            clust_name, plot_v_seg, plot_tcd_chart, plot_dens_map,
            x_star, y_star, m_obs, bv_o, ub_o, bv_obsrv, ub_obsrv, extin_max,
            ebv_sig, dm_sig, extin_fix, dm_fix, bv_obs_uniq, ub_obs_uniq,
            bv_int_uniq, ub_int_uniq, extin_list, dist, d_max, e_max, hist,
            xedges, yedges, E_BV, dist_mod, x_prob, y_prob, m_prob, bv_prob,
            ub_prob)

    print('End.')


if __name__ == '__main__':
    main()
