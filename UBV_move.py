
from os import listdir
from os.path import isfile, join
from modules import zams_interp
from modules import read_input
from modules import filt_cent_rad
from modules import ext_solutions
from modules import zams_solutions
from modules import ext_dist_histo
from modules import prob_membs
from modules import output_file
from modules import make_plots
from modules._version import __version__


def params_input():
    """
    Read data from file.
    """
    with open('params_input.dat', "r") as f_dat:
        # Iterate through each line in the file.
        for l, line in enumerate(f_dat):
            if not line.startswith("#") and line.strip() != '':
                reader = line.split()
                if reader[0] == 'CN':
                    col_numbers = list(map(int, reader[1:]))
                if reader[0] == 'FL':
                    cx, cy, r = list(map(float, reader[1:-1]))
                    cent, rad = [cx, cy], r
                    inout_flag = reader[-1]
                if reader[0] == 'EM':
                    extin_max = float(reader[1])
                if reader[0] == 'SG':
                    ebv_sig, dm_sig = float(reader[1]), float(reader[2])
                if reader[0] == 'FI':
                    extin_fix, dst_fix = float(reader[1]), float(reader[2])
                if reader[0] == 'SO':
                    sols_write = reader[1]
                if reader[0] == 'PL':
                    plot_v_seg = True if reader[1] is 'y' else False
                    plot_dens_map = True if reader[2] is 'y' else False
                    dens_flag = reader[3]
                    plot_tcd_uniq = True if reader[4] is 'y' else False
                    plot_tcd_prob = True if reader[5] is 'y' else False

    if extin_max < extin_fix:
        print("\nMaximum E(B-V) value can not be smaller than "
              "the fixed E(B-V) value.\n")
        raise SystemExit

    return col_numbers, cent, rad, inout_flag, extin_max, ebv_sig, dm_sig,\
        extin_fix, dst_fix, sols_write, plot_v_seg, plot_dens_map, dens_flag,\
        plot_tcd_uniq, plot_tcd_prob


def get_files():
    '''
    Store the paths and names of all the input clusters stored in the
    input folder.
    '''
    cl_files = [join('input/', f) for f in listdir('input/') if
                isfile(join('input/', f))]

    return cl_files


def main():
    '''
    Take a cluster data file with photometric (U-B) and (B-V) colors and move
    each star over a given ZAMS thus obtaining its parameters: extinction,
    distance, absolute magnitude and spectral type.
    '''

    print('\n-------------------------------------------')
    print('             [UBV-move {}]'.format(__version__))
    print('-------------------------------------------\n')

    col_numbers, cent, rad, inout_flag, extin_max, ebv_sig, dm_sig,\
        extin_fix, dst_fix, sols_write, plot_v_seg, plot_dens_map, dens_flag,\
        plot_tcd_uniq, plot_tcd_prob = params_input()

    zams_inter, bv_o, ub_o, M_abs, sp_type = zams_interp.main()

    # Process all files inside 'input/' folder.
    cl_files = get_files()
    for clust_file in cl_files:

        clust_name = clust_file.split('/')[1].split('.')[0]
        print("\nProcessing: {}".format(clust_name))

        # Get input data from file.
        in_data = read_input.main(col_numbers, clust_file)

        # Filter by center and radius.
        id_star, data_clean = filt_cent_rad.main(
            in_data, cent, rad, inout_flag)
        x_star, y_star, m_obs, e_m, bv_obsrv, e_bv, ub_obsrv, e_ub = data_clean

        # Resolve stars on ZAMS.
        extin_list, zams_indxs, max_extin = ext_solutions.main(
            extin_max, len(id_star), bv_obsrv, ub_obsrv, zams_inter)

        # Assign intrinsic mag/colors, distance, spectral types for each
        # solution; identify stars with  unique solutions.
        extin_list, ext_dist_all, M_abs_final, bv_final, ub_final, dist,\
            sp_type_final, id_uniq, x_uniq, y_uniq, m_uniq, d_uniq, ext_unq,\
            bv_obs_uniq, ub_obs_uniq, bv_int_uniq, ub_int_uniq =\
            zams_solutions.main(
                id_star, x_star, y_star, extin_list, zams_indxs,
                zams_inter, M_abs, sp_type, m_obs, bv_obsrv, e_bv,
                ub_obsrv, e_ub)

        # Generate 2D histogram of found valid distances and extinctions.
        # Find the maximum density value, identified as the most probable
        # distance and extinction.
        if dens_flag == 'all':
            d_max, e_max, hist, xedges, yedges = ext_dist_histo.main(
                ext_dist_all)
        else:
            d_max, e_max, hist, xedges, yedges = ext_dist_histo.main(
                [ext_unq, d_uniq])

        # Identify most probable members using either the extinction and
        # distance values found above, or fixed ones given by the user.
        E_BV, dist_kpc, id_prob, x_prob, y_prob, m_prob, bv_prob, ub_prob =\
            prob_membs.main(
                ebv_sig, dm_sig, extin_fix, dst_fix, id_star, x_star, y_star,
                m_obs, bv_obsrv, ub_obsrv, extin_list, dist, hist, xedges,
                yedges)

        output_file.main(
            clust_name, id_star, x_star, y_star, m_obs, e_m, bv_obsrv, e_bv,
            ub_obsrv, e_ub, extin_list, M_abs_final, bv_final, ub_final, dist,
            sp_type_final, id_uniq, id_prob, sols_write)

        make_plots.main(
            clust_name, plot_v_seg, plot_dens_map, plot_tcd_uniq,
            plot_tcd_prob, x_star, y_star, m_obs, bv_o, ub_o, bv_obsrv,
            ub_obsrv, max_extin, ebv_sig, dm_sig, x_uniq, y_uniq, m_uniq,
            d_uniq, bv_obs_uniq, ub_obs_uniq, bv_int_uniq, ub_int_uniq, dist,
            d_max, e_max, hist, xedges, yedges, E_BV, dist_kpc, x_prob,
            y_prob, m_prob, bv_prob, ub_prob)

    print('End.')


if __name__ == '__main__':
    main()
