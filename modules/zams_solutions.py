
from .ext_solutions import intrsc_values


def main(id_star, x_star, y_star, extin_list, zams_indxs, zams_inter, M_abs,
         sp_type, m_obs, bv_obsrv, e_bv, ub_obsrv, e_ub):
    """
    For each solution assigned to each observed star, find its absolute
    magnitude, intrinsic colors, spectral types, and distances (in kpc).
    """

    # Generate final lists for writing to file.
    M_abs_final = [[] for _ in range(len(id_star))]
    bv_final = [[] for _ in range(len(id_star))]
    ub_final = [[] for _ in range(len(id_star))]
    sp_type_final = [[] for _ in range(len(id_star))]
    dist = [[] for _ in range(len(id_star))]

    # Store *all* extinction and distance values. Used to generate the
    # extinction-distance density maps.
    ext_dist_all = [[], []]

    # Store unique solutions.
    x_uniq, y_uniq, m_uniq, d_uniq = [], [], [], []
    id_uniq, bv_obs_uniq, ub_obs_uniq, bv_int_uniq, ub_int_uniq =\
        [], [[], [], [], []], [[], [], [], []], [], []

    for indx, star_indxs in enumerate(zams_indxs):

        if star_indxs and len(star_indxs) <= 4:
            for indx2, ind in enumerate(star_indxs):
                # Store absolute magnitudes.
                M_abs_final[indx].append(round(zams_inter[2][ind], 3))
                # Store intrinsic colors
                bv_final[indx].append(round(zams_inter[0][ind], 3))
                ub_final[indx].append(round(zams_inter[1][ind], 3))
                # Get spectral type.
                sp_in = min(list(range(len(M_abs))),
                            key=lambda i: abs(
                                M_abs[i] - zams_inter[2][ind]))
                sp_type_final[indx].append(sp_type[sp_in])
                # Calculate distance.
                E_BV = extin_list[indx][indx2]
                A_v = 3.1 * E_BV
                dist_mod = m_obs[indx] - zams_inter[2][ind]
                d_kpc = round((10 ** (0.2 * (dist_mod + 5 - A_v))) / 1000., 3)
                dist[indx].append(d_kpc)

                # Store all extinction and dist values.
                ext_dist_all[0].append(E_BV)
                ext_dist_all[1].append(d_kpc)

        # Identify stars with a unique solution for plotting.
        if len(star_indxs) == 1:
            id_uniq.append(id_star[indx])
            x_uniq.append(x_star[indx])
            y_uniq.append(y_star[indx])
            m_uniq.append(m_obs[indx])
            bv_obs_uniq[0].append(bv_obsrv[indx])
            bv_obs_uniq[1].append(e_bv[indx])
            ub_obs_uniq[0].append(ub_obsrv[indx])
            ub_obs_uniq[1].append(e_ub[indx])
            # Distances.
            E_BV = extin_list[indx][0]
            A_v = 3.1 * E_BV
            dist_mod = m_obs[indx] - zams_inter[2][star_indxs[0]]
            d_kpc = round((10 ** (0.2 * (dist_mod + 5 - A_v))) / 1000., 3)
            d_uniq.append(d_kpc)
            # Corrected values.
            bv_intrsc, ub_intrsc = intrsc_values(
                bv_obsrv[indx], ub_obsrv[indx], extin_list[indx][0])
            bv_int_uniq.append(bv_intrsc)
            ub_int_uniq.append(ub_intrsc)

        # Fill empty solutions with '--'.
        if len(star_indxs) < 4:
            for _ in range(4 - len(star_indxs)):
                M_abs_final[indx].append('--')
                bv_final[indx].append('--')
                ub_final[indx].append('--')
                sp_type_final[indx].append('--')
                dist[indx].append('--')
                extin_list[indx].append('--')

        if len(star_indxs) > 4:
            print(
                'Star with too many solutions (>4):\n', int(id_star[indx]),
                x_star[indx], y_star[indx], extin_list[indx])
            for i in range(4):
                M_abs_final[indx].append('--')
                bv_final[indx].append('--')
                ub_final[indx].append('--')
                sp_type_final[indx].append('--')
                dist[indx].append('--')
                extin_list[indx][i] = 'ERR'

    print('Distances, spectral types and intrinsic mags/colors obtained.')
    print("N (stars w/ unique solutions) = {}".format(len(bv_int_uniq)))

    return extin_list, ext_dist_all, M_abs_final, bv_final, ub_final, dist,\
        sp_type_final, id_uniq, x_uniq, y_uniq, m_uniq, d_uniq, bv_obs_uniq,\
        ub_obs_uniq, bv_int_uniq, ub_int_uniq
