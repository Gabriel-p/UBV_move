
from os.path import exists
from os import makedirs


def main(clust_name, id_star, x_star, y_star, m_obs, e_m, bv_obsrv, e_bv,
         ub_obsrv, e_ub, extin_list, M_abs_final, dist, sp_type_final):
    """
    Generate output data file.
    """
    # Generate output dir/subdir if it doesn't exist.
    if not exists('output'):
        makedirs('output')

    with open('output/' + clust_name + '_stars_out.dat', "w") as f_out:
        header = [
            '#ID', 'x', 'y', 'V', 'ev', 'BV', 'ebv', 'UB', 'eub', 'E(B-V)',
            'Mv', 'd(kpc)', 'SP', 'E(B-V)', 'Mv', 'd(kpc)', 'SP', 'E(B-V)',
            'Mv', 'd(kpc)', 'SP', 'E(B-V)', 'Mv', 'd(kpc)', 'SP']
        f_out.write(
            '{:<8}{:>11}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}'
            '{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}'
            '{:>10}{:>10}{:>10}{:>10}{:>10}'.format(*header))
        f_out.write('\n')

    lines = [[] for _ in range(len(id_star))]
    for indx, id_st in enumerate(id_star):
        line = [id_star[indx], x_star[indx], y_star[indx], m_obs[indx],
                e_m[indx], bv_obsrv[indx], e_bv[indx], ub_obsrv[indx],
                e_ub[indx]]
        for elem in line:
            lines[indx].append(elem)
        for i in range(4):
            lines[indx].append(extin_list[indx][i])
            lines[indx].append(M_abs_final[indx][i])
            lines[indx].append(dist[indx][i])
            lines[indx].append(sp_type_final[indx][i])

    with open('output/' + clust_name + '_stars_out.dat', "a") as f_out:
        for line in lines:
            f_out.write(
                '{:<9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} '
                '{:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} '
                '{:>9} {:>9} {:>9} {:>9} {:>9}'.format(*line))
            f_out.write('\n')

    print('Output file generated.')
