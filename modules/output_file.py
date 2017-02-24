
from os.path import exists
from os import makedirs


def main(clust_name, id_star, x_star, y_star, m_obs, e_m, bv_obsrv, e_bv,
         ub_obsrv, e_ub, extin_list, M_abs_final, bv_final, ub_final, dist,
         sp_type_final, id_uniq, id_prob, sols_write):
    """
    Generate output data file.
    """
    # Generate output dir/subdir if it doesn't exist.
    if not exists('output'):
        makedirs('output')

    if sols_write == 'unique':
        with open('output/' + clust_name + '_stars_out.dat', "w") as f_out:
            header = [
                '#ID', 'x', 'y', 'V', 'ev', 'BV', 'ebv', 'UB', 'eub',
                'Mv', 'BVo', 'UBo', 'E(B-V)', 'd(kpc)', 'SP']
            f_out.write(
                '{:<8}{:>11}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}'
                '{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}'.format(*header))
            f_out.write('\n')
    else:
        with open('output/' + clust_name + '_stars_out.dat', "w") as f_out:
            header = [
                '#ID', 'x', 'y', 'V', 'ev', 'BV', 'ebv', 'UB', 'eub',
                'Mv', 'BVo', 'UBo', 'E(B-V)', 'd(kpc)', 'SP',
                'Mv', 'BVo', 'UBo', 'E(B-V)', 'd(kpc)', 'SP',
                'Mv', 'BVo', 'UBo', 'E(B-V)', 'd(kpc)', 'SP',
                'Mv', 'BVo', 'UBo', 'E(B-V)', 'd(kpc)', 'SP']
            f_out.write(
                '{:<8}{:>11}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}'
                '{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}'
                '{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}'
                '{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}'
                '{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}'.format(*header))
            f_out.write('\n')

    lines = []
    for indx, id_st in enumerate(id_star):
        line = [id_star[indx], x_star[indx], y_star[indx], m_obs[indx],
                e_m[indx], bv_obsrv[indx], e_bv[indx], ub_obsrv[indx],
                e_ub[indx]]
        # See if only unique stars should be printed to file.
        write_flag = True
        if sols_write == 'unique':
            if id_st not in id_uniq:
                write_flag = False
        if sols_write == 'probable':
            if id_st not in id_prob:
                write_flag = False
        if write_flag:
            temp = []
            for elem in line:
                temp.append(elem)
            for i in range(4):
                temp.append(M_abs_final[indx][i])
                temp.append(bv_final[indx][i])
                temp.append(ub_final[indx][i])
                temp.append(extin_list[indx][i])
                temp.append(dist[indx][i])
                temp.append(sp_type_final[indx][i])
            lines.append(temp)

    with open('output/' + clust_name + '_stars_out.dat', "a") as f_out:
        for line in lines:
            if sols_write == 'unique':
                f_out.write(
                    '{:<9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} '
                    '{:>9} {:>9} {:>9} {:>9} {:>9} {:>9}'.format(*line))
            else:
                f_out.write(
                    '{:<9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} '
                    '{:>9} {:>9} {:>9} {:>9} {:>9} {:>9} '
                    '{:>9} {:>9} {:>9} {:>9} {:>9} {:>9} '
                    '{:>9} {:>9} {:>9} {:>9} {:>9} {:>9} '
                    '{:>9} {:>9} {:>9} {:>9} {:>9} {:>9}'.format(*line))
            f_out.write('\n')

    print('Output file generated.')
