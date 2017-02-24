
import numpy as np


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
    """
    """
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

    return zams_inter, bv_o, ub_o, M_abs, sp_type
