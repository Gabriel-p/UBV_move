
import numpy as np
from scipy.ndimage.filters import gaussian_filter


def main(ebv_sig, dm_sig, extin_fix, dst_fix, id_star, x_star, y_star, m_obs,
         bv_obsrv, ub_obsrv, extin_list, dist, hist, xedges, yedges):
    """
    Obtain probable member stars.
    """
    # H_g is the 2D histogram with a Gaussian filter applied, using a large
    # sigma value.
    h_g = gaussian_filter(hist, 7., mode='constant')
    x_max, y_max = np.unravel_index(h_g.argmax(), h_g.shape)
    d_m, ebv_m = np.average(xedges[x_max:x_max + 2]), \
        np.average(yedges[y_max:y_max + 2])

    # Use fixed values if these were set.
    if extin_fix >= 0.:
        E_BV = extin_fix
    else:
        E_BV = ebv_m
    if dst_fix >= 0.:
        dist_kpc = dst_fix
    else:
        dist_kpc = d_m

    # Probable cluster stars.
    id_prob, x_prob, y_prob, m_prob, bv_prob, ub_prob = [], [], [], [], [], []
    # For each observed star.
    for id_st, x_st, y_st, m_star, bv_star, ub_star, ext_star, dist_star in\
        zip(*[id_star, x_star, y_star, m_obs, bv_obsrv, ub_obsrv, extin_list,
              dist]):
        # For each extinction + distance solution assigned to this star.
        for e, d in zip(*[ext_star, dist_star]):
            # Skip dummy '--' values.
            if isinstance(e, float) and isinstance(d, float):
                # If both extinction and distance values are within the
                # accepted ranges, the star is a probable member.
                if abs(e - E_BV) <= ebv_sig and abs(d - dist_kpc) <= dm_sig:
                    if id_st not in id_prob:
                        id_prob.append(id_st)
                        x_prob.append(x_st)
                        y_prob.append(y_st)
                        m_prob.append(m_star)
                        bv_prob.append(bv_star)
                        ub_prob.append(ub_star)

    print("N (probable members) = {}".format(len(id_prob)))

    return E_BV, dist_kpc, id_prob, x_prob, y_prob, m_prob, bv_prob, ub_prob
