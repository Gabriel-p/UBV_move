
import numpy as np


def main(ext_dist_all):
    """
    Histogram for distance versus extinction diagram.
    """
    # Maximum distance value in axis.
    d_max = np.sort(ext_dist_all[1])[int(0.9 * len(ext_dist_all[1]))]
    e_max = np.sort(ext_dist_all[0])[int(0.9 * len(ext_dist_all[0]))]

    hist, xedges, yedges = np.histogram2d(
        ext_dist_all[1], ext_dist_all[0],
        bins=[int(d_max / 0.002), int(max(ext_dist_all[0]) / 0.005)])

    return d_max, e_max, hist, xedges, yedges
