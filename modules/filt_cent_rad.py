
import numpy as np
from scipy.spatial.distance import cdist


def main(in_data, cent, rad, inout_flag):
    """
    Filter stars by center and radius given.
    """

    ids, x, y = in_data[:3]
    if rad > 0.:
        print("Filter data by center and radius given.")

        dists = cdist(list(zip(*[x, y])), [cent])
        # Store indexes of stars that should be removed.
        if inout_flag == 'in':
            del_indexes = [i for i, d in enumerate(dists) if d > rad]
        else:
            del_indexes = [i for i, d in enumerate(dists) if d <= rad]

        # Remove stars from id list first since this are strings.
        id_clean = np.delete(np.array(ids), del_indexes)
        # Remove stars from the remaining lists.
        clean = np.delete(np.array(in_data[1:]), del_indexes, axis=1)
        print("Stars left after filtering: {}".format(len(id_clean)))
    else:
        id_clean, clean = ids, in_data[1:]

    return id_clean, clean
