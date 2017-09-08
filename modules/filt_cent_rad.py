
import numpy as np
from scipy.spatial.distance import cdist


def main(in_data, cent, rad):
    """
    """

    ids, x, y = in_data[:3]
    if rad > 0.:
        print("Filter data by center and radius given.")

        dists = cdist(list(zip(*[x, y])), [cent])
        # Store indexes of stars that should be removed.
        del_indexes = [i for i, d in enumerate(dists) if d > rad]

        # Remove stars from id list first since this are strings.
        id_clean = np.delete(np.array(ids), del_indexes)
        # Remove stars from the remaining lists.
        clean = np.delete(np.array(in_data[1:]), del_indexes, axis=1)
    else:
        id_clean, clean = ids, in_data[1:]

    return id_clean, clean
