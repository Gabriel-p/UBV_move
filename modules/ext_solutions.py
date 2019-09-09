
import numpy as np
import scipy.spatial.distance as sp


def main(extin_max, N_st, bv_obsrv, ub_obsrv, zams_inter):
    """
    Find valid E(B-V) solutions for observed stars.
    """
    # Initialize dummy large points in the ZAMS.
    zams_old = [[9999.] * N_st, [9999.] * N_st]
    # Holds the extinction values assigned to each star.
    extin_list = [[] for _ in range(N_st)]
    # List that holds indexes for points in the interpolated ZAMS.
    zams_indxs = [[] for _ in range(N_st)]

    print('Obtaining extinction solutions.')
    # Define range for E(B-V) solutions.
    extin_range = np.arange(0, extin_max, 0.01)
    # Loop through all extinction values in range.
    for extin in extin_range:

        # Get intrinsic values for both colors with this extinction value.
        bv_intrsc, ub_intrsc = intrsc_values(bv_obsrv, ub_obsrv, extin)

        # For each star in its intrinsic position find the point in the
        # ZAMS closest to it and store that point's index in the ZAMS and
        # the distance star-point. The point must fall inside a small
        # rectangle around the star.
        min_indxs = track_distance(zams_inter, bv_intrsc, ub_intrsc)

        # Check if any solution was found. Otherwise exit.
        if not np.isnan(min_indxs).all():

            # For each star, check that an intersection point with the ZAMS was
            # found.
            for indx, n_indx in enumerate(min_indxs):
                # If the star was positioned over the ZAMS, see that this
                # new position is significantly different from the previous
                # one. If it is, store this extinction value as a new
                # solution for the star. Otherwise the new value is
                # discarded as an interpolating point in the ZAMS, very
                # close to the previous one.
                if not np.isnan(n_indx):
                    dist = np.sqrt(
                        (zams_old[0][indx] - zams_inter[0][int(n_indx)]) ** 2 +
                        (zams_old[1][indx] - zams_inter[1][int(n_indx)]) ** 2)
                    if dist > 0.1:
                        # Append new value.
                        zams_indxs[indx].append(int(n_indx))
                        # Update value.
                        zams_old[0][indx] = zams_inter[0][int(n_indx)]
                        zams_old[1][indx] = zams_inter[1][int(n_indx)]
                        # Store extinction value solution.
                        extin_list[indx].append(round(extin, 2))
        else:
            print("No solutions found beyond E(B-V)={:.2f}".format(extin))
            break

    print('Extinction range processed.')
    return extin_list, zams_indxs, extin


def intrsc_values(bv_obsrv, ub_obsrv, e_bv):
    '''
    Takes *observed* color and magnitude lists and returns corrected/intrinsic
    lists according to a given extinction law and value.
    '''
    # For UBVI system.
    #
    # E(U-B) = 0.72*E(B-V) + 0.05*E(B-V)**2
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
    bv_ub = np.array([bv_intrsc, ub_intrsc]).T

    # Get distances of *every* star to *every* point in this
    # ZAMS. The list is made of N sub-lists where N is the
    # number of cluster stars and M items in each sub-list where M is
    # the number of points in the ZAMS/track/isochrone.
    # dist_m = [[star_1], [star_2], ..., [star_N]]
    # [star_i] = [dist_1, dist_2, ..., dist_M]
    dist_m = sp.cdist(bv_ub, zams_inter[:2].T, 'euclidean')

    # Identify the distances of the closest point in the track for each star.
    min_dists = dist_m.min(axis=1)
    # Identify the IDs of the same closest points.
    min_dist_indxs = np.argmin(dist_m, axis=1)

    # Mask those stars with star-ZAMS_point distance greater than this value.
    # This is considered an 'intersection' between the extinction vector and
    # the ZAMS.
    msk = min_dists >= 0.01

    # Change to float type, otherwise the mask can not be applied.
    min_dist_indxs = min_dist_indxs.astype(float)
    # Assign nan values to all masked points.
    min_dist_indxs[msk] = np.nan

    return min_dist_indxs
