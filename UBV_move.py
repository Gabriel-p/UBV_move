# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 17:02:23 2013

@author: gabriel
"""


'''
Take a star data file with photometric (U-B) and (B-V) colors and place each
one over a given ZAMS obtaining the parameters for all.
'''

import numpy as np
import scipy.spatial.distance as sp
import matplotlib.pyplot as plt


def intrsc_values(bv_obsrv, ub_obsrv, e_bv):
    '''
    Takes *observed* color and magnitude lists and returns corrected or
    intrinsic lists. Depends on the system selected/used.
    '''
    # For UBVI system.
    #
    # E(U-B) = 0.72*E(B-V) +0.05*E(B-V)**2
    #
    # E(B-V) = (B-V) - (B-V)o
    # E(U-B) = (U-B) - (U-B)o
    #
    bv_intrsc = np.array(bv_obsrv) - e_bv
    e_ub = 0.72*e_bv + 0.05*e_bv**2
    ub_intrsc = np.array(ub_obsrv) - e_ub

    return bv_intrsc, ub_intrsc
    


def track_distance(track, bv_intrsc, ub_intrsc):
    '''
    Function that takes as input two arrays of (U-B)o and (B-V)o colors for
    a group of stars and returns for each one the point in the ZAMS closest
    to it.
    '''
    # Create color-color list.
    bv_ub = [bv_intrsc, ub_intrsc]

    # Get distances of *every* star to *every* point in this
    # ZAMS. The list is made of N sub-lists where N is the
    # number of cluster stars and M items in each sub-list where M is
    # the number of points in the function/track/isochrone.
    # dist_m = [[star_1], [star_2], ..., [star_N]]
    # [star_i] = [dist_1, dist_2, ..., dist_M]
    dist_m = np.array(sp.cdist(zip(*bv_ub), zip(*track), 'euclidean'))
                
    # Identify the position of the closest point in track for each star and
    # save the index pointing to it.
    min_dists = dist_m.min(axis=1)
    
    min_dist_indxs = []
    for indx, item in enumerate(min_dists):
        min_dist_indxs.append(dist_m[indx].tolist().index(item))
    
    return min_dist_indxs


# Read the files with the ZAMS to be used. Options are S-K, Girardi-Marigo 
# and Girardi PARSEC with a given metallicity.
ZAMS_file = 'zams_SK'

bv_o, ub_o, M_abs, sp_type = [], [], [], []
with open(ZAMS_file, mode="r") as z_f:
    for line in z_f:
        li=line.strip()
        # Jump comments.
        if not li.startswith("#"):
            reader = li.split()           
            bv_o.append(float(reader[2]))
            ub_o.append(float(reader[1]))
            M_abs.append(float(reader[3]))
            sp_type.append(str(reader[0]))
            
track = [bv_o, ub_o]


# Read the file that stores the photometric data for all stars.
data_file = 'data_input'
    
# Loads the data in 'myfile' as a list of N lists where N is the number of
# columns. Each of the N lists contains all the data for the column.
data = np.loadtxt(data_file, unpack=True)
id_star, m_obs, e_m, bv_obsrv, e_bv, ub_obsrv, e_ub, vi_obs, e_vi = data[0], \
data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8]


# Ask for E(B-V) value.
print 'Input extinction value to be used.'
e_bv = float(raw_input('E(B-V): '))


# Get intrinsec values for both colors.
bv_intrsc, ub_intrsc = intrsc_values(bv_obsrv, ub_obsrv, e_bv)


# For each star in its intrinsic position find the point in the ZAMS closest
# to it and assign that star to that position. This way we obtain the
# spectral type and the absolute magnitude. 
min_dist_indxs = track_distance(track, bv_intrsc, ub_intrsc)
    
print min_dist_indxs
  

M_abs_filt = [M_abs[indx] for indx in min_dist_indxs]
sp_type_stars = [sp_type[indx] for indx in min_dist_indxs]

A_v = 3.1*e_bv
dist_mod = (np.array(m_obs) - M_abs_filt)
dist = 10**(0.2*(dist_mod+5-A_v))

print dist
print sp_type_stars

f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

plt.ylim(max(ub_obsrv)+0.5, min(ub_obsrv)-0.5)
ax1.scatter(bv_obsrv, ub_obsrv, c='b')
ax1.plot(track[0], track[1], c='k', ls='-', lw=2.)

plt.ylim(max(ub_intrsc)+0.5, min(ub_intrsc)-0.5)
ax2.scatter(bv_intrsc, ub_intrsc, c='r')
ax2.plot(track[0], track[1], c='k', ls='-', marker='x', lw=0.5)

plt.show()