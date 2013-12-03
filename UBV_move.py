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
    intrinsic lists according to a given extinction value.
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
    


def track_distance(track, bv_intrsc, ub_intrsc, e_bv, e_ub):
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
    # the number of points in the ZAMS/track/isochrone.
    # dist_m = [[star_1], [star_2], ..., [star_N]]
    # [star_i] = [dist_1, dist_2, ..., dist_M]
    dist_m = np.array(sp.cdist(zip(*bv_ub), zip(*track), 'euclidean'))
                
    # Identify the position of the closest point in track for each star and
    # save the distance to it.
    min_dists = dist_m.min(axis=1)
    
    # Save index of closest point in ZAMS to each star along with the
    # distance value.
    min_dist_indxs = []
    for indx, dist in enumerate(min_dists):
        # Check if distance star-ZAMS_point falls inside the star's error
        # rectangle.
        err_max = np.sqrt(e_bv[indx]**2+e_ub[indx]**2)
        if dist <= err_max:
            min_dist_indxs.append([dist_m[indx].tolist().index(dist), dist])
        else:
            min_dist_indxs.append([99999, 9999.])
    
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
data_file = 'data_input.dat'
    
# Loads the data in 'myfile' as a list of N lists where N is the number of
# columns. Each of the N lists contains all the data for the column.
data = np.loadtxt(data_file, unpack=True)
id_star, x_star, y_star, m_obs, e_m, bv_obsrv, e_bv, ub_obsrv, e_ub, vr_obs, \
e_vr, vi_obs, e_vi, r_star = data[0], data[1], data[2], data[3], data[4], data[5], \
data[6], data[7], data[8], data[9], data[10], data[11], data[12], data[13]


# Define range for E(B-V) value.
extin_max = 2.
extin_range = np.arange(0, extin_max, 0.01)


zams_indx_dist = [[99999, 9999.] for _ in range(len(id_star))]
extin_list = [99.9 for _ in range(len(id_star))]
# Loop through all e_bv values in range.
for extin in extin_range:
    
    # Get intrinsec values for both colors.
    bv_intrsc, ub_intrsc = intrsc_values(bv_obsrv, ub_obsrv, extin)
    
    # For each star in its intrinsic position find the point in the ZAMS closest
    # to it and store that point's index in the ZAMS and the distance star-point.
    # The point must fall inside the star's error rectangle.
    min_dist_indxs = track_distance(track, bv_intrsc, ub_intrsc, e_bv, e_ub)

    
    for indx, item in enumerate(min_dist_indxs):
        dist_new = item[1]
        dist_old = zams_indx_dist[indx][1]
        if dist_new < dist_old:
            zams_indx_dist[indx][0] = item[0]
            zams_indx_dist[indx][1] = item[1]
            extin_list[indx] = extin

ub_intrsc, bv_intrsc, M_abs_final, sp_type_final, dist = [], [], [], [], []
for indx, star in enumerate(zams_indx_dist):
    index = star[0]
    if index != 99999:
        ub_intrsc.append(track[1][index])
        bv_intrsc.append(track[0][index])
        M_abs_final.append(M_abs[index])
        sp_type_final.append(sp_type[index])
        # Calculate distance.
        A_v = 3.1*extin
        dist_mod = m_obs[indx] - M_abs[index]
        dist.append(10**(0.2*(dist_mod+5-A_v)))
    else:
        ub_intrsc.append(ub_obsrv[indx])
        bv_intrsc.append(bv_obsrv[indx])


# Plots.
import matplotlib.gridspec as gridspec
# figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 = 
# y1/y2 = 2.5 
fig = plt.figure(figsize=(30, 20)) # create the top-level container
gs1 = gridspec.GridSpec(12, 12)  # create a GridSpec object


ax1 = plt.subplot(gs1[0:6, 0:4])
plt.xlim(-0.7, 2.5)
plt.ylim(1.7, -1.3)
plt.xlabel('(B-V)', fontsize=12)
plt.ylabel('(U-B)', fontsize=12)
ax1.minorticks_on()
ax1.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
# Plot ZAMS.
ax1.scatter(track[0], track[1], c='k', marker='x', lw=0.5)
ax1.plot(track[0], track[1], c='k', ls='-')
# Plot stars.
for indx, star in enumerate(bv_intrsc):
    if bv_obsrv[indx] < 0.5:
        ax1.scatter(bv_intrsc[indx], ub_intrsc[indx], c='r', lw=0.5)
        ax1.errorbar(bv_obsrv[indx], ub_obsrv[indx], yerr=e_ub[indx],
                     xerr=e_bv[indx], fmt='.', ms=2., lw=0.5, c='b')    
        ax1.plot([bv_intrsc[indx], bv_obsrv[indx]],
                 [ub_intrsc[indx], ub_obsrv[indx]], lw=0.2, ls='-')


ax2 = plt.subplot(gs1[0:6, 4:8])
plt.xlim(-0.7, 2.5)
plt.ylim(1.7, -1.3)
plt.xlabel('(B-V)', fontsize=12)
plt.ylabel('(U-B)', fontsize=12)
ax2.minorticks_on()
ax2.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
# Plot ZAMS.
ax2.scatter(track[0], track[1], c='k', marker='x', lw=0.5)
ax2.plot(track[0], track[1], c='k', ls='-')
# Plot stars.
for indx, star in enumerate(bv_intrsc):

    if 0.5 <= bv_obsrv[indx] < 0.75:
        ax2.scatter(bv_intrsc[indx], ub_intrsc[indx], c='r', lw=0.5)
        ax2.errorbar(bv_obsrv[indx], ub_obsrv[indx], yerr=e_ub[indx],
                     xerr=e_bv[indx], fmt='.', ms=2., lw=0.5, c='b')    
        ax2.plot([bv_intrsc[indx], bv_obsrv[indx]],
                 [ub_intrsc[indx], ub_obsrv[indx]], lw=0.2, ls='-')
                 
                 
ax3 = plt.subplot(gs1[0:6, 8:12])
plt.xlim(-0.7, 2.5)
plt.ylim(1.7, -1.3)
plt.xlabel('(B-V)', fontsize=12)
plt.ylabel('(U-B)', fontsize=12)
ax3.minorticks_on()
ax3.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
# Plot ZAMS.
ax3.scatter(track[0], track[1], c='k', marker='x', lw=0.5)
ax3.plot(track[0], track[1], c='k', ls='-')
# Plot stars.
for indx, star in enumerate(bv_intrsc):
    if 0.75 <= bv_obsrv[indx] < 1.:
        ax3.scatter(bv_intrsc[indx], ub_intrsc[indx], c='r', lw=0.5)
        ax3.errorbar(bv_obsrv[indx], ub_obsrv[indx], yerr=e_ub[indx],
                     xerr=e_bv[indx], fmt='.', ms=2., lw=0.5, c='b')    
        ax3.plot([bv_intrsc[indx], bv_obsrv[indx]],
                 [ub_intrsc[indx], ub_obsrv[indx]], lw=0.2, ls='-')                 
                 


ax4 = plt.subplot(gs1[6:12, 0:4])
plt.xlim(-0.7, 2.5)
plt.ylim(1.7, -1.3)
plt.xlabel('(B-V)', fontsize=12)
plt.ylabel('(U-B)', fontsize=12)
ax4.minorticks_on()
ax4.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
# Plot ZAMS.
ax4.scatter(track[0], track[1], c='k', marker='x', lw=0.5)
ax4.plot(track[0], track[1], c='k', ls='-')
# Plot stars.
for indx, star in enumerate(bv_intrsc):
    if 1. <= bv_obsrv[indx] < 1.25:
        ax4.scatter(bv_intrsc[indx], ub_intrsc[indx], c='r', lw=0.5)
        ax4.errorbar(bv_obsrv[indx], ub_obsrv[indx], yerr=e_ub[indx],
                     xerr=e_bv[indx], fmt='.', ms=2., lw=0.5, c='b')    
        ax4.plot([bv_intrsc[indx], bv_obsrv[indx]],
                 [ub_intrsc[indx], ub_obsrv[indx]], lw=0.2, ls='-')  


ax5 = plt.subplot(gs1[6:12, 4:8])
plt.xlim(-0.7, 2.5)
plt.ylim(1.7, -1.3)
plt.xlabel('(B-V)', fontsize=12)
plt.ylabel('(U-B)', fontsize=12)
ax5.minorticks_on()
ax5.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
# Plot ZAMS.
ax5.scatter(track[0], track[1], c='k', marker='x', lw=0.5)
ax5.plot(track[0], track[1], c='k', ls='-')
# Plot stars.
for indx, star in enumerate(bv_intrsc):
    if 1.25 <= bv_obsrv[indx] < 1.5:
        ax5.scatter(bv_intrsc[indx], ub_intrsc[indx], c='r', lw=0.5)
        ax5.errorbar(bv_obsrv[indx], ub_obsrv[indx], yerr=e_ub[indx],
                     xerr=e_bv[indx], fmt='.', ms=2., lw=0.5, c='b')    
        ax5.plot([bv_intrsc[indx], bv_obsrv[indx]],
                 [ub_intrsc[indx], ub_obsrv[indx]], lw=0.2, ls='-') 
                 
                 
ax6 = plt.subplot(gs1[6:12, 8:12])
plt.xlim(-0.7, 2.5)
plt.ylim(1.7, -1.3)
plt.xlabel('(B-V)', fontsize=12)
plt.ylabel('(U-B)', fontsize=12)
ax6.minorticks_on()
ax6.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
# Plot ZAMS.
ax6.scatter(track[0], track[1], c='k', marker='x', lw=0.5)
ax6.plot(track[0], track[1], c='k', ls='-')
# Plot stars.
for indx, star in enumerate(bv_intrsc):
    if 1.5 <= bv_obsrv[indx]:
        ax6.scatter(bv_intrsc[indx], ub_intrsc[indx], c='r', lw=0.5)
        ax6.errorbar(bv_obsrv[indx], ub_obsrv[indx], yerr=e_ub[indx],
                     xerr=e_bv[indx], fmt='.', ms=2., lw=0.5, c='b')    
        ax6.plot([bv_intrsc[indx], bv_obsrv[indx]],
                 [ub_intrsc[indx], ub_obsrv[indx]], lw=0.2, ls='-') 
                 

fig.tight_layout()

# Generate output file.
plt.savefig('output_CMD.png', dpi=300)

print 'End.'
