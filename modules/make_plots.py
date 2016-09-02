
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter


def main(clust_name, ext_dist_all, bv_o, ub_o, bv_obsrv, ub_obsrv, extin_max,
         bv_obs_uniq, ub_obs_uniq, bv_int_uniq, ub_int_uniq, extin_list, dist):
    """
    Generate final plots.
    """
    print('Creating output plots.')

    fig = plt.figure(figsize=(20, 10))  # create the top-level container

    # Maximum distance value in axis.
    d_max = np.sort(ext_dist_all[1])[int(0.9 * len(ext_dist_all[1]))]
    e_max = np.sort(ext_dist_all[0])[int(0.9 * len(ext_dist_all[0]))]

    hist, xedges, yedges = np.histogram2d(
        ext_dist_all[1], ext_dist_all[0],
        bins=[int(d_max / 0.002), int(max(ext_dist_all[0]) / 0.005)])

    # Plot density map.
    ax1 = fig.add_subplot(231)
    # Axis limits.
    plt.xlim(0, d_max)
    plt.ylim(0, e_max)
    plt.xlabel('Dist (kpc)', fontsize=12)
    plt.ylabel('$E_{(B-V)}$', fontsize=14)
    # H_g is the 2D histogram with a gaussian filter applied
    h_g = gaussian_filter(hist, 1.5, mode='constant')
    x_max, y_max = np.unravel_index(h_g.argmax(), h_g.shape)
    d_m, ebv_m = np.average(xedges[x_max:x_max + 2]), \
        np.average(yedges[y_max:y_max + 2])
    text1 = '$E_{(B-V)}^{max}\,=\,%0.2f$' '\n' % ebv_m
    text2 = '$dist^{max}\,=\,%0.2f$' % d_m
    text = text1 + text2
    ax1.text(0.7, 0.85, text, transform=ax1.transAxes,
             bbox=dict(facecolor='white', alpha=0.85), fontsize=14)
    plt.imshow(h_g.transpose(), origin='lower',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
               cmap=plt.get_cmap('jet'), aspect='auto')

    # Plot density map.
    ax2 = fig.add_subplot(232)
    # Axis limits.
    plt.xlim(0, d_max)
    plt.ylim(0, e_max)
    plt.xlabel('Dist (kpc)', fontsize=12)
    plt.ylabel('$E_{(B-V)}$', fontsize=14)
    # H_g is the 2D histogram with a gaussian filter applied
    h_g = gaussian_filter(hist, 2., mode='constant')
    x_max, y_max = np.unravel_index(h_g.argmax(), h_g.shape)
    d_m, ebv_m = np.average(xedges[x_max:x_max + 2]), \
        np.average(yedges[y_max:y_max + 2])
    text1 = '$E_{(B-V)}^{max}\,=\,%0.2f$' '\n' % ebv_m
    text2 = '$dist^{max}\,=\,%0.2f$' % d_m
    text = text1 + text2
    ax2.text(0.7, 0.85, text, transform=ax2.transAxes,
             bbox=dict(facecolor='white', alpha=0.85), fontsize=14)
    plt.imshow(h_g.transpose(), origin='lower',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
               cmap=plt.get_cmap('jet'), aspect='auto')

    # Plot density map.
    ax3 = fig.add_subplot(233)
    # Axis limits.
    plt.xlim(0, d_max)
    plt.ylim(0, e_max)
    plt.xlabel('Dist (kpc)', fontsize=12)
    plt.ylabel('$E_{(B-V)}$', fontsize=14)
    # H_g is the 2D histogram with a gaussian filter applied
    h_g = gaussian_filter(hist, 2.5, mode='constant')
    x_max, y_max = np.unravel_index(h_g.argmax(), h_g.shape)
    d_m, ebv_m = np.average(xedges[x_max:x_max + 2]), \
        np.average(yedges[y_max:y_max + 2])
    text1 = '$E_{(B-V)}^{max}\,=\,%0.2f$' '\n' % ebv_m
    text2 = '$dist^{max}\,=\,%0.2f$' % d_m
    text = text1 + text2
    ax3.text(0.7, 0.85, text, transform=ax3.transAxes,
             bbox=dict(facecolor='white', alpha=0.85), fontsize=14)
    plt.imshow(h_g.transpose(), origin='lower',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
               cmap=plt.get_cmap('jet'), aspect='auto')

    # Plot density map.
    ax4 = fig.add_subplot(234)
    # Axis limits.
    plt.xlim(0, d_max)
    plt.ylim(0, e_max)
    plt.xlabel('Dist (kpc)', fontsize=12)
    plt.ylabel('$E_{(B-V)}$', fontsize=14)
    # H_g is the 2D histogram with a gaussian filter applied
    h_g = gaussian_filter(hist, 3., mode='constant')
    x_max, y_max = np.unravel_index(h_g.argmax(), h_g.shape)
    d_m, ebv_m = np.average(xedges[x_max:x_max + 2]), \
        np.average(yedges[y_max:y_max + 2])
    text1 = '$E_{(B-V)}^{max}\,=\,%0.2f$' '\n' % ebv_m
    text2 = '$dist^{max}\,=\,%0.2f$' % d_m
    text = text1 + text2
    ax4.text(0.7, 0.85, text, transform=ax4.transAxes,
             bbox=dict(facecolor='white', alpha=0.85), fontsize=14)
    plt.imshow(h_g.transpose(), origin='lower',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
               cmap=plt.get_cmap('jet'), aspect='auto')

    # Plot density map.
    ax5 = fig.add_subplot(235)
    # Axis limits.
    plt.xlim(0, d_max)
    plt.ylim(0, e_max)
    plt.xlabel('Dist (kpc)', fontsize=12)
    plt.ylabel('$E_{(B-V)}$', fontsize=14)
    # H_g is the 2D histogram with a gaussian filter applied
    h_g = gaussian_filter(hist, 3.5, mode='constant')
    x_max, y_max = np.unravel_index(h_g.argmax(), h_g.shape)
    d_m, ebv_m = np.average(xedges[x_max:x_max + 2]), \
        np.average(yedges[y_max:y_max + 2])
    text1 = '$E_{(B-V)}^{max}\,=\,%0.2f$' '\n' % ebv_m
    text2 = '$dist^{max}\,=\,%0.2f$' % d_m
    text = text1 + text2
    ax5.text(0.7, 0.85, text, transform=ax5.transAxes,
             bbox=dict(facecolor='white', alpha=0.85), fontsize=14)
    plt.imshow(h_g.transpose(), origin='lower',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
               cmap=plt.get_cmap('jet'), aspect='auto')

    # Plot density map.
    ax6 = fig.add_subplot(236)
    # Axis limits.
    plt.xlim(0, d_max)
    plt.ylim(0, e_max)
    plt.xlabel('Dist (kpc)', fontsize=12)
    plt.ylabel('$E_{(B-V)}$', fontsize=14)
    # H_g is the 2D histogram with a gaussian filter applied
    h_g = gaussian_filter(hist, 4., mode='constant')
    x_max, y_max = np.unravel_index(h_g.argmax(), h_g.shape)
    d_m, ebv_m = np.average(xedges[x_max:x_max + 2]), \
        np.average(yedges[y_max:y_max + 2])
    text1 = '$E_{(B-V)}^{max}\,=\,%0.2f$' '\n' % ebv_m
    text2 = '$dist^{max}\,=\,%0.2f$' % d_m
    text = text1 + text2
    ax6.text(0.7, 0.85, text, transform=ax6.transAxes,
             bbox=dict(facecolor='white', alpha=0.85), fontsize=14)
    plt.imshow(h_g.transpose(), origin='lower',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
               cmap=plt.get_cmap('jet'), aspect='auto')

    fig.tight_layout()
    # Generate output plot file.
    plt.savefig(clust_name + '_dens_map.png', dpi=150)
    print('Plot 1 created.')

    fig = plt.figure(figsize=(30, 10))  # create the top-level container
    ax1 = fig.add_subplot(131)
    plt.xlim(-0.7, 2.5)
    plt.ylim(1.7, -1.3)
    plt.xlabel('$(B-V)$', fontsize=18)
    plt.ylabel('$(U-B)$', fontsize=18)
    ax1.minorticks_on()
    ax1.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
    # Plot ZAMS.
    plt.plot(bv_o, ub_o, c='k', ls='-')
    # Plot extinction line.
    plt.plot([-0.33, 1.17], [-1.2, -0.0075], c='k', lw=1.5, ls='--')
    # Plot all observed stars.
    plt.scatter(bv_obsrv, ub_obsrv, c='b', lw=0.5, s=10.)

    ax2 = fig.add_subplot(132)
    plt.xlim(-0.7, 2.5)
    plt.ylim(1.7, -1.3)
    plt.xlabel('$(B-V)$', fontsize=18)
    plt.ylabel('$(U-B)$', fontsize=18)
    ax2.minorticks_on()
    ax2.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
    ax2.text(0.67, 0.93, '$E_{(B-V)}^{max}\,=\,%0.2f$' % extin_max,
             transform=ax2.transAxes, bbox=dict(facecolor='white', alpha=0.85),
             fontsize=18)
    # Plot ZAMS.
    # plt.scatter(track[0], track[1], c='k', marker='x', lw=0.5, s=30.)
    plt.plot(bv_o, ub_o, c='k', ls='-')
    # Plot extinction line.
    plt.plot([-0.33, 1.17], [-1.2, -0.0075], c='k', lw=1.5, ls='--')
    # Plot stars with unique solutions.
    # for indx, star in enumerate(bv_int_uniq):
    # Plot error bars.
    plt.errorbar(
        bv_obs_uniq[0], ub_obs_uniq[0], yerr=ub_obs_uniq[1],
        xerr=bv_obs_uniq[1], fmt='.', ms=2., lw=0.3, c='b')
    # Plot corrected stars with unique solutions.
    plt.scatter(bv_int_uniq, ub_int_uniq, c='r', lw=0.5, s=30.)
    # Plot extinction lines.
    plt.plot([bv_int_uniq, bv_obs_uniq[0]],
             [ub_int_uniq, ub_obs_uniq[0]], lw=0.3, ls='-')

    ax3 = fig.add_subplot(133)
    plt.xlim(-0.7, 2.5)
    plt.ylim(1.7, -1.3)
    plt.xlabel('$(B-V)$', fontsize=18)
    plt.ylabel('$(U-B)$', fontsize=18)
    ax3.minorticks_on()
    ax3.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
    text1 = '$E_{(B-V)}\,=\,%0.2f\pm 0.2$' '\n' % ebv_m
    text2 = '$dist\,=\,%0.2f\pm 0.2$' % d_m
    text = text1 + text2
    ax3.text(0.67, 0.93, text, transform=ax3.transAxes,
             bbox=dict(facecolor='white', alpha=0.85), fontsize=18)
    # Plot ZAMS.
    plt.plot(bv_o, ub_o, c='k', ls='-')
    # Plot extinction line.
    plt.plot([-0.33, 1.17], [-1.2, -0.0075], c='k', lw=1.5, ls='--')
    # Plot probable cluster stars.
    temp = [[], []]
    for indx, bv_star in enumerate(bv_obsrv):
        for i in range(4):
            if type(extin_list[indx][i]) is float and\
                    type(dist[indx][i]) is float:
                if abs(extin_list[indx][i] - ebv_m) <= 0.2 and \
                        abs(dist[indx][i] - d_m) <= 0.2:
                    temp[0].append(bv_star)
                    temp[1].append(ub_obsrv[indx])

    plt.scatter(temp[0], temp[1], c='b', lw=0.5, s=20.)

    fig.tight_layout()
    # Generate output plot file.
    plt.savefig(clust_name + '_CMD.png', dpi=150)
    print('Plot 2 created.')
