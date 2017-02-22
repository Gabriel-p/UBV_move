
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter


def star_size(mag, min_mag):
    '''
    Convert magnitudes into intensities and define sizes of stars in
    finding chart.
    '''
    # Scale factor.
    factor = 500. * (1 - 1 / (1 + 150 / len(mag) ** 0.85))
    return 0.1 + factor * 10 ** ((np.array(mag) - min_mag) / -2.5)


def dens_map():
    """
    """
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
    plt.xlabel(r'Dist (kpc)', fontsize=12)
    plt.ylabel(r'$E_{(B-V)}$', fontsize=14)
    # H_g is the 2D histogram with a gaussian filter applied
    h_g = gaussian_filter(hist, 1.5, mode='constant')
    x_max, y_max = np.unravel_index(h_g.argmax(), h_g.shape)
    d_m, ebv_m = np.average(xedges[x_max:x_max + 2]), \
        np.average(yedges[y_max:y_max + 2])
    text1 = r'$E_{(B-V)}^{max}\,=\,%0.2f$' '\n' % ebv_m
    text2 = r'$dist^{max}\,=\,%0.2f$' % d_m
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
    plt.savefig('output/' + clust_name + '_dens_map.png', dpi=300)
    print('Density maps plotted.')


def tcd_chart():
    """
    """
    fig = plt.figure(figsize=(20, 20))  # create the top-level container
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
    plt.scatter(bv_obsrv, ub_obsrv, c='b', lw=0.5, edgecolors='k', s=10.)

    ax2 = fig.add_subplot(222)
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
    plt.plot(bv_o, ub_o, c='k', ls='-')
    # Plot extinction line.
    plt.plot([-0.33, 1.17], [-1.2, -0.0075], c='k', lw=1.5, ls='--')
    # Plot error bars.
    plt.errorbar(
        bv_obs_uniq[0], ub_obs_uniq[0], yerr=ub_obs_uniq[1],
        xerr=bv_obs_uniq[1], ms=2., lw=0.3, c='b', fmt='.')
    # Plot corrected stars with unique solutions.
    plt.scatter(bv_int_uniq, ub_int_uniq, c='r', lw=0.5, edgecolors='k', s=20.)
    # Plot extinction lines.
    plt.plot([bv_int_uniq, bv_obs_uniq[0]],
             [ub_int_uniq, ub_obs_uniq[0]], lw=0.3, ls='-')

    ax3 = fig.add_subplot(223)
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
    for bv_star, ub_star, ext_star, dist_star in zip(*[
            bv_obsrv, ub_obsrv, extin_list, dist]):
        for e, d in zip(*[ext_star, dist_star]):
            if isinstance(e, float) and isinstance(d, float):
                if abs(e - ebv_m) <= 0.2 and abs(d - d_m) <= 0.2:
                    temp[0].append(bv_star)
                    temp[1].append(ub_star)
    plt.scatter(temp[0], temp[1], c='b', lw=0.5, s=20.)

    ax4 = fig.add_subplot(224)
    plt.xlabel('x', fontsize=18)
    plt.ylabel('y', fontsize=18)
    ax4.minorticks_on()
    ax4.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
             zorder=1)
    min_mag = min(m_obs)
    st_sizes_arr = star_size(m_obs, min_mag)
    plt.scatter(x_star, y_star, marker='o', c='black', s=st_sizes_arr,
                zorder=4)
    st_sizes_arr_uniq = star_size(m_obs_uniq, min_mag)
    plt.scatter(x_uniq, y_uniq, marker='o', facecolors='none', edgecolor='r',
                s=st_sizes_arr_uniq + 10., lw=.5, zorder=4)
    ax4.set_aspect('equal')

    fig.tight_layout()
    # Generate output plot file.
    plt.savefig('output/' + clust_name + '_CMD.png', dpi=300)
    print('Two-color diagrams and xy chart plotted.')


def v_segmented(clust_name, m_obs, bv_obsrv, ub_obsrv, bv_o, ub_o):
    """
    """
    #
    fig = plt.figure(figsize=(20, 20))
    v_max = math.floor(min(m_obs))
    v_range = [float(v_max + _) for _ in range(0, 9)]

    segments = []
    for _ in range(0, 9):
        idx = np.where((np.array(m_obs) >= float(v_max + _)) &\
            (np.array(m_obs) < float(v_max + _ + 1.)))
        segments.append([bv_obsrv[idx], ub_obsrv[idx]])
        # Store last V value.
        v_min = float(v_max + _ + 1.)
    v_range.append(math.ceil(max(m_obs)))

    # Add last segment
    idx = np.where((np.array(m_obs) >= v_min) & (np.array(m_obs) <= 1.e6))
    segments.append([bv_obsrv[idx], ub_obsrv[idx]])

    for i in range(1, 10):
        ax = fig.add_subplot(3, 3, i)

        plt.xlim(-0.7, 2.5)
        plt.ylim(1.7, -1.3)
        plt.xlabel('$(B-V)$', fontsize=14)
        plt.ylabel('$(U-B)$', fontsize=14)
        ax.minorticks_on()
        ax.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
        ax.text(0.7, 0.93, r'${:.0f} \leq \,V_{{obs}}\, < {:.0f}$'.format(
            v_range[i - 1], v_range[i]), transform=ax.transAxes,
            bbox=dict(facecolor='white', alpha=0.85), fontsize=14)
        plt.scatter(segments[i - 1][0], segments[i - 1][1], edgecolor='k',
                    lw=.5, zorder=4)
        # Plot ZAMS.
        plt.plot(bv_o, ub_o, c='k', ls='-', zorder=1)
        # Plot extinction line.
        plt.plot([-0.33, 1.17], [-1.2, -0.0075], c='k', lw=1.5, ls='--')

    fig.tight_layout()
    # Generate output plot file.
    plt.savefig('output/' + clust_name + '_segment.png', dpi=300)
    print('Two-color diagrams segmented by V intervals plotted.')


def main(clust_name, x_star, y_star, x_uniq, y_uniq, m_obs, m_obs_uniq,
         ext_dist_all, bv_o, ub_o, bv_obsrv, ub_obsrv, extin_max,
         bv_obs_uniq, ub_obs_uniq, bv_int_uniq, ub_int_uniq, extin_list, dist):
    """
    Generate final plots.
    """
    print('Creating output plots.')

    dens_map()
    # tcd_chart()
    v_segmented(clust_name, m_obs, bv_obsrv, ub_obsrv, bv_o, ub_o)
