
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


def v_segmented(clust_name, m_obs, bv_obsrv, ub_obsrv, bv_o, ub_o):
    """
    """
    #
    fig = plt.figure(figsize=(20, 20))

    v_max = math.floor(min(m_obs))
    v_1st_seg = float(math.ceil(sorted(m_obs)[:10][-1]))

    segments = []

    # Add first segment
    idx = np.where((np.array(m_obs) >= v_max) & (np.array(m_obs) <= v_1st_seg))
    segments.append([bv_obsrv[idx], ub_obsrv[idx]])

    for _ in range(0, 8):
        idx = np.where((np.array(m_obs) >= float(v_1st_seg + _)) &
                       (np.array(m_obs) < float(v_1st_seg + _ + 1.)))
        segments.append([bv_obsrv[idx], ub_obsrv[idx]])
        # Store last V value.
        v_min = float(v_1st_seg + _ + 1.)

    # Add last segment
    idx = np.where((np.array(m_obs) >= v_min) & (np.array(m_obs) <= 1.e6))
    segments.append([bv_obsrv[idx], ub_obsrv[idx]])

    v_range = [float(v_max)]
    v_range += [float(v_1st_seg + _) for _ in range(0, 8)]
    v_range.append(math.ceil(max(m_obs)))

    for i in range(1, 10):
        if segments[i - 1][0].any():
            ax = fig.add_subplot(3, 3, i)

            plt.xlim(-0.4, 2.6)
            plt.ylim(1.7, -1.3)
            plt.xlabel('$(B-V)$', fontsize=14)
            plt.ylabel('$(U-B)$', fontsize=14)
            ax.minorticks_on()
            ax.grid(b=True, which='major', color='gray', linestyle='-',
                    zorder=1)
            text = r'${:.0f} \leq \,V_{{obs}}\, < {:.0f}$' + '\nN={}'
            ax.text(
                .7, 0.9, text.format(
                    v_range[i - 1], v_range[i], len(segments[i - 1][0])),
                transform=ax.transAxes,
                bbox=dict(facecolor='white', alpha=0.85), fontsize=14)
            plt.scatter(segments[i - 1][0], segments[i - 1][1], edgecolor='k',
                        lw=.5, zorder=4)
            # Plot ZAMS.
            plt.plot(bv_o, ub_o, c='k', ls='-', zorder=2)
            # Plot extinction line.
            plt.plot([-0.33, 1.17], [-1.2, -0.0075], c='k', lw=1.5, ls='--')

    fig.tight_layout()
    # Generate output plot file.
    plt.savefig('output/' + clust_name + '_A.png', dpi=300)
    print('Two-color diagrams segmented by V intervals plotted.')


def dens_map(clust_name, d_max, e_max, hist, xedges, yedges):
    """
    """
    fig = plt.figure(figsize=(20, 10))  # create the top-level container

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
    plt.savefig('output/' + clust_name + '_B.png', dpi=300)
    print('Density maps plotted.')


def tcd_chart(clust_name, bv_o, ub_o, x_star, y_star, m_obs, bv_obsrv,
              ub_obsrv, extin_max, x_uniq, y_uniq, m_uniq, d_uniq, bv_obs_uniq,
              ub_obs_uniq, bv_int_uniq, ub_int_uniq):
    """
    """
    fig = plt.figure(figsize=(20, 20))

    ax1 = fig.add_subplot(221)
    ax1.set_title("TCD of all stars")
    plt.xlim(-0.4, 2.1)
    plt.ylim(1.7, -1.3)
    plt.xlabel('$(B-V)$', fontsize=18)
    plt.ylabel('$(U-B)$', fontsize=18)
    ax1.minorticks_on()
    ax1.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
    ax1.text(0.75, 0.95, 'N={}'.format(len(bv_obsrv)), transform=ax1.transAxes,
             bbox=dict(facecolor='white', alpha=0.85), fontsize=18)
    # Plot ZAMS.
    plt.plot(bv_o, ub_o, c='k', ls='-', zorder=2)
    # Plot extinction line.
    plt.plot([-0.33, 1.17], [-1.2, -0.0075], c='k', lw=1.5, ls='--')
    # Plot all observed stars.
    plt.scatter(bv_obsrv, ub_obsrv, c='b', lw=0.5, edgecolors='k', s=10.,
                zorder=4)

    ax2 = fig.add_subplot(222)
    ax2.set_title("TCD of stars with unique solutions on ZAMS")
    plt.xlim(-0.4, 2.1)
    plt.ylim(1.7, -1.3)
    plt.xlabel('$(B-V)$', fontsize=18)
    plt.ylabel('$(U-B)$', fontsize=18)
    ax2.minorticks_on()
    ax2.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
    text1 = r'$E_{{(B-V)}}^{{max}}\,=\,{:0.2f}$'.format(extin_max) + '\n'
    text2 = "N={}".format(len(bv_int_uniq))
    text = text1 + text2
    ax2.text(0.75, 0.92, text, transform=ax2.transAxes,
             bbox=dict(facecolor='white', alpha=0.85), fontsize=18)
    # Plot ZAMS.
    plt.plot(bv_o, ub_o, c='k', ls='-', zorder=2)
    # Plot extinction line.
    plt.plot([-0.33, 1.17], [-1.2, -0.0075], c='k', lw=1.5, ls='--')
    # Plot error bars.
    plt.errorbar(
        bv_obs_uniq[0], ub_obs_uniq[0], yerr=ub_obs_uniq[1],
        xerr=bv_obs_uniq[1], ms=2., lw=0.3, c='b', fmt='.')
    # Plot corrected stars with unique solutions.
    plt.scatter(bv_int_uniq, ub_int_uniq, c='r', lw=0.5, edgecolors='k', s=20.,
                zorder=4)
    # Plot extinction lines.
    plt.plot([bv_int_uniq, bv_obs_uniq[0]],
             [ub_int_uniq, ub_obs_uniq[0]], lw=0.3, ls='-')

    ax3 = fig.add_subplot(223)
    ax3.set_title("Stars w unique solutions colored according to d (kpc)")
    plt.xlim(max(min(bv_obsrv) - .2, -1.), max(bv_obsrv) + .2)
    plt.ylim(max(m_obs) + .5, min(m_obs) - .5)
    plt.xlabel('$(B-V)$', fontsize=18)
    plt.ylabel('$V$', fontsize=18)
    ax3.minorticks_on()
    ax3.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
    # Plot all stars
    plt.scatter(bv_obsrv, m_obs, c='grey', s=20., edgecolor='w', lw=.5,
                zorder=3)
    # Plot unique solution cluster stars.
    cm = plt.cm.get_cmap('RdYlBu')
    im = plt.scatter(bv_obs_uniq[0], m_uniq, c=d_uniq, lw=0.5, edgecolors='k',
                     s=40., cmap=cm, zorder=4)

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)

    # cbar = plt.colorbar()
    cbar.set_label("Dist (kpc)")
    # ax3.set_aspect('auto')

    ax4 = fig.add_subplot(224)
    ax4.set_title("Positional chart, stars w/ unique solutions in red")
    plt.xlabel('x', fontsize=18)
    plt.ylabel('y', fontsize=18)
    ax4.minorticks_on()
    ax4.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
             zorder=1)
    # All stars
    min_mag = min(m_obs)
    st_sizes_arr = star_size(m_obs, min_mag)
    plt.scatter(x_star, y_star, marker='o', c='black', s=st_sizes_arr,
                zorder=4)
    # Probable cluster members
    st_sizes_arr_uniq = star_size(m_uniq, min_mag)
    plt.scatter(x_uniq, y_uniq, marker='o', facecolors='none', edgecolor='r',
                s=st_sizes_arr_uniq + 10., lw=.5, zorder=4)
    ax4.set_aspect('equal')

    fig.tight_layout()
    # Generate output plot file.
    plt.savefig('output/' + clust_name + '_C.png', dpi=300)
    print('Stars with unique solutions plotted.')


def tcd_probs(clust_name, bv_o, ub_o, ebv_sig, dm_sig, dist, x_star, y_star,
              m_obs, bv_obsrv, E_BV, dist_kpc, x_prob, y_prob, m_prob,
              bv_prob, ub_prob):
    """
    """
    fig = plt.figure(figsize=(20, 20))

    ax1 = fig.add_subplot(221)
    ax1.set_title("TCD of probable member stars")
    plt.xlim(-0.4, 2.1)
    plt.ylim(1.7, -1.3)
    plt.xlabel('$(B-V)$', fontsize=18)
    plt.ylabel('$(U-B)$', fontsize=18)
    ax1.minorticks_on()
    ax1.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
    # Plot ZAMS.
    plt.plot(bv_o, ub_o, c='k', ls='-', zorder=2)
    # Plot moved ZAMS
    E_UB = 0.72 * E_BV + 0.05 * E_BV ** 2
    plt.plot(np.array(bv_o) + E_BV, np.array(ub_o) + E_UB, c='g', ls='--',
             zorder=2)
    # Plot extinction line.
    plt.plot([-0.33, 1.17], [-1.2, -0.0075], c='k', lw=1.5, ls='--')
    # Plot probable cluster stars.
    plt.scatter(bv_prob, ub_prob, c='b', lw=0.5, edgecolors='k', s=20.,
                zorder=4)
    text1 = '$E_{{(B-V)}}\,=\, {:0.2f} \pm {}$'.format(E_BV, ebv_sig) + '\n'
    text2 = '$dist\,=\, {:0.2f} \pm {}$'.format(dist_kpc, dm_sig) + '\n'
    text3 = "N={}".format(len(bv_prob))
    text = text1 + text2 + text3
    ax1.text(0.67, 0.9, text, transform=ax1.transAxes,
             bbox=dict(facecolor='white', alpha=0.85), fontsize=18)

    ax2 = fig.add_subplot(222)
    ax2.set_title("CMD of probable member stars")
    plt.xlim(max(min(bv_obsrv) - .2, -1.), max(bv_obsrv) + .2)
    plt.ylim(max(m_obs) + .5, min(m_obs) - .5)
    plt.xlabel('$(B-V)$', fontsize=18)
    plt.ylabel('$V$', fontsize=18)
    ax2.minorticks_on()
    ax2.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
    # Plot all stars
    plt.scatter(bv_obsrv, m_obs, c='grey', s=20., edgecolor='w', lw=.5,
                zorder=3)
    # Plot probable cluster stars.
    plt.scatter(bv_prob, m_prob, c='b', lw=0.5, edgecolors='k', s=20.,
                zorder=4)

    ax4 = fig.add_subplot(224)
    ax4.set_title("Positional chart, probable members in blue")
    plt.xlabel('x', fontsize=18)
    plt.ylabel('y', fontsize=18)
    ax4.minorticks_on()
    ax4.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
             zorder=1)
    # All stars
    min_mag = min(m_obs)
    st_sizes_arr = star_size(m_obs, min_mag)
    plt.scatter(x_star, y_star, marker='o', c='black', s=st_sizes_arr,
                zorder=4)
    # Probable cluster members
    st_sizes_arr_uniq = star_size(m_prob, min_mag)
    plt.scatter(x_prob, y_prob, marker='o', facecolors='none', edgecolor='b',
                s=st_sizes_arr_uniq + 10., lw=.5, zorder=4)
    ax4.set_aspect('equal')

    fig.tight_layout()
    # Generate output plot file.
    plt.savefig('output/' + clust_name + '_D.png', dpi=300)
    print('Probable member stars plotted.')


def main(clust_name, plot_v_seg, plot_dens_map, plot_tcd_uniq, plot_tcd_prob,
         x_star, y_star, m_obs, bv_o, ub_o, bv_obsrv, ub_obsrv, extin_max,
         ebv_sig, dm_sig, x_uniq, y_uniq, m_uniq, d_uniq, bv_obs_uniq,
         ub_obs_uniq, bv_int_uniq, ub_int_uniq, dist, d_max, e_max, hist,
         xedges, yedges, E_BV, dist_kpc, x_prob, y_prob, m_prob, bv_prob,
         ub_prob):
    """
    Generate final plots.
    """
    print('Creating output plots.')

    if plot_v_seg:
            v_segmented(clust_name, m_obs, bv_obsrv, ub_obsrv, bv_o, ub_o)
    if plot_dens_map:
        dens_map(clust_name, d_max, e_max, hist, xedges, yedges)
    if plot_tcd_uniq:
        tcd_chart(
            clust_name, bv_o, ub_o, x_star, y_star, m_obs, bv_obsrv,
            ub_obsrv, extin_max, x_uniq, y_uniq, m_uniq, d_uniq, bv_obs_uniq,
            ub_obs_uniq, bv_int_uniq, ub_int_uniq)
    if plot_tcd_prob:
        tcd_probs(clust_name, bv_o, ub_o, ebv_sig, dm_sig, dist, x_star,
                  y_star, m_obs, bv_obsrv, E_BV, dist_kpc, x_prob, y_prob,
                  m_prob, bv_prob, ub_prob)
