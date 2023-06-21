#!/usr/bin/env python3
#
#
# ==========================================================================
# 
# figure_tomography_peaks.py
# 
# ==========================================================================
#  11/11/2020   sij    created 
#
#-----------------------------------------------------------------------------
""" makes a figure showing the identification of peaks in electron density data slice """

import argparse
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
import calc_backbone_metric


# electron density peak-finding parameters
mindiff = 0.0001
min_series_max = 10.
prom_scale = 30.
height_scale = 8.
min_peak_width = 8.


def make_tomography_peaks_figure(tomofile, outfile = \
               '/home/sjonesme/Desktop/meeting_image/peak_detection_figure.jpg', \
               ind = 105):

    # read tomography file
    tomo_lon, tomo_lat, e_density_map = calc_backbone_metric.read_tomo(tomofile)
    emap = np.asarray(e_density_map)
    emap[emap < 0.] = 0.
    alt = tomofile[-7:-4]

    # find peaks in electron density slice
    lon_strip = emap[:, ind]
    peaks,_ = find_peaks(lon_strip, prominence = lon_strip.max() / prom_scale, \
                               height = emap.max() / height_scale, \
                               width = min_peak_width)


    # create figure
    fig = plt.figure(figsize=(15,8))
    gs = fig.add_gridspec(ncols=2,nrows=1,width_ratios=[15, 1])

    # plot data slice
    ax = fig.add_subplot(gs[0])
    im = ax.imshow(emap, origin='lower', extent=[tomo_lon[0], tomo_lon[-1], \
                               tomo_lat[0], tomo_lat[-1]], aspect="auto")
    ax.plot((tomo_lon[ind], tomo_lon[ind]), (tomo_lat[0], tomo_lat[-1]), color='orange')
    ax.set_ylabel('Latitude', fontsize=26)
    ax.set_xlabel('Carrington Longitude', fontsize=26)
    ax.set_title('Electron Density Slice at ' + alt + ' $R_{sun}$', fontsize=28)
    ax.tick_params(labelsize=18)

    # plot detected peaks
    ax1 = fig.add_subplot(gs[1])
    ax1.plot(lon_strip, tomo_lat)
    ax1.plot([lon_strip[peak] for peak in peaks], [tomo_lat[peak] for peak in peaks], 'b*')
    ax1.set_xticklabels([])
    ax1.tick_params(labelsize=18)
    ax1.yaxis.tick_right()
    ax1.set_ylim(tomo_lat[0], tomo_lat[-1])
    plt.subplots_adjust(wspace=0.05, hspace=0)

    ## add colorbar
    axs=[ax, ax1]
    cbar = fig.colorbar(im, ax=axs, location='bottom', pad=0.15)
    cbar.set_label('Density ($10^{7}cm^{-3}$)', fontsize=20)
    cbar.ax.tick_params(labelsize=18)

    # save results
    plt.savefig(outfile)
    print(f'Wrote {outfile}')

    return

 

if __name__ == '__main__':

    ARG_PARSER = argparse.ArgumentParser()
    ARG_PARSER.add_argument('-t', '--tomofile', action='store')
    ARG_PARSER.add_argument('-o', '--outfile', action='store', \
              default='/home/sjonesme/Desktop/meeting_image/peak_detection_figure.jpg')
    ARG_PARSER.add_argument('-i', '--ind', action='store', type=int, default=105)
    ARGS = ARG_PARSER.parse_args()

    make_tomography_peaks_figure(ARGS.tomofile, ARGS.outfile, ARGS.ind)
