#!/usr/bin/env python3
#
#
# ==========================================================================
# 
# figure_backbone_2098comp.py
# 
# ==========================================================================
# 11/18/2020    sij    created 
#
#-----------------------------------------------------------------------------
""" make figures for results section of AGU 2020 poster """

import numpy as np
from calc_backbone_metric import read_tomo
from calc_backbone_metric import calc_backbone_metric
from calc_backbone_metric import find_neutral_line_coords
import os
import math
from math import sin
from math import cos
import argparse
from astropy.io import fits
from scipy.linalg import norm
from scipy.signal import find_peaks
import csv
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt



def identify_tomography_peaks(tomofile):
    """ identify ridges of higher electron density in tomography-based slices """
    # using this version instead of calc_prediction_metric.py version to record the parameters used in paper

    # read tomography file
    tomo_lon, tomo_lat, e_density_map = read_tomo(tomofile)
    emap = np.asarray(e_density_map)
    emap[emap < 0.] = 0.

    # find peaks in electron density slice
    streamer_lats=[]
    streamer_lons=[]
    for lon_ind in range(len(tomo_lon)):
        lon_strip = emap[:,lon_ind]
        biggest_diff = np.max(abs(np.diff(lon_strip)))
        if biggest_diff >= 0.00000001 and lon_strip.max() > emap.max() / 10.:
            peaks,_ = find_peaks(lon_strip, prominence = lon_strip.max() / 30., height = emap.max() / 8., \
                               width = 8.)
            #peaks,_ = find_peaks(lon_strip)
            streamer_lats.extend([tomo_lat[peak] for peak in peaks])
            streamer_lons.extend([tomo_lon[lon_ind] for peak in peaks])

    for lat_ind in range(len(tomo_lat)):
        wrap_end = round(10. / 360. * len(tomo_lon)) #tack on some extra data so we can find peaks near the edges
        lat_strip = np.append(emap[lat_ind, :], emap[lat_ind, : wrap_end])
        #lat_strip.append(emap[lat_ind, : wrap_end])
        biggest_diff = np.max(abs(np.diff(lat_strip)))
        if biggest_diff >= 0.0000001 and lat_strip.max() > emap.max() / 10.:
             peaks,_ = find_peaks(lat_strip, prominence = lat_strip.max() / 30., \
                               height = emap.max() / 8., width = 8.)
             peaks[peaks >= len(tomo_lon)] -= len(tomo_lon)
             streamer_lats.extend([tomo_lat[lat_ind] for peak in peaks])
             streamer_lons.extend([tomo_lon[peak] for peak in peaks])


    # consolidate peaks to grid resolution
    consolidated_lats = []
    consolidated_lons = []
    emap_blank = np.zeros((len(tomo_lat),len(tomo_lon)))
    grid_res = 360. / (len(tomo_lon) + 1)
    for slat,slon in zip(streamer_lats, streamer_lons):
        closest_lat = np.argmin([abs(slat - tlat) for tlat in tomo_lat])
        closest_lon = np.argmin([abs(slon - tlon) for tlon in tomo_lon])
        emap_blank[closest_lat, closest_lon] = 1
    nonzeroes = emap_blank.nonzero()
    consolidated_lats = [tomo_lat[lat_ind] for lat_ind in nonzeroes[0]]
    consolidated_lons = [tomo_lon[lon_ind] for lon_ind in nonzeroes[1]]

    return consolidated_lons, consolidated_lats


def comparison_figure(cs_xs, cs_ys, streamer_xs, streamer_ys, tomofile, altitude, ax, title=False):

    # calculate backbone metric value
    backbone_metric = calc_backbone_metric(cs_xs, cs_ys, streamer_xs, streamer_ys)
    backbone_str = 'M = ' + str(round(backbone_metric, 2))
    nneutral_str = f'N = {len(cs_xs)}'
    alt_str = 'Alt. = ' + str(altitude)

    # get electron density map
    tomo_lon, tomo_lat, e_density_map = read_tomo(tomofile)
    emap = np.asarray(e_density_map)

    # add image and plots
    ax.plot(streamer_xs, streamer_ys,'b.', label='Tomo.-based Streamer')
    ax.plot(cs_xs,cs_ys,'rx',label='WSA Neutral Line')
    ax.imshow(emap, origin='lower', extent=[tomo_lon[0], tomo_lon[-1], tomo_lat[0], \
              tomo_lat[-1]], aspect = 'auto')
    ax.set_ylabel('Latitude (deg)',fontsize=18.)
    ax.set_xlabel('Carrington Longitude (deg)', fontsize=18.)
    ax.annotate(alt_str, (30., -80), color='white', fontsize=18.)
    ax.annotate(nneutral_str, (130., -80), color='white', fontsize=18.)
    ax.annotate(backbone_str, (230., -80), color='white', fontsize=18.)
    ax.set_xticks([0, 45, 90., 135, 180., 225, 270., 315])
    ax.tick_params(labelsize=16.)
    if title:
        ax.set_title(title, fontsize=22)

    return ax


def figure_2098comp(tomofile1, tomofile2, tomofile3, bcbfile1, bcbfile2, altitude1, \
             altitude2, altitude3, outfile):

    # find coordinates of current sheet from all bcbfiles
    cs_xs1_1, cs_ys1_1 = find_neutral_line_coords(bcbfile1, altitude1)
    cs_xs2_1, cs_ys2_1 = find_neutral_line_coords(bcbfile2, altitude1)
    cs_xs1_2, cs_ys1_2 = find_neutral_line_coords(bcbfile1, altitude2)
    cs_xs2_2, cs_ys2_2 = find_neutral_line_coords(bcbfile2, altitude2)
    cs_xs1_3, cs_ys1_3 = find_neutral_line_coords(bcbfile1, altitude3)
    cs_xs2_3, cs_ys2_3 = find_neutral_line_coords(bcbfile2, altitude3)

    # find coordinates of neutral lines from tomofiles
    streamer_xs1, streamer_ys1 = identify_tomography_peaks(tomofile1)
    streamer_xs2, streamer_ys2 = identify_tomography_peaks(tomofile2)
    streamer_xs3, streamer_ys3 = identify_tomography_peaks(tomofile3)

    # get title strings
    basefile1=os.path.basename(bcbfile1)
    basefile2=os.path.basename(bcbfile2)
    title1=basefile1[4:12] + ' ' + basefile1[-9:-5]
    title2=basefile2[4:12] + ' ' + basefile2[-9:-5]

    # create figure 
    fig, axs = plt.subplots(3, 2, figsize = [18,13], squeeze=False, sharex='col')

    # add figure for all tomofiles, bcbfiles
    ax = comparison_figure(cs_xs1_1, cs_ys1_1, streamer_xs1, streamer_ys1, \
                      tomofile1, altitude1, axs[0][0], title=title1)
    ax = comparison_figure(cs_xs2_1, cs_ys2_1, streamer_xs1, streamer_ys1, \
                      tomofile1, altitude1, axs[0][1], title=title2)
    ax = comparison_figure(cs_xs1_2, cs_ys1_2, streamer_xs2, streamer_ys2, \
                      tomofile2, altitude2, axs[1][0])
    ax = comparison_figure(cs_xs2_2, cs_ys2_2, streamer_xs2, streamer_ys2, \
                      tomofile2, altitude2, axs[1][1])
    ax = comparison_figure(cs_xs1_3, cs_ys1_3, streamer_xs3, streamer_ys3, \
                      tomofile3, altitude3, axs[2][0])
    ax = comparison_figure(cs_xs2_3, cs_ys2_3, streamer_xs3, streamer_ys3, \
                      tomofile3, altitude3, axs[2][1])

    # add figure decorations
    fig.tight_layout()

    # save figure
    fig.savefig(outfile)
    print(f'Wrote: {outfile}')

    return



if __name__ == '__main__':

    ARG_PARSER = argparse.ArgumentParser()
    ARG_PARSER.add_argument('-t1', '--tomofile1', action='store', default='', dest='tomofile1')
    ARG_PARSER.add_argument('-t2', '--tomofile2', action='store', default='', dest='tomofile2')
    ARG_PARSER.add_argument('-t3', '--tomofile3', action='store', default='', dest='tomofile3')
    ARG_PARSER.add_argument('-a1', '--altitude1', action='store', default='', dest='altitude1', type=float)
    ARG_PARSER.add_argument('-a2', '--altitude2', action='store', default='', dest='altitude2', type=float)
    ARG_PARSER.add_argument('-a3', '--altitude3', action='store', default='', dest='altitude3', type=float)
    ARG_PARSER.add_argument('-b1', '--bcbfile1', action='store', default='', dest='bcbfile1')
    ARG_PARSER.add_argument('-b2', '--bcbfile2', action='store', default='', dest='bcbfile2')
    ARG_PARSER.add_argument('-o', '--outfile', action='store')

    ARGS = ARG_PARSER.parse_args()
    tomofile1 = ARGS.tomofile1
    tomofile2 = ARGS.tomofile2
    tomofile3 = ARGS.tomofile3
    bcbfile1 = ARGS.bcbfile1
    bcbfile2 = ARGS.bcbfile2
    altitude1 = ARGS.altitude1
    altitude2 = ARGS.altitude2
    altitude3 = ARGS.altitude3

    # each tomofile corresponds to the same tomography run, but at a different altitude;
    #    the altitudes need to correspond to those altitudes.  the tomofiles are created 
    #    with the IDL program digest_tomography.py and I usually save them in the TOMO_COMP/DATA/
    #    directory

    figure_2098comp(tomofile1, tomofile2, tomofile3, bcbfile1, bcbfile2, altitude1, altitude2, \
                     altitude3, ARGS.outfile)


    exit(0)
