#!/usr/bin/env python3
#
#
# ==========================================================================
# 
# compare_streamer_lats.py
# 
# ==========================================================================
# 12/17/2020    sij    copied from compare_streamer_lats.py 
#
#-----------------------------------------------------------------------------
""" routines to calculate value of Backbone metric """

import numpy as np
import os
import math
import argparse
from astropy.io import fits
from scipy.linalg import norm
from scipy.signal import find_peaks
import csv
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import read_wsa_bfield

# electron density peak-finding parameters
mindiff = 0.0001 
min_series_max = 10.
prom_scale = 30.
height_scale = 8.
min_peak_width = 8.


def angular_separation(phi0, theta0, phi1, theta1, degrees=False, cossep=False):
    """ calculates distance between two points on a sphere, in radians """


    if degrees:    # convert to radians, if needed
        first_term = math.sin(math.radians(theta0))*math.sin(math.radians(theta1))
        second_term = math.cos(math.radians(theta0)) * math.cos(math.radians(theta1)) * \
                  math.cos(math.radians(phi0-phi1))
    else:
        first_term = math.sin(theta0)*math.sin(theta1)
        second_term = math.cos(theta0) * math.cos(theta1) * math.cos(phi0-phi1)

    try:
        dist = np.arccos(first_term + second_term)
    except:
        print(first_term + second_term)

    if cossep:     # return the cosine of the separation instead
        return first_term + second_term

    return dist


def find_neutral_line_coords(bcbfile, altitude):
    """ finds neutral lines in WSA magnetic field model (bcb_* file) slices """

    # read bcb_* WSA output file, shift to carrington frame
    wsa_bfield = read_wsa_bfield.read_wsa_bfield(bcbfile)
    wsa_bfield.carrington_frame()

    # find altitude slice closest to desired altitude
    closest_alt = np.argmin([abs(altitude - rad) for rad in wsa_bfield.radii])
    B_alt = wsa_bfield.br[closest_alt]

    # translate to relative strength of radial field component
    B_tot = norm(B_alt, axis = 0)
    relative_Br_alt = np.divide(B_alt[0], B_tot)

    # use matplotlib contour routine to get contour coordinates
    # Br is the first slice (dimension 2 of four), radial coordinate is dimension 1
    xs = []
    ys = []
    contour_obj = plt.contour(wsa_bfield.lons, wsa_bfield.lats, relative_Br_alt, levels=[0.0])
    for segs in contour_obj.allsegs[0]:
        xs.extend(segs[:,0])
        ys.extend(segs[:,1])

    return(xs,ys)


def read_tomo(filename):
    """ reads .csv files with slice at fixed radius created by digest_tomography.pro """

    tomo_lat = []
    e_density_map = []
    index = 0
    with open(filename) as tomo_file:
        ne_datablock = csv.reader(tomo_file)
        for row in tomo_file:
            row_array = row.split(',')
            if index == 0:
                tomo_lon = [float(val) for val in row_array[1:-1]]
            else:
                tomo_lat.append(float(row_array[0]))
                e_density_map.append([float(val) for val in row_array[1:-1]])
            index += 1
    return(tomo_lon, tomo_lat, e_density_map)


def identify_tomography_peaks(tomofile):
    """ identify ridges of higher electron density in tomography-based slices """

    # read tomography file
    tomo_lon, tomo_lat, e_density_map = read_tomo(tomofile)
    emap = np.asarray(e_density_map)
    emap[emap < 0.] = 0.

    # find peaks in electron density slice
    streamer_lats=[]
    streamer_lons=[]
    # loop over every column of image
    for lon_ind in range(len(tomo_lon)):
        lon_strip = emap[:,lon_ind]
        biggest_diff = np.max(abs(np.diff(lon_strip)))
        if biggest_diff >= mindiff and lon_strip.max() > emap.max() / min_series_max:
            peaks,_ = find_peaks(lon_strip, prominence = lon_strip.max() / prom_scale, \
                               height = emap.max() / height_scale, \
                               width = min_peak_width)
            #peaks,_ = find_peaks(lon_strip)
            streamer_lats.extend([tomo_lat[peak] for peak in peaks])
            streamer_lons.extend([tomo_lon[lon_ind] for peak in peaks])

    # loop over every row of image
    for lat_ind in range(len(tomo_lat)):
        wrap_end = round(10. / 360. * len(tomo_lon)) #tack on some extra data so we can find peaks near the edges
        lat_strip = np.append(emap[lat_ind, :], emap[lat_ind, : wrap_end])
        #lat_strip.append(emap[lat_ind, : wrap_end])
        biggest_diff = np.max(abs(np.diff(lat_strip)))
        if biggest_diff >= mindiff and lat_strip.max() > emap.max() / min_series_max:
             peaks,_ = find_peaks(lat_strip, prominence = lat_strip.max() / prom_scale, \
                               height = emap.max() / height_scale, width = min_peak_width)
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


def combined_figure(cs_xs, cs_ys, streamer_xs, streamer_ys, tomofile, figure_outfile, metric_val):
    """ creates a figure showing electron density peaks and model current sheet locations """

    tomo_lon, tomo_lat, e_density_map = read_tomo(tomofile)
    emap = np.asarray(e_density_map)
    # temporary plotting instructions - for debugging
    plt.imshow(emap, origin='lower', extent=[tomo_lon[0], tomo_lon[-1], tomo_lat[0], tomo_lat[-1]])
    plt.plot(streamer_xs, streamer_ys,'b.', label='Tomo.-based Streamer')
    plt.plot(cs_xs,cs_ys,'rx',label='WSA Neutral Line')
    plt.title('Tomography Result vs. Model')
    plt.ylabel('Latitude (deg)')
    plt.xlabel('Carrington Longitude (deg)')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    plt.annotate("B= " + str(metric_val), (230., -80), color='white', fontsize=18.)
    plt.savefig(figure_outfile)
    plt.close()

    return



def proximity_function(phi, theta, gamma, kappa):
    """ calculates the value of the proximity function at the point phi, theta """

    prefactor = kappa / np.sinh(kappa) / 4. / math.pi
    #print(len(phi))
    gammadotx = angular_separation(phi, theta, gamma[0], gamma[1], \
             cossep=True, degrees=True)
    kent_val = prefactor * math.exp(kappa * gammadotx)
    #kent_val = math.exp(kappa*gammadotx)

    return kent_val



def calc_backbone_metric(cs_xs, cs_ys, streamer_xs, streamer_ys, test=False):
    """ calculates the value of the backbone metric for the given streamer, current sheet locations """

    # streamer_xs are the phi locations of peaks in the electron density, streamer_ys = theta locations
    gammas = np.column_stack((streamer_xs, streamer_ys))
    kappa = 100.

    # find proximity function value at each point in current sheet
    backbone = []
    for phi,theta in zip(cs_xs, cs_ys):
        angdist = []
        for slon,slat in gammas:
            angdist.append(angular_separation(phi, theta, slon, slat, degrees=True))
        gamma_ind = np.argmin(angdist)  # index of minimum angular distance
        backbone.append(proximity_function(phi, theta, gammas[gamma_ind, :], kappa))

    # find proximity function value at each lat,lon in shell
    xs = np.linspace(0,360., num=360)
    ys = np.linspace(-90.,90.,num=181)
    cs_xs2, cs_ys2 = np.meshgrid(xs, ys)
    cs_ys2 = np.reshape(cs_ys2, (360*181))
    cs_xs2 = np.reshape(cs_xs2, (360*181))
    shell_backbone = []
    for phi, theta in gammas:
        angdist = []
        for slon, slat in gammas:
            angdist.append(angular_separation(phi,theta, slon, slat,degrees=True))
        gamma_ind = np.argmin(angdist)
        shell_backbone.append(proximity_function(phi, theta, gammas[gamma_ind,:], kappa))
    

    # average over entire current sheet
    metric = sum(backbone) / len(backbone) / sum(shell_backbone) * len(shell_backbone)
    #metric = np.mean(backbone)

    if test:
        return backbone

    return metric


def backbone_metric(tomofile, bcbfile, altitude, figure_outfile=False):
    """ main module """

    cs_xs, cs_ys = find_neutral_line_coords(bcbfile, altitude)
    streamer_xs, streamer_ys = identify_tomography_peaks(tomofile)
    metric = calc_backbone_metric(cs_xs, cs_ys, streamer_xs, streamer_ys)
    if figure_outfile:
        combined_figure(cs_xs, cs_ys, streamer_xs, streamer_ys, tomofile, figure_outfile, metric)
    return metric



if __name__ == '__main__':

    ARG_PARSER = argparse.ArgumentParser()
    ARG_PARSER.add_argument('-t', '--tomofile', action='store', default='', dest='tomofile')
    ARG_PARSER.add_argument('-a', '--altitude', action='store', default='', dest='altitude', type=float)
    ARG_PARSER.add_argument('-b', '--bcbfile', action='store', default='', dest='bcbfile')
    ARG_PARSER.add_argument('-fo', '--figure_outfile', action='store', default=False)
    ARGS = ARG_PARSER.parse_args()

    #figure_outfile = '/home/sjonesme/Desktop/meeting_image/tomography_currentsheet_comparison.jpg')
    print(backbone_metric(ARGS.tomofile, ARGS.bcbfile, ARGS.altitude, ARGS.figure_outfile))
