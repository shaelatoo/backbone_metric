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


def find_neutral_line_coords(bcbfile, altitude):
    """ finds neutral lines in WSA magnetic field model (bcb_* file) slices """

    # read bcb_* WSA output file
    #     FORTRAN CODE EXCERPT:  steps = INT((bcube_outer_rad - 1.0d0)/delta_R)
    data,header = fits.getdata(bcbfile, header=True)
    nsteps = round((header['RADOUT'] - 1.) / header['DELTA_R'])
    radii = np.linspace(1.0, header['RADOUT'], nsteps)
    lats = np.linspace(-90, 90, data.shape[2])     ###   NOT SURE ABOUT THIS - SHOULD IT START AT -90 + GRID/2?
    lons = np.linspace(0, 360. - header['grid'], round(360./header['grid'])) + header['CARRLONG']

    # find altitude slice closest to desired altitude
    closest_alt = np.argmin([abs(altitude - rad) for rad in radii])
    B_alt = data[closest_alt]

    # translate to relative strength of radial field component
    B_tot = norm(B_alt, axis = 0)
    relative_Br_alt = np.divide(B_alt[0], B_tot)

    # shift image so that longitude ~ 0 is on left
    indzero=np.argmin(abs(lons-360.))
    if lons[indzero] < 360.:
        indzero += 1
    shifted_relative_Br_alt = np.roll(relative_Br_alt, len(lons) - indzero, axis = 1)
    shifted_lons = np.roll(lons, len(lons) - indzero)
    shifted_lons[shifted_lons >= 360.] -= 360.

    # use matplotlib contour routine to get contour coordinates
    # Br is the first slice (dimension 2 of four), radial coordinate is dimension 1
    xs = []
    ys = []
    contour_obj = plt.contour(shifted_lons, lats, shifted_relative_Br_alt, levels=[0.0])
    for segs in contour_obj.allsegs[0]:
        xs.extend(segs[:,0])
        ys.extend(segs[:,1])

    return(xs,ys)


def read_tomo(filename):
    """ reads .sav files with slice at fixed radius created by digest_tomography.pro """

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


def angular_separation(phi0, theta0, phi1, theta1, degrees=False, cossep=False):
    """ calculates distance between coordinate locations coords0 and coords1 """

    # convert to radians, if needed
    if degrees:
        first_term = sin(math.radians(theta0))*sin(math.radians(theta1))
        second_term = cos(math.radians(theta0)) * cos(math.radians(theta1)) * \
                  cos(math.radians(phi0-phi1))
    else:
        first_term = sin(theta0)*sin(theta1)
        second_term = cos(theta0) * cos(theta1) * cos(phi0-phi1)
    try:
        dist = np.arccos(first_term + second_term)
    except:
        print(first_term + second_term)

    if cossep:
        return first_term + second_term

    return dist


def proximity_function(phi, theta, gamma, kappa):
    """ calculates the value of the proximity function at the point phi, theta """

    #prefactor = kappa / np.sinh(kappa) / 4. / math.pi
    #print(len(phi))
    gammadotx = angular_separation(phi, theta, gamma[0], gamma[1], cossep=True, degrees=True)
    #kent_val = prefactor * math.exp(kappa * gammadotx)
    kent_val = math.exp(kappa*gammadotx)

    return kent_val



def calc_backbone_metric(cs_xs, cs_ys, streamer_xs, streamer_ys, test=False):
    """ calculates the value of the backbone metric for the given streamer, current sheet locations """

    # streamer_xs are the phi locations of peaks in the electron density, streamer_ys = theta locations
    gammas = np.column_stack((streamer_xs, streamer_ys))
    kappa = 100.

    # find proximity function value at each point in current sheet
    backbone = []
    for phi,theta in zip(cs_xs, cs_ys):
        #print(phi.shape)
        angdist = []
        for slon,slat in gammas:
            angdist.append(angular_separation(phi, theta, slon, slat, degrees=True))
        gamma_ind = np.argmin(angdist)
        val = proximity_function(phi, theta, gammas[gamma_ind, :], kappa)
        backbone.append(val)

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



#def testing(tomofile):
#   """ test procedure created when I was writing the file - may not work anymore """

    #streamer_xs, streamer_ys = identify_tomography_peaks(tomofile)
    #xs = np.linspace(0,360., num=360)
    #ys = np.linspace(-90.,90.,num=181)
    #cs_xs2, cs_ys2 = np.meshgrid(xs, ys)
    #cs_ys2 = np.reshape(cs_ys2, (1, 360*181))
    #cs_xs2 = np.reshape(cs_xs2, (1, 360*181))
    #backbone = calc_backbone_metric(cs_xs2, cs_ys2, streamer_xs, streamer_ys, True)
    #fig = plt.figure(figsize=(10, 6))
    #plt.imshow(np.reshape(backbone,(181,360)), extent=(cs_xs2.min(), cs_xs2.max(), \
    #                    cs_ys2.min(), cs_ys2.max()), origin="lower")
    #plt.colorbar(fraction=0.026, pad=0.04)
    #plt.plot(streamer_xs, streamer_ys, 'b.', label = 'e- density peaks')
    #plt.title('Proximity Function Value, kappa=100.', fontsize=16)
    #plt.ylabel('Latitude', fontsize=14)
    #plt.xlabel('Longitude', fontsize=14)
    #tfilename = os.path.basename(tomofile)
    #tfilestem = os.path.splitext(tfilename)[0]
    #plt.savefig('/home/sjonesme/Desktop/meeting_image/backbone_vals_test_' + \
    #              tfilestem + '.jpg')

    #return


if __name__ == '__main__':

    ARG_PARSER = argparse.ArgumentParser()
    ARG_PARSER.add_argument('-t', '--tomofile', action='store', default='', dest='tomofile')
    ARG_PARSER.add_argument('-a', '--altitude', action='store', default='', dest='altitude', type=float)
    ARG_PARSER.add_argument('-b', '--bcbfile', action='store', default='', dest='bcbfile')
    ARG_PARSER.add_argument('-fo', '--figure_outfile', action='store', default=False)
    ARG_PARSER.add_argument('-test', '--test', action='store_true')
    ARGS = ARG_PARSER.parse_args()
    tomofile = ARGS.tomofile
    bcbfile = ARGS.bcbfile
    altitude = ARGS.altitude
    figure_outfile = ARGS.figure_outfile
    test = ARGS.test

    #figure_outfile = '/home/sjonesme/Desktop/meeting_image/tomography_currentsheet_comparison.jpg')
    if test:
        testing(tomofile)
    else:
        print(backbone_metric(tomofile, bcbfile, altitude, figure_outfile))
        

    exit(0)
