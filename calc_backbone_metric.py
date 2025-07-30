#!/usr/bin/env python3
#
#
# ==========================================================================
# 
# calc_backbone_metric.py
# 
# ==========================================================================
# 12/17/2020    sij    copied from compare_streamer_lats.py 
# March 2024    sij    added s-factor comparison capability
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
import pandas as pd
from get_streamerbelt_byhand import rfits_n3d_py
from pathlib import Path


# electron density peak-finding parameters
mindiff = 0.0001 
min_series_max = 10.
prom_scale = 30.
height_scale = 8.
min_peak_width = 8.

# squashing factor peak-finding parameters
mindiff = 0.0001 
min_series_max = 10.
prom_scale = 30.
height_scale = 3
#height_scale = 2
#min_peak_width = 2.
min_peak_width = 1.5

# parameters
#tomodir='/home/sjonesme/stereo_planner/tomography_data/N3D_COR1B_2024/fits/'


def angular_separation(phi0, theta0, phi1, theta1, degrees=False, cossep=False):
    """ calculates distance between two points on a sphere, in radians """

    if degrees:    # convert to radians, if needed
        first_term = math.sin(math.radians(theta0))*math.sin(math.radians(theta1))
        second_term = math.cos(math.radians(theta0)) * math.cos(math.radians(theta1)) * \
                  math.cos(math.radians(phi0-phi1))
    else:
        first_term = math.sin(theta0)*math.sin(theta1)
        second_term = math.cos(theta0) * math.cos(theta1) * math.cos(phi0-phi1)

    dist = np.arccos(first_term + second_term)

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

    # translate to relative strength of radial field component
    Br = np.array(wsa_bfield.br[closest_alt])
    Bphi =  np.array(wsa_bfield.bphi[closest_alt])
    Btheta =  np.array(wsa_bfield.btheta[closest_alt])
    B_tot = norm(np.array([Bphi, Btheta, Br]),axis=0)  # norm is the same in cartesian or spherical coords
    relative_Br = np.divide(Br, B_tot)

    return segment_bfield_slice(relative_Br, wsa_bfield.lons, wsa_bfield.lats)


def segment_bfield_slice(B_slice, lons, lats):

    # use matplotlib contour routine to get contour coordinates
    # Br is the first slice (dimension 2 of four), radial coordinate is dimension 1
    xs = []
    ys = []
    contour_obj = plt.contour(lons, lats, B_slice, levels=[0.0])
    for segs in contour_obj.allsegs[0]:
        xs.extend(segs[:,0])
        ys.extend(segs[:,1])

    return (xs,ys)


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


def locate_ridges(lons, lats, ridgemap):

    ridge_lats=[]
    ridge_lons=[]
    # loop over every column of image
    for lon_ind in range(len(lons)):
        lon_strip = ridgemap[:,lon_ind]
        biggest_diff = np.max(abs(np.diff(lon_strip)))
        if biggest_diff >= mindiff and lon_strip.max() > ridgemap.max() / min_series_max:
            peaks,_ = find_peaks(lon_strip, prominence = lon_strip.max() / prom_scale, \
                               height = ridgemap.max() / height_scale, \
                               width = min_peak_width)
            ridge_lats.extend([lats[peak] for peak in peaks])
            ridge_lons.extend([lons[lon_ind] for peak in peaks])

    # loop over every row of image
    for lat_ind in range(len(lats)):
        wrap_end = round(10. / 360. * len(lons)) #tack on some extra data so we can find peaks near the edges
        lat_strip = np.append(ridgemap[lat_ind, :], ridgemap[lat_ind, : wrap_end])
        biggest_diff = np.max(abs(np.diff(lat_strip)))
        if biggest_diff >= mindiff and lat_strip.max() > ridgemap.max() / min_series_max:
             peaks,_ = find_peaks(lat_strip, prominence = lat_strip.max() / prom_scale, \
                               height = ridgemap.max() / height_scale, width = min_peak_width)
             peaks[peaks >= len(lons)] -= len(lons)
             ridge_lats.extend([lats[lat_ind] for peak in peaks])
             ridge_lons.extend([lons[peak] for peak in peaks])

    # consolidate peaks to grid resolution
    consolidated_lats = []
    consolidated_lons = []
    map_blank = np.zeros((len(lats),len(lons)))
    grid_res = 360. / (len(lons) + 1)
    for slat,slon in zip(ridge_lats, ridge_lons):
        closest_lat = np.argmin([abs(slat - tlat) for tlat in lats])
        closest_lon = np.argmin([abs(slon - tlon) for tlon in lons])
        map_blank[closest_lat, closest_lon] = 1
    nonzeroes = map_blank.nonzero()
    consolidated_lats = [lats[lat_ind] for lat_ind in nonzeroes[0]]
    consolidated_lons = [lons[lon_ind] for lon_ind in nonzeroes[1]]

    return consolidated_lats,consolidated_lons


def identify_tomography_peaks(tomofile):
    """ identify ridges of higher electron density in tomography-based slices """

    # read tomography file
    tomo_lon, tomo_lat, e_density_map = read_tomo(tomofile)
    emap = np.asarray(e_density_map)
    emap[emap < 0.] = 0.

    # get peak locations
    density_lats, density_lons = locate_ridges(tomo_lon, tomo_lat, emap)

    return density_lons, density_lats


def get_emap(tomofile, tomodir):
    # gets electron density slice associated with csv file containing byhand streamer
    #     belt coordinates

    orig_tomofile_basename = Path(tomofile).stem
    alt = float(orig_tomofile_basename[orig_tomofile_basename.rfind('_r')+2:])
    orig_tomofile_basename = orig_tomofile_basename[:orig_tomofile_basename.rfind('_')]
    orig_tomofile = os.path.join(tomodir, orig_tomofile_basename)
    r, th, ph, density = rfits_n3d_py(orig_tomofile)
    closest_rind = np.argmin(np.abs(r-alt))
    tomo_slice = density[closest_rind]  # correct dimension???
    return ph, th, tomo_slice    


def combined_figure(cs_xs, cs_ys, streamer_xs, streamer_ys, tomofile, figure_outfile,
                             metric_val, byhand, tomodir):
    """ creates a figure showing electron density peaks and model current sheet locations """

    print(f'making file name: {figure_outfile}')
    # read density map and grid coordinates
    if byhand:
        tomo_lon, tomo_lat, emap = get_emap(tomofile, tomodir)
    else:
        tomo_lon, tomo_lat, e_density_map = read_tomo(tomofile)
        emap = np.asarray(e_density_map)

    # make figure
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
    backbone_weights = np.cos(np.radians(cs_ys))  # weighting to de-emphasize higher-
                             # latitude points, which are disproportionately represented

    # find proximity function value at each lat,lon in shell
    xs = np.linspace(0,360., num=180)
    ys = np.linspace(-90.,90.,num=91)
    cs_xs2, cs_ys2 = np.meshgrid(xs, ys)
    cs_ys2 = np.reshape(cs_ys2, (180*91))    # note, should be grid-dependent
    cs_xs2 = np.reshape(cs_xs2, (180*91))
    shell_backbone = []
    for phi, theta in zip(cs_xs2, cs_ys2):
        angdist = []
        for slon, slat in gammas:
            angdist.append(angular_separation(phi,theta, slon, slat,degrees=True))
        gamma_ind = np.argmin(angdist)
        shell_backbone.append(proximity_function(phi, theta, gammas[gamma_ind,:], kappa))
    shell_weights = np.cos(np.radians(cs_ys2))

    # average over entire current sheet
    metric = np.average(backbone, weights=backbone_weights) / np.average(shell_backbone,
                             weights=shell_weights)

    if test:
        return backbone

    return metric


def find_wsa_peaks(wsafile, altitude, squashing=False):
    """ find peaks in the squashing factor layer of a standard WSA output file:
           filenames like: wsa_yyyymmddhhmmR<realization>_<observatory>.fits
    """

    # read wsa output file
    slice_of_interest,modelheader = fits.getdata(wsafile, header=True)

    # define grid
    grid = modelheader['GRID']
    nlats = round(180./grid)
    lats = np.linspace(-90. + 0.5*grid, 90. - 0.5*grid, nlats)
    lons = np.linspace(0.5*grid, 360. - 0.5*grid, 2*nlats)

    # find peaks
    if squashing:
        # Squashing factor at RADOUT
        if not altitude == modelheader["RADOUT"]:
            print(f'Warning: User indicated tomographic slice is at {altitude} R_sun.  ' 
                             f'Squashing factor is measured at {modelheader["RADOUT"]} R_sun.')
        slice_of_interest = np.abs(slice_of_interest[9])
        sf_lats, sf_lons = locate_ridges(lons, lats, slice_of_interest)
    else:
        # Coronal field at RADUSER (nT)
        if not altitude == modelheader["RADUSER"]:
            print(f'Warning: User indicated tomographic slice is at {altitude} R_sun.  '
                             f'B-field is measured at {modelheader["RADUSER"]} R_sun.')
        slice_of_interest = slice_of_interest[8]
        sf_lons,sf_lats = segment_bfield_slice(slice_of_interest, lons, lats)

    return sf_lons, sf_lats


def backbone_metric(tomofile, bcb_file=None, squashing_file=None, wsa_file=None, 
                             altitude=2.5, figure_outfile=False, byhand=False,
                             tomodir=None):
    """ calls appropriate calculation routine, depending on desired comparison """

    # get model slice that user wants to compare with tomography results
    if wsa_file:
        # user wants to use B-field slice from WSA output file (wsa_*)
        print(f'Doing neutral line comparison with wsa file: {wsa_file}.')
        cs_xs,cs_ys = find_wsa_peaks(wsa_file, altitude)
    elif squashing_file:
        # user wants to compare squashing factor from WSA output file (wsa_*)
        print(f'Doing squashing factor comparison with  wsa file: {squashing_file}.')
        cs_xs,cs_ys = find_wsa_peaks(squashing_file, altitude, squashing=True)
    elif bcb_file:
        # user wants to find B-field slice in WSA B-field output "bcb" file
        print(f'Doing neutral line comparison with bcb file: {bcb_file}.')
        cs_xs,cs_ys = find_neutral_line_coords(bcb_file, altitude)
    else:
        print("User must specify a bcb_file, squashing_file, or wsa_file")
        return

    # get tomography peaks for comparison
    if byhand:
        streamer_xs, streamer_ys = get_byhand_coords(tomofile)
    else:
        streamer_xs, streamer_ys = identify_tomography_peaks(tomofile)

    # calculate metric value
    metric = calc_backbone_metric(cs_xs, cs_ys, streamer_xs, streamer_ys)

    # optionally, create figure illustrating comparison
    if figure_outfile:
        combined_figure(cs_xs, cs_ys, streamer_xs, streamer_ys, tomofile, figure_outfile,
                             metric, byhand, tomodir)

    return metric


def get_byhand_coords(tomofile):
    # reads longitudes and latitudes from a csv file created by tracing
    #      features by hand

    df = pd.read_csv(tomofile)
    return df.lons, df.lats



if __name__ == '__main__':

    ARG_PARSER = argparse.ArgumentParser()
    ARG_PARSER.add_argument('-t', '--tomofile', action='store', default='', dest='tomofile')
    ARG_PARSER.add_argument('-s', '--squashing_file', action='store', default=None, 
                             dest='squashing_file')
    ARG_PARSER.add_argument('-w', '--wsa_file', action='store', default=None, dest='wsa_file')
    ARG_PARSER.add_argument('-b', '--bcb_file', action='store', default=None, dest='bcb_file')
    ARG_PARSER.add_argument('-a', '--altitude', action='store', default=2.5, dest='altitude',
                             type=float)
    ARG_PARSER.add_argument('-o', '--figure_outfile', action='store', default=False)
    ARG_PARSER.add_argument('-by', '--byhand', action='store_true')
    ARG_PARSER.add_argument('-td', '--tomodir', default=False)
    ARGS = ARG_PARSER.parse_args()

    if ARGS.byhand and ARGS.figure_outfile:
        if not ARGS.tomodir:
            print(colored('If byhand constraints are used, tomodir mush be provided.',
                             'red', attrs={'dark'}))
    else:
        metric = backbone_metric(ARGS.tomofile, bcb_file=ARGS.bcb_file, 
                             squashing_file=ARGS.squashing_file,
                             wsa_file=ARGS.wsa_file, altitude=ARGS.altitude, figure_outfile= 
                             ARGS.figure_outfile, byhand=ARGS.byhand, tomodir=ARGS.tomodir)
        print(f'Metric value: {metric}')
