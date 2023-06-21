#!/usr/bin/env python3
#
#
# ==========================================================================
# 
# figure_backbone_timeline
# 
# ==========================================================================
# 06/17/2021    sij    created 
#
#-----------------------------------------------------------------------------
""" makes a figure showing the average ADAPT-GONG model performance
    over the course of mid-June to mid-July 2010 """

import path
import numpy as np
import argparse
import os
from astropy.io import fits
from astropy.time import Time
from matplotlib import pyplot as plt
import calc_backbone_metric



def main(listfile, tomofile, outfile, altitude, savefile=False):
    """ main module """

    # read model filenames from listfile
    directory = os.path.dirname(listfile)
    bcbfiles = []
    with open(listfile) as file:
        lines = file.readlines()
        for line in lines:
            bcbfiles.append(line[:-1])

    # check that no files seem to be missing
    if not (len(bcbfiles) % 12 == 0):
        print(f'Number of files listed in {listfile} is not divisible by 12')
        return


    # get all realization 0 files
    real0files = filter(lambda s: 'R000' in s, bcbfiles)

    # calculate backbone metric as a function of time, averaging over all realizations
    mean_backbones = []
    std_backbones = []
    ranges = []
    times = []
    for file in real0files:
        fullfile = f'{directory}/{file}'
        print(f'Running {fullfile}')

        # find mean score
        metric_vals = []
        for realization in range(0,12):    
            thisfile = fullfile.replace("R000", f"R{str(realization).zfill(3)}")
            thismetric = calc_backbone_metric.backbone_metric(tomofile, thisfile, altitude)
            metric_vals.append(thismetric)
        mean_backbones.append(np.mean(metric_vals))
        std_backbones.append(np.std(metric_vals))
        ranges.append([np.mean(metric_vals) - np.min(metric_vals), np.max(metric_vals)-np.mean(metric_vals)])

        # find observation time
        header = fits.getheader(fullfile)
        obstime = header['obstime']
        obstime = obstime.replace('s', '')
        obstime = obstime.replace('m', '')
        obstime = obstime.replace('h', '')
        obstime = obstime.replace('_', ' ')
        obstime = obstime[0:4] + '-' + obstime[5:7] + '-' + obstime[8:]
        times.append(Time(obstime, format = 'iso'))

    # make figure showing performance of mean realization as a function of time
    make_figure(mean_backbones, std_backbones, ranges, times, outfile)

    # save calculation results to file
    if savefile:
        np.savez(savefile, times=times, backbones=backbones_by_realization)




def make_figure(backbone_scores, backbone_stds, ranges, times, outfile):
    """ figure maker """

    # create figure
    fig,ax = plt.subplots()

    # add time series
    ax.errorbar([x.datetime for x in times], backbone_scores, fmt='bo', yerr=np.array(ranges).T.tolist(), \
              linestyle='-', ecolor='red')
    ax.errorbar([x.datetime for x in times], backbone_scores, fmt='bo', yerr=backbone_stds, \
              linestyle='-', ecolor='blue')
    ax.set_title('Model Performance Based on Maps From CR2098')
    ax.set_ylabel('Backbone Metric')
    ax.set_xlabel('Photospheric Map Date')

    # save figure
    fig.autofmt_xdate()
    fig.savefig(outfile)
    print(f'Wrote {outfile}')

    return



if __name__ == '__main__':

    ARG_PARSER = argparse.ArgumentParser()
    ARG_PARSER.add_argument('-t', '--tomofile', action='store', dest='tomofile')
    ARG_PARSER.add_argument('-l', '--listfile', action='store', dest='listfile')
    ARG_PARSER.add_argument('-o', '--outfile', action='store')
    ARG_PARSER.add_argument('-a', '--altitude', action='store', type=float)
    ARG_PARSER.add_argument('-s', '--savefile', action='store', default = False)
    ARGS = ARG_PARSER.parse_args()

    main(ARGS.listfile, ARGS.tomofile, ARGS.outfile, ARGS.altitude, ARGS.savefile)
