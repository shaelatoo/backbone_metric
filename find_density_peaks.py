#!/usr/bin/env python3
#
#
# ==========================================================================
# find_density_peaks.py
#
# this was all discarded code from plot_timeline.py - maybe it doesn't work or 
#    is obsolete - it certainly looks incomplete
# ==========================================================================
# 10/15/2020    sij    created 
#
#-----------------------------------------------------------------------------
""" looks for peaks in tomography-based electron density slices """

import numpy as np

def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def multi_gauss_function(x, a1, a2, x01, x02, sigma1, sigma2, a3=False, x03=False, sigma3=False)
    if a3:
        return gauss_function(x, a1, x01, sigma1) + gauss_function(x, a2, x02, sigma2) + \
               gauss_function(x, a3, x03, sigma3)
    return gauss_function(x, a1, x01, sigma1) + gauss_function(x, a2, x02, sigma2)


def fit_multi_gaussian(x, y):

    # trying fitting to a single gaussian curve
    p0 = [max(y), 0., stdev(x)]
    fit_params, foo = curve_fit(gauss_function, x, y, p0 = p0)

    # try fitting to two gaussian curves
    # smooth curve
    # identify regions of high y values
    # actually, try just separating them by a reasonable value, say
    #p0 = [max(y), max(y),
    #return(streamer_lats)


    """description of function purpose"""


if __name__ == '__main__':

    ARG_PARSER = argparse.ArgumentParser()
    ARG_PARSER.add_argument('-f', '--folder', action='store', default='', dest='folder', \
                             help='')


    FOLDER = ARGS.folder.upper()


    exit(0)
