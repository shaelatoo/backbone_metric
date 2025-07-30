#!/usr/bin/env python3
#
#
# ==========================================================================
# 
# get_streamerbelt_byhand.py
# 
# ==========================================================================
# 22 July 2025    sij    created 
#
#-----------------------------------------------------------------------------
"""description of purpose"""

import numpy as np
from matplotlib import pyplot as plt
import argparse
from astropy.io import fits
import os
import pandas as pd


class streamer_belt():
    # much logic blatantly copied from: https://stackoverflow.com/questions/28758079/python-how-to-get-coordinates-on-mouse-click-using-matplotlib-canvas

    def __init__(self, tomo_slice, phis, thetas):
        self.img = tomo_slice
        self.phis = phis
        self.thetas = thetas
        self.xs = []
        self.ys = []
        self.these_xs = []
        self.these_ys = []
        self.ax = None
        self.fig = None

    def __onclick__(self,click):
        self.these_xs.append(round(click.xdata, ndigits=1))
        self.these_ys.append(round(click.ydata, ndigits=1))
        self.ax.plot(round(click.xdata, ndigits=1), round(click.ydata, ndigits=1), 'g*')
        self.fig.canvas.draw()

    def getCoords(self):
        fig,ax = plt.subplots()
        ax.set_title('Click to mark; Close to continue')
        im = ax.imshow(self.img, origin='lower', extent=[self.phis[0], self.phis[-1], \
                self.thetas[0], self.thetas[-1]])
        if len(self.xs) > 0:
            pts = ax.plot(self.xs, self.ys, 'r.')
        cid = fig.canvas.mpl_connect('button_press_event', self.__onclick__)
        self.ax = ax
        self.fig = fig
        plt.show()



def trace_streamer_belt(fitsfilename, outstem, alt=2.5):
    """ used to trace streamer belt features in e- density slices at height

    Args
       fitsfilename: tomography fits file downloaded from SSC with e- density data
       alt: altitude at which to slice tomography data for analysis
    Returns
       streamer_xs: longitudes (degrees) of points that trace out streamer belt
       streamer_ys: latitudes (degrees) of points that trace out streamer belt
    """

    # read tomography file, get slice at alt
    r, thetas, phis, density = rfits_n3d_py(fitsfilename)
    closest_rind = np.argmin(np.abs(r-alt))
    print(f'User-specified altitude: {alt}, tomographic slice altitude: {r[closest_rind]}')
    tomo_slice = density[closest_rind]  # correct dimension???

    # display image and get user clicks
    streamer_structure = streamer_belt(tomo_slice, phis, thetas)
    while True:
        streamer_structure.getCoords()
        choice = 'foo'
        while choice != 'Y' and choice != 'N':
            choice = input("Keep these points? Y or N: ")
        if choice == 'Y':
            streamer_structure.xs = streamer_structure.xs + streamer_structure.these_xs
            streamer_structure.ys = streamer_structure.ys + streamer_structure.these_ys
        streamer_structure.these_xs = []
        streamer_structure.these_ys = []
        quit_or_not = input("All done? Y or N:")
        if quit_or_not == 'Y':
            break

    # save data to file
    outfile = f'{outstem}{os.path.basename(fitsfilename)}_r{alt}.csv'
    df = pd.DataFrame({'lons':streamer_structure.xs, 'lats':streamer_structure.ys})
    df.to_csv(outfile)
    print(f'wrote: {outfile}')

    # return streamer belt coordinates
    return streamer_structure.xs, streamer_structure.ys


def rfits_n3d_py(fitsfilename):
    #based on Tongjiang's IDL code rfits_n3d

    # read data file
    data,header = fits.getdata(fitsfilename, header=True)
    ss = np.shape(data)
    print(ss)
    iph = np.arange(ss[2]) + 1
    ith = np.arange(ss[1]) + 1
    ir = np.arange(ss[0]) + 1

    # data array coordinates
    ph = (iph - header['crpix1']) * header['cdelt1'] + header['crval1']
    th = (ith - header['crpix2']) * header['cdelt2'] + header['crval2']
    r = (ir - header['crpix3']) * header['cdelt3'] + header['crval3']

    return r, th, ph, data



if __name__ == '__main__':

    ARG_PARSER = argparse.ArgumentParser()
    ARG_PARSER.add_argument('fitsfilename', action='store', \
                             help='name of fits file containing tomography data')
    ARG_PARSER.add_argument('outstem', action='store', \
                             help='path or path+prefix of output file')
    ARG_PARSER.add_argument('-a', '--altitude', action='store', type=float, default=2.5, \
                             help='altitude at which to slice electron density data')

    args = ARG_PARSER.parse_args()

    trace_streamer_belt(args.fitsfilename, args.outstem, args.altitude)
