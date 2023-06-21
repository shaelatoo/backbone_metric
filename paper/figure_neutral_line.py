#!/usr/bin/env python3
#
#
# ==========================================================================
# 
# figure_neutral_line.py
# 
# ==========================================================================
# 09/20/2021    sij    created 
#
#-----------------------------------------------------------------------------
""" create figure showing location of neutral line in coronal model slice """

import read_wsa_bfield
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse
import numpy as np
from matplotlib import pyplot as plt
from calc_backbone_metric import find_neutral_line_coords

def make_figure(bcbfile, altitude, outfile):
    """description of function purpose"""

    # read bcb_* WSA output file, shift to carrington frame
    wsa_bfield = read_wsa_bfield.read_wsa_bfield(bcbfile)
    wsa_bfield.carrington_frame()

    # find altitude slice closest to desired altitude
    closest_alt = np.argmin([abs(altitude - rad) for rad in wsa_bfield.radii])
    Br_alt = (wsa_bfield.br[closest_alt])[0]

    # create figure
    fig,ax = plt.subplots()
    im = ax.imshow(Br_alt/100000., origin='lower', extent=[wsa_bfield.lons[0], wsa_bfield.lons[-1], \
                               wsa_bfield.lats[0], wsa_bfield.lats[-1]])
    ax.set_ylabel('Latitude (deg)')
    ax.set_xlabel('Carrington Longitude (deg)')
    xs,ys = find_neutral_line_coords(bcbfile, altitude)
    ax.plot(xs,ys,'rx')
    ax.set_title('WSA Neutral Line at ' + str(altitude) +  '$R_{sun}$')

    # add colorbar
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=0.05)
    fig.add_axes(ax_cb)
    cb = plt.colorbar(im, cax=ax_cb)
    cb.set_label('Mag. Flux (G)', rotation=90.)
    ax_cb.yaxis.set_label_position("right")
    ax_cb.yaxis.tick_right()
    #ax_cb.yaxis.set_tick_params(labelright=False)

    # write figure to file
    fig.savefig(outfile)
    plt.close(fig)
    print(f'Wrote {outfile}')
 


if __name__ == '__main__':

    ARG_PARSER = argparse.ArgumentParser()
    ARG_PARSER.add_argument('bcbfile', action='store',
                             help='WSA B-field file for which to make figure')
    ARG_PARSER.add_argument('altitude', action='store', default=2.5, type=float,
                             help='Altitude of spherical shell at which to slice model B field')
    ARG_PARSER.add_argument('figure_outfile', action='store', help='Name of file to which to save figure')

    ARGS = ARG_PARSER.parse_args()

    make_figure(ARGS.bcbfile, ARGS.altitude, ARGS.figure_outfile)

    exit(0)
