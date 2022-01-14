#!/usr/bin/env python3
#
#
# ==========================================================================
# 
#
# 
# ==========================================================================
# 06/16/2021    sij    created 
#
#-----------------------------------------------------------------------------
""" read WSA B-field file and define relevant coordinate systems """

from astropy.io import fits
import numpy as np 

class read_wsa_bfield:

    # Assign values to self parameter at initialization
    def __init__(self, bcbfile):
        """ object that contains data and corresponding coordinate axes from WSA bcb* file """

        data,header = fits.getdata(bcbfile, header=True)
        # FORTRAN CODE EXCERPT:  steps = INT((bcube_outer_rad - 1.0d0)/delta_R)
        nsteps = round((header['RADOUT'] - 1.) / header['DELTA_R'])
        self.radii = np.linspace(1.0, header['RADOUT'], nsteps)
        ###   NOT SURE ABOUT THIS - SHOULD IT START AT -90 + GRID/2?
        self.lats = np.linspace(-90, 90, data.shape[2])
        self.lons = np.linspace(0, 360. - header['grid'], round(360./header['grid'])) + header['CARRLONG']
        self.br = data
        self.bcbfile = bcbfile


    def carrington_frame(self):
        """ shift model bfield so that carrington longitude 0 corresponds to index 0 """

        indzero=np.argmin(abs(self.lons-360.))
        if self.lons[indzero] < 360.:
            indzero += 1
        self.br = np.roll(self.br, len(self.lons) - indzero, axis = 3)
        #self.br = np.roll(self.br, indzero, axis = 3)
        shifted_lons = np.roll(self.lons, len(self.lons) - indzero)
        #shifted_lons = np.roll(self.lons, indzero)
        shifted_lons[shifted_lons >= 360.] -= 360.
        self.lons = shifted_lons


# testing code
#bcbfile='/home/sjonesme/wsa/WSA_CAT/WSA/TOMOGRAPHY_COMP/DATA/WSA_OUT/AGONG/ACE/UPDATED/2.0deg/L90R5.0/bcb_201006160200R000.fits'
#wsa_bfield = read_wsa_bfield.read_wsa_bfield(bcbfile)
#br = wsa_bfield.br
#wsa_bfield.carrington_frame()
#newbr = wsa_bfield.br
#fig,ax=plt.subplots(nrows=2, ncols=1)
#ax[0].imshow(br[0,0,:,:], origin='lower')
#ax[1].imshow(newbr[0,0,:,:], origin='lower')
#ax[0].set_title('unrotated')
#ax[1].set_title('rotated')
#fig.savefig('/home/sjonesme/Desktop/meeting_image/carrframe_check.jpg')

