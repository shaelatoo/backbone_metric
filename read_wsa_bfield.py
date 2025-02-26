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
        self.header = header
        self.radii = fits.getdata(bcbfile, ext=1)
        # I'm going to make these as similar as possible to the construction in corona.f90: 
            # theta(J) = PI - 0.5d0*deltaTheta - DBLE(J-1)*deltaTheta
            # phi(I) = ((I*gridRes)-(gridRes/2))*d2r      ! PHI's
            # but working in latitude (not co-latitude) and in degrees
        # not sure why, but the way Fortran writes the data to file and the way Python
        #     reads it from the file seem to be reversed in the theta direction - not
        #     just talking about the way imshow needs the origin set to 0; printing lats[0]
        #     as calculated below yields latitude -90, and the 0th elements of self.br as calculated
        #     below corresponds to the bottom of the image produced by plot_analysis_phfld.py
        gridRes = 360./data.shape[3]
        self.lats = -(90. - gridRes/2. - np.linspace(0, data.shape[2]-1, data.shape[2])*gridRes)
        self.lons = np.linspace(1, data.shape[3], data.shape[3])*gridRes - gridRes/2 + \
                             header['CARRLONG']
        self.br = data[:,0,:,:]
        self.bphi = data[:,1,:,:]
        self.btheta = data[:,2,:,:]
        self.bcbfile = bcbfile


    def carrington_frame(self):
        """ shift model bfield so that carrington longitude 0 corresponds to index 0 """

        if self.header['CARRLONG'] != 0:
            indzero=np.argmin(abs(self.lons-360.))
            if self.lons[indzero] < 360.:
                indzero += 1
            self.br = np.roll(self.br, len(self.lons) - indzero, axis = 2)
            self.btheta = np.roll(self.btheta, len(self.lons) - indzero, axis = 2)
            self.bphi = np.roll(self.bphi, len(self.lons) - indzero, axis = 2)
            shifted_lons = np.roll(self.lons, len(self.lons) - indzero)
            shifted_lons[shifted_lons >= 360.] -= 360.
            self.lons = shifted_lons

    def pad_longitudes(self, npixels):
        """ pad b-field value arrays by npixels in longitudinal direction """

        new_br = []
        new_bphi = []
        new_btheta = []
        for radind,br in enumerate(self.br):
            foo = np.pad(br, pad_width=npixels, mode='wrap')
            new_br.append(foo[npixels:-npixels,:])  # remove bad padding from latitudinal direction
            foo = np.pad(self.bphi[radind], pad_width=npixels, mode='wrap')
            new_bphi.append(foo[npixels:-npixels,:])
            foo = np.pad(self.btheta[radind], pad_width=npixels, mode='wrap')
            new_btheta.append(foo[npixels:-npixels,:])
        self.br = np.array(new_br)
        self.bphi = np.array(new_bphi)
        self.btheta = np.array(new_btheta)
        self.lons = np.pad(self.lons, pad_width=npixels, mode='wrap')

    def pad_latitudes(self, npixels):
        """ pad b-field value arrays by npixels in latitudinal direction """

        print('Warning!  B-field data must be padded first in the longitudinal '
                     'direction, then in latitude.')
        new_br = []
        new_bphi = []
        new_btheta = []
        roll_cnt = int(self.br.shape[2]/2)
        for radind,br in enumerate(self.br):
            foo = np.pad(br, pad_width=npixels, mode='symmetric')
            foo = foo[:,npixels:-npixels]  # remove unneeded longitudinal padding
            foo[:npixels,:] = np.roll(foo[:npixels,:], roll_cnt, axis=1)
            foo[-npixels:,:] = np.roll(foo[-npixels:,:], roll_cnt, axis=1)
            new_br.append(foo)
            foo = np.pad(self.bphi[radind], pad_width=npixels, mode='symmetric')
            foo = foo[:,npixels:-npixels]  # remove unneeded longitudinal padding
            foo[:npixels,:] = np.roll(foo[:npixels,:], roll_cnt, axis=1)
            foo[-npixels:,:] = np.roll(foo[-npixels:,:], roll_cnt, axis=1)
            new_bphi.append(foo)
            foo = np.pad(self.btheta[radind], pad_width=npixels, mode='symmetric')
            foo = foo[:,npixels:-npixels]  # remove unneeded longitudinal padding
            foo[:npixels,:] = np.roll(foo[:npixels,:], roll_cnt, axis=1)
            foo[-npixels:,:] = np.roll(foo[-npixels:,:], roll_cnt, axis=1)
            new_btheta.append(foo)
        self.br = np.array(new_br)
        self.bphi = np.array(new_bphi)
        self.btheta = np.array(new_btheta)
        self.lats = np.pad(self.lats, pad_width=npixels, mode='symmetric')

