import numpy as N
import os
from npfile import *
from setup import MagicSetup
import time as T


class MagicGraph(MagicSetup):

    def __init__(self, ivar=-1, datadir='.', format='B', quiet=True):
        """
        :param format: format of bynary output: 'n' (native), 'B' (big endian)
                       or 'l' (low endian)
        :type format: string
        """
        MagicSetup.__init__(self, datadir=datadir, quiet=True)
        self.precision = 'f'

        if ivar < 0:
            varfile = 'G_1' + '.' + self.tag
        else:
            varfile = 'G_' + str(ivar) + '.' + self.tag

        #read data
        filename = os.path.join(datadir, varfile)
        inline = npfile(filename, endian=format)

        # read the header
        st = inline.fort_read('|S20')
        runID = inline.fort_read('|S64')[0]
        self.time, n_r_max, n_theta_max, self.nphi, n_r_ic_max, minc, nThetaBs, \
              self.ra, self.ek, self.pr, self.prmag, self.radratio, \
              self.sigma_ratio = inline.fort_read(self.precision)

        self.nr = int(n_r_max)
        self.ntheta = int(n_theta_max)
        self.nphi = int(self.nphi)
        self.nr_ic = n_r_ic_max - 1
        self.nThetaBs = int(nThetaBs)
        self.minc = int(minc)

        if not quiet:
            print 'Rayleigh = %.1e, Ekman = %.1e, Prandtl = %.1e' % (self.ra, 
                          self.ek, self.pr)
            print 'nr = %i, nth = %i, nphi = %i' % (self.nr, self.ntheta, 
                          self.nphi)

        self.colatitude = inline.fort_read('f')
        self.radius = N.zeros((self.nr), self.precision)

        self.entropy = N.zeros((self.nphi, self.ntheta, self.nr), self.precision)
        self.vr = N.zeros_like(self.entropy)
        self.vtheta = N.zeros_like(self.entropy)
        self.vphi = N.zeros_like(self.entropy)
        self.Br = N.zeros_like(self.entropy)
        self.Btheta = N.zeros_like(self.entropy)
        self.Bphi = N.zeros_like(self.entropy)

        for k in range(n_r_max*nThetaBs):
            # radius and Thetas in this block
            ir, rad, ilat1, ilat2 = inline.fort_read(self.precision)
            ir = int(ir)
            ilat1 = int(ilat1) - 1
            ilat2 = int(ilat2) - 1
            self.radius[ir] = rad
            self.entropy[:,ilat1:ilat2+1, ir] = self.readOneField(inline, 
                                         ilat1, ilat2)
            self.vr[:, ilat1:ilat2+1, ir] = self.readOneField(inline, 
                                         ilat1, ilat2)
            self.vtheta[:, ilat1:ilat2+1, ir] = self.readOneField(inline, 
                                         ilat1, ilat2)
            self.vphi[:, ilat1:ilat2+1, ir] = self.readOneField(inline, 
                                         ilat1, ilat2)
            self.Br[:, ilat1:ilat2+1, ir] = self.readOneField(inline, 
                                         ilat1, ilat2)
            self.Btheta[:, ilat1:ilat2+1, ir] = self.readOneField(inline, 
                                         ilat1, ilat2)
            self.Bphi[:, ilat1:ilat2+1, ir] = self.readOneField(inline, 
                                         ilat1, ilat2)

        # Sorting of data (strange hemispherical way that the magic
        # code uses
        self.entropy = self.rearangeLat(self.entropy)
        self.vr = self.rearangeLat(self.vr)
        self.vtheta = self.rearangeLat(self.vtheta)
        self.vphi = self.rearangeLat(self.vphi)
        self.Br = self.rearangeLat(self.Br)
        self.Btheta = self.rearangeLat(self.Btheta)
        self.Bphi = self.rearangeLat(self.Bphi)

        # Normalise r
        self.radius = self.radius/(1.-self.radratio)
        r_cmb = 1./(1.-self.radratio)
        r_icb = r_cmb - 1.

        # Full solution by repeating minc times the structure
        if self.minc != 1:
            self.nphi = self.minc * self.nphi +1
            self.entropy = N.zeros((self.nphi, self.ntheta, self.nr), 
                                   self.precision)
            self.vr = N.zeros_like(self.entropy)
            self.vtheta = N.zeros_like(self.entropy)
            self.vphi = N.zeros_like(self.entropy)
            self.Br = N.zeros_like(self.entropy)
            self.Btheta = N.zeros_like(self.entropy)
            self.Bphi = N.zeros_like(self.entropy)

            for i in range(self.minc):


    def readOneField(self, inline, lat1, lat2):
        """
        Parse a longitude line...
        """
        nth_loc = lat2 - lat1 + 1
        data = N.zeros((self.nphi, nth_loc), self.precision)
        #print data.shape
        for k, lat in enumerate(range(nth_loc)):
            dat = inline.fort_read(self.precision, shape=(self.nphi))
            data[:, lat] = dat
        return data

    def rearangeLat(self, field):
        even = field[:, ::2, :] 
        odd = field[:, 1::2, :]
        return N.concatenate((even, odd[:, ::-1, :]), axis=1)

