# -*- coding: utf-8 -*-
import numpy as N
import os, re
from setup import MagicSetup
from libmagic import scanDir
import glob
#f2py -c -m greader sub.f90

__author__  = "$Author$"
__date__   = "$Date$"
__version__ = "$Revision$"


try:
    import greader as G
    os.environ['F_UFMTENDIAN'] = 'big'
    os.system('export GFORTRAN_CONVERT_UNIT=big_endian')
    #os.environ['GFORTRAN_CONVERT_UNIT']='big_endian'
    lect = 'f2py'
except ImportError:
    from npfile import *
    lect = 'python'

class MagicGraph(MagicSetup):

    def __init__(self, ivar=None, datadir='.', format='B', quiet=True, 
                 ave=False, tag=None):
        """
        :param format: format of binary output: 'n' (native), 'B' (big endian)
                       or 'l' (little endian)
        :param ave: in case of the average G file G_ave.tag
        :param ivar: the number of the G file
        """
        self.precision = 'f'

        if ave:
            self.name = 'G_ave'
        else:
            self.name = 'G_'

        if tag is not None:
            if ivar is not None:
                file = '%s%i.%s' % (self.name, ivar, tag)
                filename = os.path.join(datadir, file)
            else:
                files = scanDir('%s*%s' % (self.name, tag))
                if len(files) != 0:
                    filename = os.path.join(datadir, files[-1])
                else:
                    print 'No such tag... try again'
                    return

            if os.path.exists('log.%s' % tag):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
        else:
            if ivar is not None:
                files = scanDir('%s%i*' % (self.name, ivar))
                filename = os.path.join(datadir, files[-1])
            else:
                files = scanDir('%s*' % self.name)
                filename = os.path.join(datadir, files[-1])
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists('log.%s' % ending):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % ending)

        if not os.path.exists(filename):
            print 'No such file'
            return

        if lect != 'python':
            G.greader.readg(filename)
            self.nr = G.greader.nr
            self.ntheta = G.greader.nt
            self.npI = G.greader.np
            self.minc = G.greader.minc
            self.time = G.greader.time
            self.ra = G.greader.ra
            self.ek = G.greader.ek
            self.pr = G.greader.pr
            self.prmag = G.greader.prmag
            self.radratio = G.greader.radratio
            self.sigma = G.greader.sigma
            if self.npI == self.ntheta*2:
                self.npI = self.npI/self.minc
            self.nphi = self.npI*self.minc +1
            self.radius = G.greader.radius
            self.colatitude = G.greader.colat
            self.entropy = G.greader.entropy
            self.vr = G.greader.vr
            self.vtheta = G.greader.vt
            self.vphi = G.greader.vp
            if self.prmag != 0:
                self.Br = G.greader.br
                self.Btheta = G.greader.bt
                self.Bphi = G.greader.bp
        else:
            #read data
            inline = npfile(filename, endian=format)

            # read the header
            version = inline.fort_read('|S20')
            runID = inline.fort_read('|S64')[0]
            self.time, n_r_max, n_theta_max, self.npI, n_r_ic_max, minc, nThetaBs, \
                  self.ra, self.ek, self.pr, self.prmag, self.radratio, \
                  self.sigma_ratio = inline.fort_read(self.precision)

            self.nr = int(n_r_max)
            self.ntheta = int(n_theta_max)
            self.npI = int(self.npI)
            self.minc = int(minc)
            if self.npI == self.ntheta*self.minc:
                self.npI = self.npI/self.minc
            self.nr_ic = n_r_ic_max - 1
            self.nThetaBs = int(nThetaBs)

            if not quiet:
                print 'Rayleigh = %.1e, Ekman = %.1e, Prandtl = %.1e' % (self.ra, 
                              self.ek, self.pr)
                print 'nr = %i, nth = %i, nphi = %i' % (self.nr, self.ntheta, 
                              self.npI)

            self.colatitude = inline.fort_read('f')
            self.radius = N.zeros((self.nr), self.precision)

            entropy = N.zeros((self.npI, self.ntheta, self.nr), self.precision)
            vr = N.zeros_like(entropy)
            vtheta = N.zeros_like(entropy)
            vphi = N.zeros_like(entropy)
            if self.prmag != 0:
                Br = N.zeros_like(entropy)
                Btheta = N.zeros_like(entropy)
                Bphi = N.zeros_like(entropy)

            for k in range(n_r_max*nThetaBs):
                # radius and Thetas in this block
                ir, rad, ilat1, ilat2 = inline.fort_read(self.precision)
                ir = int(ir)
                ilat1 = int(ilat1) - 1
                ilat2 = int(ilat2) - 1
                self.radius[ir] = rad
                nth_loc = ilat2 - ilat1 + 1

                if version == 'Graphout_Version_9  ':
                    if self.prmag != 0:
                        data = inline.fort_read(self.precision, shape=(nth_loc,self.npI))
                        entropy[:,ilat1:ilat2+1,ir] = data.T
                        data = inline.fort_read(self.precision, shape=(nth_loc,self.npI))
                        vr[:,ilat1:ilat2+1,ir] = data.T
                        data = inline.fort_read(self.precision, shape=(nth_loc,self.npI))
                        vtheta[:,ilat1:ilat2+1,ir] = data.T
                        data = inline.fort_read(self.precision, shape=(nth_loc,self.npI))
                        vphi[:,ilat1:ilat2+1,ir] = data.T
                        data = inline.fort_read(self.precision, shape=(nth_loc,self.npI))
                        Br[:,ilat1:ilat2+1,ir] = data.T
                        data = inline.fort_read(self.precision, shape=(nth_loc,self.npI))
                        Btheta[:,ilat1:ilat2+1,ir] = data.T
                        data = inline.fort_read(self.precision, shape=(nth_loc,self.npI))
                        Bphi[:,ilat1:ilat2+1,ir] = data.T
                    else:
                        # vectorize
                        # data = inline.fort_read('f',shape=(4*(nth_loc*npI+2)-2)) 
                        data = inline.fort_read(self.precision, shape=(nth_loc,self.npI))
                        entropy[:,ilat1:ilat2+1,ir] = data.T
                        data = inline.fort_read(self.precision, shape=(nth_loc,self.npI))
                        vr[:,ilat1:ilat2+1,ir] = data.T
                        data = inline.fort_read(self.precision, shape=(nth_loc,self.npI))
                        vtheta[:,ilat1:ilat2+1,ir] = data.T
                        data = inline.fort_read(self.precision, shape=(nth_loc,self.npI))
                        vphi[:,ilat1:ilat2+1,ir] = data.T
                else:
                    if self.prmag != 0:
                        # To vectorize, one must also read the end-of-line symbol that
                        # accounts for 2 additionnal floats
                        # For the last line of the chunk, we should not read the last symbol
                        data = inline.fort_read(self.precision,
                                                shape=(7*(self.npI+2)*nth_loc-2))
                        # Add 2 zeros for the last lines
                        data = N.append(data, 0.)
                        data = N.append(data, 0.)
                        data = data.reshape((7, nth_loc, self.npI+2))
                        data = data[:,:,:-2:]
                        entropy[:,ilat1:ilat2+1, ir] = data[0,...].T
                        vr[:, ilat1:ilat2+1, ir] = data[1,...].T
                        vtheta[:, ilat1:ilat2+1, ir] = data[2,...].T
                        vphi[:, ilat1:ilat2+1, ir] = data[3,...].T
                        Br[:, ilat1:ilat2+1, ir] = data[4, ...].T
                        Btheta[:, ilat1:ilat2+1, ir] = data[5, ...].T
                        Bphi[:, ilat1:ilat2+1, ir] = data[6, ...].T
                    else:
                        data = inline.fort_read(self.precision,
                                                shape=(4*(self.npI+2)*nth_loc-2))
                        # Add 2 zeros for the last lines
                        data = N.append(data, 0.)
                        data = N.append(data, 0.)
                        data = data.reshape((4, nth_loc, self.npI+2))
                        data = data[:,:,:-2:]
                        entropy[:,ilat1:ilat2+1, ir] = data[0,...].T
                        vr[:, ilat1:ilat2+1, ir] = data[1,...].T
                        vtheta[:, ilat1:ilat2+1, ir] = data[2,...].T
                        vphi[:, ilat1:ilat2+1, ir] = data[3,...].T

            inline.close()
            # Sorting of data (strange hemispherical way that the magic
            # code uses
            entropy = self.rearangeLat(entropy)
            vr = self.rearangeLat(vr)
            vtheta = self.rearangeLat(vtheta)
            vphi = self.rearangeLat(vphi)
            if self.prmag != 0:
                Br = self.rearangeLat(Br)
                Btheta = self.rearangeLat(Btheta)
                Bphi = self.rearangeLat(Bphi)

            # Normalise r
            self.radius = self.radius/(1.-self.radratio)
            r_cmb = 1./(1.-self.radratio)
            r_icb = r_cmb - 1.

            # Full solution by repeating minc times the structure
            self.nphi = self.minc * self.npI +1
            self.entropy = entropy
            self.vr = vr
            self.vtheta = vtheta
            self.vphi = vphi
            if self.prmag != 0:
                self.Br = Br
                self.Btheta = Btheta
                self.Bphi = Bphi

            """
            for i in range(self.minc):
                self.entropy[i*self.npI:(i+1)*self.npI, ...] = entropy
                self.vr[i*self.npI:(i+1)*self.npI, ...] = vr
                self.vtheta[i*self.npI:(i+1)*self.npI, ...] = vtheta
                self.vphi[i*self.npI:(i+1)*self.npI, ...] = vphi
                if self.mode == 0:
                    self.Br[i*self.npI:(i+1)*self.npI, ...] = Br
                    self.Btheta[i*self.npI:(i+1)*self.npI, ...] = Btheta
                    self.Bphi[i*self.npI:(i+1)*self.npI, ...] = Bphi

            self.entropy[-1, ...] = self.entropy[0, ...]
            self.vr[-1, ...] = self.vr[0, ...]
            self.vtheta[-1, ...] = self.vtheta[0, ...]
            self.vphi[-1, ...] = self.vphi[0, ...]
            if self.mode == 0:
                self.Br[-1, ...] = self.Br[0, ...]
                self.Btheta[-1, ...] = self.Btheta[0, ...]
                self.Bphi[-1, ...] = self.Bphi[0, ...]
            """

    def readOneField(self, inline, lat1, lat2):
        """
        Parse a longitude line...
        Useless now...
        """
        nth_loc = lat2 - lat1 + 1
        # To vectorize, one must also read the end-of-line symbol that
        # accounts for 2 additionnal floats
        # For the last line of the chunk, we should not read the last symbol
        data = inline.fort_read(self.precision, shape=((self.npI+2)*nth_loc-2))
        # Add 2 zeros for the last lines
        data = N.append(data, 0.)
        data = N.append(data, 0.)
        # reshape 
        data = data.reshape((nth_loc, self.npI+2))
        # remove the end-of-line...
        return data[:,:-2:].T

    def rearangeLat(self, field):
        even = field[:, ::2, :] 
        odd = field[:, 1::2, :]
        return N.concatenate((even, odd[:, ::-1, :]), axis=1)

