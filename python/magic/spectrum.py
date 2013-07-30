# -*- coding: utf-8 -*-
import os, re
import pylab as P
import numpy as N
from setup import MagicSetup
from libmagic import scanDir

__author__  = "$Author$"
__date__   = "$Date$"
__version__ = "$Revision$"


class MagicSpectrum(MagicSetup):

    def __init__(self, datadir='.', field='e_kin', iplot=True, ispec=None, 
                 ave=False, gather=False, tag=None):
        """
        :param field: the spectrum you want to plot, 'e_kin' for kinetic
                      energy, 'e_mag' for magnetic
        :param iplot: output plot, default is True
        :param ispec: the number of the spectrum you want to plot
        :param ave: in case you want to plot an average spectrum, 
                    then use ave=True
        :param gather: if you want to gather the spectra on the same figure,
                       then use gather=True, default is False
        """
        self.gather = gather

        if field in ('eKin', 'ekin', 'e_kin', 'Ekin', 'E_kin', 'eKinR'):
            if ave:
                self.name = 'kin_spec_ave'
            else:
                self.name = 'kin_spec_'
        elif field in('eMag', 'emag', 'e_mag', 'Emag', 'E_mag', 'eMagR'):
            if ave:
                self.name = 'mag_spec_ave'
            else:
                self.name = 'mag_spec_'

        if tag is not None:
            if ispec is not None:
                file = '%s%i.%s' % (self.name, ispec, tag)
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
            if ispec is not None:
                files = scanDir('%s%i*' % (self.name, ispec))
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

        file = open(filename, 'r')
        if ave is False:
            header = file.readline()
        data = []
        for line in file.readlines():
            st = line.replace('D', 'E')
            data.append(st.split())
        data = N.array(data, dtype='Float64')
        file.close()

        self.index = data[:, 0]
        if self.name == 'kin_spec_ave' or self.name == 'kin_spec_':
            self.ekin_poll = data[:, 1]
            self.ekin_polm = data[:, 2]
            self.ekin_torl = data[:, 3]
            self.ekin_torm = data[:, 4]
        elif self.name == 'mag_spec_ave' or self.name == 'mag_spec_':
            self.emag_poll = data[:, 1]
            self.emag_polm = data[:, 2]
            self.emag_torl = data[:, 3]
            self.emag_torm = data[:, 4]
            self.emagic_poll = data[:, 5]
            self.emagic_polm = data[:, 6]
            self.emagic_torl = data[:, 7]
            self.emagic_torm = data[:, 8]
            self.emagcmb_l = data[:, 9]
            self.emagcmb_m = data[:, 10]
            self.eCMB = data[:, 11]

        if iplot:
            self.plot()

    def plot(self):
        P.rc('figure.subplot', right=0.95, top=0.98, hspace=0.24)
        P.rc('axes', labelsize=20)
        P.rc('xtick', labelsize=14)
        P.rc('ytick', labelsize=14)
        if self.name == 'kin_spec_ave' or self.name == 'kin_spec_':
            if self.gather:
                P.figure()
                P.subplot(211)
                P.loglog(self.index, self.ekin_poll, 'k-', label='poloidal')
                P.loglog(self.index, self.ekin_torl, 'b-', label='toroidal')
                P.xlabel('Degree $\ell$')
                P.ylabel('Kinetic energy')
                P.xlim(self.index.min(), self.index.max())
                P.legend(loc='upper right', frameon=False)
                P.subplot(212)
                P.loglog(self.index[::self.minc]+1, self.ekin_polm[::self.minc],
                         'k-', label='poloidal')
                P.loglog(self.index[::self.minc]+1, self.ekin_torm[::self.minc], 
                         'b-', label='toroidal')
                P.xlabel('Order $m$')
                P.ylabel('Kinetic energy')
                P.xlim(self.index.min(), self.index.max())
                P.legend(loc='upper right', frameon=False)
            else:
                P.figure()
                P.loglog(self.index, self.ekin_poll, 'k-', label='poloidal')
                P.loglog(self.index, self.ekin_torl, 'b-', label='toroidal')
                P.xlabel('Degree $\ell$')
                P.ylabel('Kinetic energy')
                P.xlim(self.index.min(), self.index.max())
                P.legend(loc='upper right', frameon=False)
                P.figure()
                P.loglog(self.index[::self.minc]+1, self.ekin_polm[::self.minc],
                         'k-', label='poloidal')
                P.loglog(self.index[::self.minc]+1, self.ekin_torm[::self.minc],
                         'b-', label='toroidal')
                P.xlabel('$m$ + 1')
                P.ylabel('Kinetic energy')
                P.xlim(self.index.min(), self.index.max())
                P.legend(loc='upper right', frameon=False)
        elif self.name == 'mag_spec_ave' or self.name == 'mag_spec_':
            if self.gather:
                P.figure()
                P.subplot(211)
                P.loglog(self.index, self.emag_poll, 'k-', label='poloidal')
                P.loglog(self.index, self.emag_torl, 'b-', label='toroidal')
                P.loglog(self.index, self.emagcmb_l, 'g-', label='cmb')
                P.xlabel('Degree $\ell$')
                P.ylabel('Magnetic Energy')
                P.xlim(self.index.min(), self.index.max())
                P.legend(loc='upper right', frameon=False)
                P.subplot(212)
                P.loglog(self.index[::self.minc], self.emag_polm[::self.minc],
                         'k-', label='poloidal')
                P.loglog(self.index[::self.minc], self.emag_torm[::self.minc],
                         'b-', label='toroidal')
                P.loglog(self.index[::self.minc], self.emagcmb_m[::self.minc],
                         'g-', label='cmb')
                P.xlabel('Order $m$')
                P.ylabel('Magnetic energy')
                P.xlim(self.index.min(), self.index.max())
                P.legend(loc='upper right', frameon=False)
            else:
                P.figure()
                P.loglog(self.index, self.emag_poll/self.emag_poll.max(), 'k-', label='poloidal')
                P.loglog(self.index, self.emag_torl/self.emag_torl.max(), 'b-', label='toroidal')
                P.loglog(self.index, self.emagcmb_l/self.emagcmb_l.max(), 'g-', label='cmb')
                P.xlabel('Degree $\ell$')
                P.ylabel('Magnetic Energy')
                P.xlim(self.index.min(), self.index.max())
                P.legend(loc='upper right', frameon=False)
                P.figure()
                P.loglog(self.index[::self.minc]+1, self.emag_polm[::self.minc],
                         'k-', label='poloidal')
                P.loglog(self.index[::self.minc]+1, self.emag_torm[::self.minc],
                         'b-', label='toroidal')
                P.loglog(self.index[::self.minc]+1, self.emagcmb_m[::self.minc],
                         'g-', label='cmb')
                P.xlabel('$m$+1')
                P.ylabel('Magnetic energy')
                P.xlim(self.index.min(), self.index.max())
                P.legend(loc='upper right', frameon=False)
