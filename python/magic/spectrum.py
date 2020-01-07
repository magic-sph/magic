# -*- coding: utf-8 -*-
import os, re
import matplotlib.pyplot as plt
import numpy as np
from .log import MagicSetup
from .libmagic import scanDir, fast_read
from .npfile import *
from magic.setup import labTex


class MagicSpectrum(MagicSetup):
    """
    This class can be used to read and display the spectra:

        * Kinetic energy spectra: :ref:`kin_spec_#.TAG <secKinSpecFile>`
        * Magnetic energy spectra: :ref:`mag_spec_#.TAG <secMagSpecFile>`
        * Spectra of the velocity square: :ref:`u2_spec_#.TAG <secu2SpecFile>`

    >>> # display the content of kin_spec_1.tag
    >>> # where tag is the most recent file in the current directory
    >>> sp = MagicSpectrum(field='e_kin', ispec=1)
    >>> # display the content of mag_spec_ave.test on one single figure
    >>> sp = MagicSpectrum(field='e_mag', tag='test', ave=True, gather=True)
    """

    def __init__(self, datadir='.', field='e_kin', iplot=True, ispec=None,
                 ave=False, gather=False, normalize=False, tag=None):
        """
        :param field: the spectrum you want to plot, 'e_kin' for kinetic
                      energy, 'e_mag' for magnetic
        :type field: str
        :param iplot: display the output plot when set to True (default is True)
        :type iplot: bool
        :param ispec: the number of the spectrum you want to plot
        :type ispec: int
        :param tag: file suffix (tag), if not specified the most recent one in
                    the current directory is chosen
        :type tag: str
        :param ave: plot a time-averaged spectrum when set to True
        :type ave: bool
        :param gather: gather the spectra on the same figure when set to True,
                       display one figure per spectrum when set to False,
                       (default is False)
        :type gather: bool
        :param datadir: current working directory
        :type datadir: str
        """
        self.gather = gather
        self.normalize = normalize
        if field in ('eKin', 'ekin', 'e_kin', 'Ekin', 'E_kin', 'eKinR', 'kin'):
            if ave:
                self.name = 'kin_spec_ave'
            else:
                self.name = 'kin_spec_'
        elif field in ('u2', 'usquare', 'u_square', 'uSquare', 'U2'):
            if ave:
                self.name = 'u2_spec_ave'
            else:
                self.name = 'u2_spec_'
        elif field in ('eMag', 'emag', 'e_mag', 'Emag', 'E_mag', 'eMagR', 'mag'):
            if ave:
                self.name = 'mag_spec_ave'
            else:
                self.name = 'mag_spec_'
        elif field in ('dtVrms'):
            self.name = 'dtVrms_spec'
        elif field in ('T','temperature','S','entropy'):
            self.name = 'T_spec_'


        if tag is not None:
            if ispec is not None:
                file = '%s%i.%s' % (self.name, ispec, tag)
                filename = os.path.join(datadir, file)
            else:
                pattern = os.path.join(datadir, '%s*%s' % (self.name, tag))
                files = scanDir(pattern)
                if len(files) != 0:
                    filename = files[-1]
                else:
                    print('No such tag... try again')
                    return

            if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
        else:
            if ispec is not None:
                pattern = os.path.join(datadir, '%s%i*' % (self.name, ispec))
                files = scanDir(pattern)
                filename = files[-1]
            else:
                pattern = os.path.join(datadir, '%s*' % self.name)
                files = scanDir(pattern)
                filename = files[-1]
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists(os.path.join(datadir, 'log.%s' % ending)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % ending)

        if not os.path.exists(filename):
            print('No such file')
            return

        if ave is False and self.name != 'dtVrms_spec':
            data = fast_read(filename, skiplines=1)
        else:
            data = fast_read(filename)

        self.index = data[:, 0]
        if self.name == 'kin_spec_':
            self.ekin_poll = data[:, 1]
            self.ekin_polm = data[:, 2]
            self.ekin_torl = data[:, 3]
            self.ekin_torm = data[:, 4]
        elif self.name == 'kin_spec_ave':
            self.ekin_poll = data[:, 1]
            self.ekin_polm = data[:, 2]
            self.ekin_torl = data[:, 3]
            self.ekin_torm = data[:, 4]
            self.ekin_poll_SD = data[:, 5]
            self.ekin_polm_SD = data[:, 6]
            self.ekin_torl_SD = data[:, 7]
            self.ekin_torm_SD = data[:, 8]
        elif self.name == 'u2_spec_ave':
            self.ekin_poll = data[:, 1]
            self.ekin_polm = data[:, 2]
            self.ekin_torl = data[:, 3]
            self.ekin_torm = data[:, 4]
        elif self.name == 'u2_spec_ave' or self.name == 'u2_spec_':
            self.ekin_poll = data[:, 1]
            self.ekin_polm = data[:, 2]
            self.ekin_torl = data[:, 3]
            self.ekin_torm = data[:, 4]
            self.ekin_poll_SD = data[:, 5]
            self.ekin_polm_SD = data[:, 6]
            self.ekin_torl_SD = data[:, 7]
            self.ekin_torm_SD = data[:, 8]
        elif self.name == 'mag_spec_':
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
        elif self.name == 'mag_spec_ave':
            self.emag_poll = data[:, 1]
            self.emag_polm = data[:, 2]
            self.emag_torl = data[:, 3]
            self.emag_torm = data[:, 4]
            self.emagcmb_l = data[:, 5]
            self.emagcmb_m = data[:, 6]
            self.emag_poll_SD = data[:, 7]
            self.emag_polm_SD = data[:, 8]
            self.emag_torl_SD = data[:, 9]
            self.emag_torm_SD = data[:, 10]
            self.emagcmb_l_SD = data[:, 11]
            self.emagcmb_m_SD = data[:, 12]
        elif self.name == 'dtVrms_spec':
            self.InerRms = data[:, 1]
            self.CorRms = data[:, 2]
            self.LFRms = data[:, 3]
            self.AdvRms = data[:, 4]
            self.DifRms = data[:, 5]
            self.BuoRms = data[:, 6]

            if data.shape[1] == 27:
                self.PreRms = data[:, 7]
                self.geos = data[:, 8] # geostrophic balance
                self.mageos = data[:, 9] # magnetostrophic balance
                self.arcMag = data[:, 10] # Pressure/Coriolis/Lorentz/Buoyancy
                self.corLor = data[:, 11] # Coriolis/Lorentz
                self.preLor = data[:, 12] # Pressure/Lorentz
                self.cia = data[:, 13] # Coriolis/Inertia/Archimedean
                self.InerRms_SD = data[:, 14]
                self.CorRms_SD = data[:, 15]
                self.LFRms_SD = data[:, 16]
                self.AdvRms_SD = data[:, 17]
                self.DifRms_SD = data[:, 18]
                self.BuoRms_SD = data[:, 19]
                self.PreRms_SD = data[:, 20]
                self.geos_SD = data[:, 21]
                self.mageos_SD = data[:, 22]
                self.arc_SD = data[:, 23]
                self.corLor_SD = data[:, 24]
                self.preLor_SD = data[:, 25]
                self.cia_SD = data[:, 26]
                self.arc = np.zeros_like(self.cia)
                self.arc_SD = np.zeros_like(self.cia)
                self.ChemRms = np.zeros_like(self.cia)
                self.ChemRms_SD = np.zeros_like(self.cia)
            elif data.shape[1] == 29:
                self.PreRms = data[:, 7]
                self.geos = data[:, 8] # geostrophic balance
                self.mageos = data[:, 9] # magnetostrophic balance
                self.arc = data[:, 10] # Pressure/Coriolis/Buoyancy
                self.arcMag = data[:, 11] # Pressure/Coriolis/Lorentz/Buoyancy
                self.corLor = data[:, 12] # Coriolis/Lorentz
                self.preLor = data[:, 13] # Pressure/Lorentz
                self.cia = data[:, 14] # Coriolis/Inertia/Archimedean
                self.InerRms_SD = data[:, 15]
                self.CorRms_SD = data[:, 16]
                self.LFRms_SD = data[:, 17]
                self.AdvRms_SD = data[:, 18]
                self.DifRms_SD = data[:, 19]
                self.BuoRms_SD = data[:, 20]
                self.PreRms_SD = data[:, 21]
                self.geos_SD = data[:, 22]
                self.mageos_SD = data[:, 23]
                self.arc_SD = data[:, 24]
                self.arcMag_SD = data[:, 25]
                self.corLor_SD = data[:, 26]
                self.preLor_SD = data[:, 27]
                self.cia_SD = data[:, 28]
                self.ChemRms = np.zeros_like(self.cia)
                self.ChemRms_SD = np.zeros_like(self.cia)
            else:
                self.ChemRms = data[:,7]
                self.PreRms = data[:, 8]
                self.geos = data[:, 9] # geostrophic balance
                self.mageos = data[:,10] # magnetostrophic balance
                self.arc = data[:, 11] # Pressure/Coriolis/Buoyancy
                self.arcMag = data[:, 12] # Pressure/Coriolis/Lorentz/Buoyancy
                self.corLor = data[:, 13] # Coriolis/Lorentz
                self.preLor = data[:, 14] # Pressure/Lorentz
                self.cia = data[:, 15] # Coriolis/Inertia/Archimedean
                self.InerRms_SD = data[:, 16]
                self.CorRms_SD = data[:, 17]
                self.LFRms_SD = data[:, 18]
                self.AdvRms_SD = data[:, 19]
                self.DifRms_SD = data[:, 20]
                self.BuoRms_SD = data[:, 21]
                self.ChemRms_SD = data[:, 22]
                self.PreRms_SD = data[:, 23]
                self.geos_SD = data[:, 24]
                self.mageos_SD = data[:, 25]
                self.arc_SD = data[:, 26]
                self.arcMag_SD = data[:, 27]
                self.corLor_SD = data[:, 28]
                self.preLor_SD = data[:, 29]
                self.cia_SD = data[:, 30]
        elif self.name == 'T_spec_':
            self.T_l = data[:,1]
            self.T_m = data[:,2]
            self.T_icb_l = data[:,3]
            self.T_icb_m = data[:,4]
            self.dT_icb_l = data[:,5]
            self.dT_icb_m = data[:,6]

        if iplot:
            self.plot()

    def plot(self):
        """
        Plotting function
        """
        if self.name == 'kin_spec_ave' or self.name == 'kin_spec_':
            if self.gather:
                fig = plt.figure()
                ax = fig.add_subplot(211)
                ax.loglog(self.index, self.ekin_poll, label='poloidal')
                ax.loglog(self.index, self.ekin_torl, label='toroidal')
                if labTex:
                    ax.set_xlabel('Degree $\ell$')
                else:
                    ax.set_xlabel('Degree l')
                ax.set_ylabel('Kinetic energy')
                ax.set_xlim(self.index.min(), self.index.max())
                ax.legend(loc='upper right', frameon=False)
                ax = fig.add_subplot(212)
                ax.loglog(self.index[::self.minc]+1, self.ekin_polm[::self.minc],
                          label='poloidal')
                ax.loglog(self.index[::self.minc]+1, self.ekin_torm[::self.minc],
                          label='toroidal')
                if labTex:
                    ax.set_xlabel('Order $m+1$')
                else:
                    ax.set_xlabel('m+1')
                ax.set_ylabel('Kinetic energy')
                ax.set_xlim(self.index.min(), self.index.max())
                ax.legend(loc='upper right', frameon=False)
            else:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                if self.normalize:
                    y = self.ekin_poll+self.ekin_torl
                    ax.loglog(self.index, y/y.max(),)
                else:
                    ax.loglog(self.index, self.ekin_poll, label='poloidal')
                    ax.loglog(self.index, self.ekin_torl, label='toroidal')
                if labTex:
                    ax.set_xlabel('Degree $\ell$')
                else:
                    ax.set_xlabel('Degree l')
                ax.set_ylabel('Kinetic energy')
                ax.set_xlim(self.index.min(), self.index.max())
                ax.legend(loc='upper right', frameon=False)

                fig = plt.figure()
                ax = fig.add_subplot(111)
                if self.normalize:
                    y = self.ekin_polm[::self.minc]+self.ekin_torm[::self.minc]
                    ax.loglog(self.index[::self.minc]+1, y/y.max())
                else:
                    ax.loglog(self.index[::self.minc]+1, self.ekin_polm[::self.minc],
                              label='poloidal')
                    ax.loglog(self.index[::self.minc]+1, self.ekin_torm[::self.minc],
                              label='toroidal')
                if labTex:
                    ax.set_xlabel('$m$ + 1')
                else:
                    ax.set_xlabel('m + 1')
                ax.set_ylabel('Kinetic energy')
                ax.set_xlim(self.index.min(), self.index.max())
                ax.legend(loc='upper right', frameon=False)
        elif self.name == 'u2_spec_ave' or self.name == 'u2_spec_':
            if self.gather:
                fig = plt.figure()
                ax = fig.add_subplot(211)
                ax.loglog(self.index, self.ekin_poll, label='poloidal')
                ax.loglog(self.index, self.ekin_torl, label='toroidal')
                if labTex:
                    ax.set_xlabel('Degree $\ell$')
                    ax.set_ylabel(r'${\cal U}^2$')
                else:
                    ax.set_xlabel('Degree l')
                    ax.set_ylabel('Velocity square')
                ax.set_xlim(self.index.min(), self.index.max())
                ax.legend(loc='upper right', frameon=False)
                ax = fig.add_subplot(212)
                ax.loglog(self.index[::self.minc]+1, self.ekin_polm[::self.minc],
                          label='poloidal')
                ax.loglog(self.index[::self.minc]+1, self.ekin_torm[::self.minc],
                          label='toroidal')
                if labTex:
                    ax.set_xlabel('Order $m+1$')
                    ax.set_ylabel(r'${\cal U}^2$')
                else:
                    ax.set_xlabel('m+1')
                    ax.set_ylabel('Velocity square')
                ax.set_xlim(self.index.min(), self.index.max())
                ax.legend(loc='upper right', frameon=False)
            else:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.loglog(self.index, self.ekin_poll, label='poloidal')
                ax.loglog(self.index, self.ekin_torl, label='toroidal')
                if labTex:
                    ax.set_xlabel('Degree $\ell$')
                    ax.set_ylabel(r'${\cal U}^2$')
                else:
                    ax.set_xlabel('Degree l')
                    ax.set_ylabel('Velocity square')
                ax.set_xlim(self.index.min(), self.index.max())
                ax.legend(loc='upper right', frameon=False)

                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.loglog(self.index[::self.minc]+1, self.ekin_polm[::self.minc],
                          label='poloidal')
                ax.loglog(self.index[::self.minc]+1, self.ekin_torm[::self.minc],
                          label='toroidal')
                if labTex:
                    ax.set_xlabel(r'$m + 1$')
                    ax.set_ylabel(r'${\cal U}^2$')
                else:
                    ax.set_xlabel('m + 1')
                    ax.set_ylabel('Velocity square')
                ax.set_xlim(self.index.min(), self.index.max())
                ax.legend(loc='upper right', frameon=False)
        elif self.name == 'mag_spec_ave' or self.name == 'mag_spec_':
            if self.gather:
                fig = plt.figure()
                ax = fig.add_subplot(211)
                ax.loglog(self.index, self.emag_poll, label='poloidal')
                ax.loglog(self.index, self.emag_torl, label='toroidal')
                ax.loglog(self.index, self.emagcmb_l, label='cmb')
                if labTex:
                    ax.set_xlabel('Degree $\ell$')
                else:
                    ax.set_xlabel('Degree l')
                ax.set_ylabel('Magnetic Energy')
                ax.set_xlim(self.index.min(), self.index.max())
                ax.legend(loc='upper right', frameon=False)

                ax = fig.add_subplot(212)
                ax.loglog(self.index[::self.minc], self.emag_polm[::self.minc],
                          label='poloidal')
                ax.loglog(self.index[::self.minc], self.emag_torm[::self.minc],
                          label='toroidal')
                ax.loglog(self.index[::self.minc], self.emagcmb_m[::self.minc],
                          label='cmb')
                if labTex:
                    ax.set_xlabel('Order $m$')
                else:
                    ax.set_xlabel('Order m')
                ax.set_ylabel('Magnetic energy')
                ax.set_xlim(self.index.min(), self.index.max())
                ax.legend(loc='upper right', frameon=False)
            else:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.loglog(self.index, self.emag_poll/self.emag_poll.max(),
                          label='poloidal')
                ax.loglog(self.index, self.emag_torl/self.emag_torl.max(),
                          label='toroidal')
                ax.loglog(self.index, self.emagcmb_l/self.emagcmb_l.max(),
                          label='cmb')
                if labTex:
                    ax.set_xlabel('Degree $\ell$')
                else:
                    ax.set_xlabel('Degree l')
                ax.set_ylabel('Magnetic Energy')
                ax.set_xlim(self.index.min(), self.index.max())
                ax.legend(loc='upper right', frameon=False)

                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.loglog(self.index[::self.minc]+1,
                          self.emag_polm[::self.minc]/self.emag_polm.max(),
                          label='poloidal')
                ax.loglog(self.index[::self.minc]+1,
                          self.emag_torm[::self.minc]/self.emag_torm.max(),
                          label='toroidal')
                ax.loglog(self.index[::self.minc]+1,
                          self.emagcmb_m[::self.minc]/self.emagcmb_m.max(),
                          label='cmb')
                if labTex:
                    ax.set_xlabel('$m$+1')
                else:
                    ax.set_xlabel('m+1')
                ax.set_ylabel('Magnetic energy')
                ax.set_xlim(self.index.min(), self.index.max())
                ax.legend(loc='upper right', frameon=False)

        elif self.name == 'dtVrms_spec':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.loglog(self.index, self.CorRms, label='Coriolis')
            ax.loglog(self.index, self.PreRms, label='Pressure')
            if self.LFRms.max() > 1e-11:
                ax.loglog(self.index, self.LFRms, label='Lorentz')
            ax.loglog(self.index, self.BuoRms + self.ChemRms, label='Buoyancy')
            ax.loglog(self.index, self.InerRms, label='Inertia')
            ax.loglog(self.index, self.DifRms, label='Viscosity')
            ax.loglog(self.index, self.geos, label='Coriolis-Pressure')
            ax.loglog(self.index, self.arcMag, ls='--',
                      label='Coriolis-Pressure-Buoyancy-Lorentz')

            if labTex:
                ax.set_xlabel('$\ell+1$')
            else:
                ax.set_xlabel('l+1')
            ax.set_ylabel('RMS forces')
            ax.set_xlim(self.index.min(), self.index.max())
            ax.legend(loc='lower right', frameon=False, ncol=2)

        elif self.name == 'T_spec_':
            if self.gather:
                fig = plt.figure()
                ax = fig.add_subplot(211)
                ax.loglog(self.index, self.T_l/self.T_l.max(),)
                ax.loglog(self.index, self.T_icb_l/self.T_icb_l.max(),label='ICB')
                if labTex:
                    ax.set_xlabel('$\ell$')
                else:
                    ax.set_xlabel('l')
                ax.set_ylabel('degree')
                ax.legend()

                ax = fig.add_subplot(212)
                ax.loglog(self.index[::self.minc]+1,
                          self.T_m[::self.minc]/self.T_m[::self.minc].max(),)
                ax.loglog(self.index[::self.minc]+1,
                          self.T_icb_m[::self.minc]/self.T_icb_m[::self.minc].max(),
                          label='ICB')

                if labTex:
                    ax.set_xlabel('Order $m+1$')
                else:
                    ax.set_xlabel('m+1')
                ax.set_ylabel('order')
                ax.legend()
                fig.tight_layout()
            else:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.loglog(self.index, self.T_l)
                ax.loglog(self.index, self.T_icb_l)
                if labTex:
                    ax.set_xlabel('$\ell$')
                else:
                    ax.set_xlabel('l')

                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.loglog(self.index[::self.minc]+1,
                          self.T_m[::self.minc],color='g')
                if labTex:
                    ax.set_xlabel('Order $m+1$')
                else:
                    ax.set_xlabel('m+1')


class MagicSpectrum2D(MagicSetup):
    """
    This class can be used to read and display 2-D spectra in the :math:`(r,\ell)`
    and in the :math:`(r,m)` planes

        * Kinetic energy spectra: :ref:`2D_kin_spec_#.TAG <sec2DSpectra>`
        * Magnetic energy spectra: :ref:`2D_mag_spec_#.TAG <sec2DSpectra>`

    >>> # display the content of 2D_kin_spec_1.tag
    >>> # where tag is the most recent file in the current directory
    >>> sp = MagicSpectrum2D(field='e_kin', ispec=1, levels=17, cm='seismic')
    >>> # display the content of 2D_mag_spec_3.test
    >>> sp = MagicSpectrum2D(field='e_mag', tag='test', ispec=3)
    """

    def __init__(self, datadir='.', field='e_mag', iplot=False, ispec=None,
                 tag=None, cm='jet', levels=33, precision=np.float64,
                 ave=False):
        """
        :param field: the spectrum you want to plot, 'e_kin' for kinetic
                      energy, 'e_mag' for magnetic
        :type field: str
        :param iplot: display the output when set to True (default is True)
        :type iplot: bool
        :param ispec: the number of the spectrum you want to plot
        :type ispec: int
        :param tag: file suffix (tag=, if not specified the most recent one
                    in the current directory is chosen
        :type tag: str
        :param cm: name of the colormap (default='jet')
        :type cm: str
        :param levels: number of contour levels (default 33)
        :type levels: int
        :param precision: single or double precision
        :type precision: str
        :param datadir: current working directory
        :type datadir: str
        :param ave: plot a time-averaged spectrum when set to True
        :type ave: bool
        """

        if field in ('eKin', 'ekin', 'e_kin', 'Ekin', 'E_kin', 'eKinR', 'kin'):
            if ave:
                self.name = '2D_kin_spec_ave'
            else:
                self.name = '2D_kin_spec_'
        elif field in('eMag', 'emag', 'e_mag', 'Emag', 'E_mag', 'eMagR', 'mag'):
            if ave:
                self.name = '2D_mag_spec_ave'
            else:
                self.name = '2D_mag_spec_'
        elif field in ('dtVrms'):
            self.name = '2D_dtVrms_spec'

        if ave:
            self.version = 'ave'
        else:
            self.version = 'snap'


        if tag is not None:
            if ispec is not None:
                file = '%s%i.%s' % (self.name, ispec, tag)
                filename = os.path.join(datadir, file)
            else:
                pattern = os.path.join(datadir, '%s*%s' % (self.name, tag))
                files = scanDir(pattern)
                if len(files) != 0:
                    filename = files[-1]
                else:
                    print('No such tag... try again')
                    return

            if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
        else:
            if ispec is not None:
                pattern = os.path.join(datadir, '%s%i*' % (self.name, ispec))
                files = scanDir(pattern)
                filename = files[-1]
            else:
                pattern = os.path.join(datadir, '%s*' % self.name)
                files = scanDir(pattern)
                filename = files[-1]
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists(os.path.join(datadir, 'log.%s' % ending)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % ending)

        if not os.path.exists(filename):
            print('No such file')
            return

        f = npfile(filename, endian='B')

        if self.name == '2D_dtVrms_spec':
            l_one = f.fort_read('i4')
            if len(l_one) == 1:
                self.version = l_one[0]
                self.n_r_max, self.l_max = f.fort_read('i4')
                self.rad = f.fort_read(precision, shape=(self.n_r_max))
                self.Cor_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.Adv_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.LF_r_l  = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.Buo_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.Chem_r_l = f.fort_read(precision,
                                            shape=(self.n_r_max, self.l_max+1))
                self.Pre_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.Dif_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.Iner_r_l = f.fort_read(precision,
                                            shape=(self.n_r_max, self.l_max+1))
                self.Geo_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.Mag_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.Arc_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.ArcMag_r_l = f.fort_read(precision,
                                              shape=(self.n_r_max, self.l_max+1))
                self.CIA_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.CLF_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.PLF_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
            else:
                self.n_r_max, self.l_max = l_one
                self.rad = f.fort_read(precision, shape=(self.n_r_max))
                self.Cor_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.Adv_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.LF_r_l  = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.Buo_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.Chem_r_l = np.zeros_like(self.Buo_r_l)
                self.Pre_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.Dif_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.Iner_r_l = f.fort_read(precision,
                                            shape=(self.n_r_max, self.l_max+1))
                self.Geo_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.Mag_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.Arc_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.ArcMag_r_l = f.fort_read(precision,
                                              shape=(self.n_r_max, self.l_max+1))
                self.CIA_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.CLF_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))
                self.PLF_r_l = f.fort_read(precision,
                                           shape=(self.n_r_max, self.l_max+1))

        else:
            if self.version == 'snap':
                if precision == np.float64:
                    out = f.fort_read('f8,3i4')[0]
                else:
                    out = f.fort_read('f4,3i4')[0]
                self.time = out[0]
                self.n_r_max, self.l_max, self.minc = out[1]
            elif self.version == 'ave':
                self.n_r_max, self.l_max, self.minc = f.fort_read('3i4')[0]
                self.time = -1.
            self.rad = f.fort_read(precision, shape=(self.n_r_max))
            self.e_pol_l = f.fort_read(precision, shape=(self.l_max, self.n_r_max))
            self.e_pol_m = f.fort_read(precision, shape=(self.l_max+1, self.n_r_max))
            self.e_tor_l = f.fort_read(precision, shape=(self.l_max, self.n_r_max))
            self.e_tor_m = f.fort_read(precision, shape=(self.l_max+1, self.n_r_max))

        self.ell = np.arange(self.l_max+1)
        f.close()

        if iplot:
            self.plot(levels, cm)


    def plot(self, levels, cm, cut=1.):
        """
        Plotting function

        :param levels: number of contour levels
        :type levels: int
        :param cm: name of the colormap
        :param cut: adjust the contour maximum to max(abs(data))*cut
        :type cut: float
        """
        if self.name == '2D_dtVrms_spec':
            vmax = np.log10(cut*self.Geo_r_l).max()
            vmin = vmax-4
            levs = np.linspace(vmin, vmax, levels)
            fig0 = plt.figure()
            ax0 = fig0.add_subplot(111)
            im = ax0.contourf(self.rad, self.ell[1:],
                              np.log10(self.Geo_r_l[:, 1:].T),
                              levs, cmap=plt.get_cmap(cm), extend='both')
            if labTex:
                ax0.set_ylabel('Degree $\ell$')
                ax0.set_xlabel('Radius $r$')
            else:
                ax0.set_ylabel('Degree l')
                ax0.set_xlabel('Radius')
            ax0.set_yscale('log')
            ax0.set_title('Coriolis - Pressure')
            plt.ylim([1,self.l_max])
            fig0.colorbar(im)

            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111, sharex=ax0, sharey=ax0)
            im = ax1.contourf(self.rad, self.ell[1:],
                              np.log10((self.Buo_r_l[:, 1:]+self.Chem_r_l[:, 1:]).T),
                              levs, cmap=plt.get_cmap(cm), extend='both')
            if labTex:
                ax1.set_ylabel('Degree $\ell$')
                ax1.set_xlabel('Radius $r$')
            else:
                ax1.set_ylabel('Degree l')
                ax1.set_xlabel('Radius')
            ax1.set_yscale('log')
            ax1.set_title('Buoyancy')
            plt.ylim([1,self.l_max])
            fig1.colorbar(im)

            if abs(self.LF_r_l).max() > 0:
                fig2 = plt.figure()
                ax2 = fig2.add_subplot(111,sharex=ax0, sharey=ax0)
                im = ax2.contourf(self.rad, self.ell[1:],
                                  np.log10(self.LF_r_l[:, 1:].T),
                                  levs, cmap=plt.get_cmap(cm), extend='both')
                if labTex:
                    ax2.set_ylabel('Degree $\ell$')
                    ax2.set_xlabel('Radius $r$')
                else:
                    ax2.set_ylabel('Degree l')
                    ax2.set_xlabel('Radius')
                ax2.set_yscale('log')
                ax2.set_title('Lorentz force')
                plt.ylim([1,self.l_max])
                fig2.colorbar(im)

            fig3 = plt.figure()
            ax3 = fig3.add_subplot(111, sharex=ax0, sharey=ax0)
            im = ax3.contourf(self.rad, self.ell[1:],
                              np.log10(self.Iner_r_l[:, 1:].T),
                              levs, cmap=plt.get_cmap(cm), extend='both')
            if labTex:
                ax3.set_ylabel('Degree $\ell$')
                ax3.set_xlabel('Radius $r$')
            else:
                ax3.set_ylabel('Degree l')
                ax3.set_xlabel('Radius')
            ax3.set_yscale('log')
            ax3.set_title('Inertia')
            plt.ylim([1,self.l_max])
            fig3.colorbar(im)

            fig4 = plt.figure()
            ax4 = fig4.add_subplot(111, sharex=ax0, sharey=ax0)
            im = ax4.contourf(self.rad, self.ell[1:],
                              np.log10(self.Dif_r_l[:, 1:].T),
                              levs, cmap=plt.get_cmap(cm), extend='both')
            if labTex:
                ax4.set_ylabel('Degree $\ell$')
                ax4.set_xlabel('Radius $r$')
            else:
                ax4.set_ylabel('Degree l')
                ax4.set_xlabel('Radius')
            ax4.set_yscale('log')
            ax4.set_title('Viscosity')
            plt.ylim([1,self.l_max])
            fig4.colorbar(im)
        else:
            fig0 = plt.figure()
            ax0 = fig0.add_subplot(111)
            vmax = np.log10(cut*self.e_pol_l).max()
            vmin = vmax-7
            levs = np.linspace(vmin, vmax, levels)
            im = ax0.contourf(self.rad, self.ell[1:], np.log10(self.e_pol_l),
                              levs, cmap=plt.get_cmap(cm), extend='both')
            fig0.colorbar(im)
            if labTex:
                ax0.set_ylabel('Degree $\ell$')
                ax0.set_xlabel('Radius $r$')
            else:
                ax0.set_ylabel('Degree l')
                ax0.set_xlabel('Radius')
            ax0.set_yscale('log')
            ax0.set_title('E pol')

            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            vmax = np.log10(self.e_tor_l).max()
            vmin = vmax-14
            levs = np.linspace(vmin, vmax, levels)
            im = ax1.contourf(self.rad, self.ell[1:], np.log10(self.e_tor_l),
                              levs, cmap=plt.get_cmap(cm), extend='both')
            fig1.colorbar(im)
            if labTex:
                ax1.set_ylabel('Degree $\ell$')
                ax1.set_xlabel('Radius $r$')
            else:
                ax1.set_ylabel('Degree l')
                ax1.set_xlabel('Radius')
            ax1.set_yscale('log')
            ax1.set_title('E tor')

            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111)
            vmax = np.log10(self.e_pol_m).max()
            vmin = vmax-14
            levs = np.linspace(vmin, vmax, levels)
            im = ax2.contourf(self.rad, self.ell[::self.minc]+1,
                              np.log10(self.e_pol_m[::self.minc,:]),
                              levs, cmap=plt.get_cmap(cm), extend='both')
            fig2.colorbar(im)
            if labTex:
                ax2.set_ylabel('Order $m+1$')
                ax2.set_xlabel('Radius $r$')
            else:
                ax2.set_ylabel('Order m+1')
                ax2.set_xlabel('Radius')
            ax2.set_yscale('log')
            ax2.set_title('E pol')

            fig3 = plt.figure()
            ax3 = fig3.add_subplot(111)
            vmax = np.log10(self.e_tor_m).max()
            vmin = vmax-14
            levs = np.linspace(vmin, vmax, levels)
            im = ax3.contourf(self.rad, self.ell[::self.minc]+1,
                              np.log10(self.e_tor_m[::self.minc,:]),
                              levs, cmap=plt.get_cmap(cm), extend='both')
            fig3.colorbar(im)
            if labTex:
                ax3.set_ylabel('Order $m+1$')
                ax3.set_xlabel('Radius $r$')
            else:
                ax3.set_ylabel('Order m+1')
                ax3.set_xlabel('Radius')
            ax3.set_yscale('log')
            ax3.set_title('E tor')
