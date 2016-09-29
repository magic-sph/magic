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
                 ave=False, gather=False, tag=None):
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
        if self.name == 'kin_spec_ave' or self.name == 'kin_spec_':
            self.ekin_poll = data[:, 1]
            self.ekin_polm = data[:, 2]
            self.ekin_torl = data[:, 3]
            self.ekin_torm = data[:, 4]
        elif self.name == 'u2_spec_ave' or self.name == 'u2_spec_':
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
        elif self.name == 'dtVrms_spec':
            self.dtVRms = data[:, 1]
            self.CorRms = data[:, 2]
            self.LFRms = data[:, 3]
            self.AdvRms = data[:, 4]
            self.DifRms = data[:, 5]
            self.BuoRms = data[:, 6]
            self.PreRms = data[:, 7]
            self.geos = data[:, 8] # geostrophic balance
            self.mageos = data[:, 9] # magnetostrophic balance
            self.arc = data[:, 10] # archimedean balance
            self.corLor = data[:, 11] # Coriolis/Lorentz
            self.preLor = data[:, 12] # Pressure/Lorentz
            self.cia = data[:, 13] # Coriolis/Inertia/Archimedean
            self.dtVRms_SD = data[:, 14]
            self.CorRms_SD = data[:, 15]
            self.LFRms_SD = data[:, 16]
            self.AdvRms_SD = data[:, 17]
            self.DifRms_SD = data[:, 18]
            self.BuoRms_SD = data[:, 19]
            self.PreRms_SD = data[:, 20]
            self.geos_SD = data[:, 21]
            self.mageos_SD = data[:, 22]
            self.arc_SD = data[:, 23]

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
                ax.loglog(self.index, self.ekin_poll, 'k-', label='poloidal')
                ax.loglog(self.index, self.ekin_torl, 'b-', label='toroidal')
                if labTex:
                    ax.set_xlabel('Degree $\ell$')
                else:
                    ax.set_xlabel('Degree l')
                ax.set_ylabel('Kinetic energy')
                ax.set_xlim(self.index.min(), self.index.max())
                ax.legend(loc='upper right', frameon=False)
                ax = fig.add_subplot(212)
                ax.loglog(self.index[::self.minc]+1, self.ekin_polm[::self.minc],
                         'k-', label='poloidal')
                ax.loglog(self.index[::self.minc]+1, self.ekin_torm[::self.minc], 
                         'b-', label='toroidal')
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
                ax.loglog(self.index, self.ekin_poll, 'k-', label='poloidal')
                ax.loglog(self.index, self.ekin_torl, 'b-', label='toroidal')
                if labTex:
                    ax.set_xlabel('Degree $\ell$')
                else:
                    ax.set_xlabel('Degree l')
                ax.set_ylabel('Kinetic energy')
                ax.set_xlim(self.index.min(), self.index.max())
                ax.legend(loc='upper right', frameon=False)

                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.loglog(self.index[::self.minc]+1, self.ekin_polm[::self.minc],
                         'k-', label='poloidal')
                ax.loglog(self.index[::self.minc]+1, self.ekin_torm[::self.minc],
                         'b-', label='toroidal')
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
                ax.loglog(self.index, self.ekin_poll, 'k-', label='poloidal')
                ax.loglog(self.index, self.ekin_torl, 'b-', label='toroidal')
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
                         'k-', label='poloidal')
                ax.loglog(self.index[::self.minc]+1, self.ekin_torm[::self.minc], 
                         'b-', label='toroidal')
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
                ax.loglog(self.index, self.ekin_poll, 'k-', label='poloidal')
                ax.loglog(self.index, self.ekin_torl, 'b-', label='toroidal')
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
                         'k-', label='poloidal')
                ax.loglog(self.index[::self.minc]+1, self.ekin_torm[::self.minc],
                         'b-', label='toroidal')
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
                ax.loglog(self.index, self.emag_poll, 'k-', label='poloidal')
                ax.loglog(self.index, self.emag_torl, 'b-', label='toroidal')
                ax.loglog(self.index, self.emagcmb_l, 'g-', label='cmb')
                if labTex:
                    ax.set_xlabel('Degree $\ell$')
                else:
                    ax.set_xlabel('Degree l')
                ax.set_ylabel('Magnetic Energy')
                ax.set_xlim(self.index.min(), self.index.max())
                ax.legend(loc='upper right', frameon=False)

                ax = fig.add_subplot(212)
                ax.loglog(self.index[::self.minc], self.emag_polm[::self.minc],
                         'k-', label='poloidal')
                ax.loglog(self.index[::self.minc], self.emag_torm[::self.minc],
                         'b-', label='toroidal')
                ax.loglog(self.index[::self.minc], self.emagcmb_m[::self.minc],
                         'g-', label='cmb')
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
                          'k-', label='poloidal')
                ax.loglog(self.index, self.emag_torl/self.emag_torl.max(), 'b-',
                          label='toroidal')
                ax.loglog(self.index, self.emagcmb_l/self.emagcmb_l.max(), 'g-',
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
                          'k-', label='poloidal')
                ax.loglog(self.index[::self.minc]+1, 
                          self.emag_torm[::self.minc]/self.emag_torm.max(),
                          'b-', label='toroidal')
                ax.loglog(self.index[::self.minc]+1,
                          self.emagcmb_m[::self.minc]/self.emagcmb_m.max(),
                          'g-', label='cmb')
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
            ax.loglog(self.index, self.dtVRms, 'k-', label='Time derivative')
            ax.loglog(self.index, self.CorRms, 'b-', label='Coriolis')
            ax.loglog(self.index, self.PreRms, 'r-', label='Pressure')
            if self.LFRms.max() > 0.:
                ax.loglog(self.index, self.LFRms, 'c-', label='Lorentz')
            ax.loglog(self.index, self.BuoRms, 'g-', label='Buoyancy')
            ax.loglog(self.index, self.AdvRms, 'y-', label='Advection')
            ax.loglog(self.index, self.DifRms, 'm-', label='Viscosity')
            ax.loglog(self.index, self.geos, 'r--', label='Coriolis-Pressure')
            #ax.loglog(self.index, self.arc, 'b--', label='Coriolis-Pressure-Buoyancy')

            if labTex:
                ax.set_xlabel('$\ell+1$')
            else:
                ax.set_xlabel('l+1')
            ax.set_ylabel('RMS forces')
            ax.set_xlim(self.index.min(), self.index.max())
            ax.legend(loc='upper right', frameon=False)


class MagicSpectrum2D(MagicSetup):
    """
    This class can be used to read and display 2-D spectra in the :math:`(r,\ell)`
    and in the :math:`(r,m)` planes

        * Kinetic energy spectra: :ref:`2D_kin_spec_#.TAG <sec2DSpectra>`
        * Velocity square spectra: :ref:`2D_u2_spec_#.TAG <sec2DSpectra>`
        * Magnetic energy spectra: :ref:`2D_mag_spec_#.TAG <sec2DSpectra>`
    
    >>> # display the content of 2D_kin_spec_1.tag
    >>> # where tag is the most recent file in the current directory
    >>> sp = MagicSpectrum2D(field='e_kin', ispec=1, levels=17, cm='seismic')
    >>> # display the content of 2D_mag_spec_3.test
    >>> sp = MagicSpectrum2D(field='e_mag', tag='test', ispec=3)
    """

    def __init__(self, datadir='.', field='e_mag', iplot=True, ispec=None, 
                 tag=None, cm='jet', levels=33, precision='Float64'):
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
        """

        if field in ('eKin', 'ekin', 'e_kin', 'Ekin', 'E_kin', 'eKinR'):
            self.name = '2D_kin_spec_'
        elif field in ('u2'):
            self.name = '2D_u2_spec_'
        elif field in('eMag', 'emag', 'e_mag', 'Emag', 'E_mag', 'eMagR'):
            self.name = '2D_mag_spec_'

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

        file = npfile(filename, endian='B')

        out = file.fort_read('%s,3i4' % precision)[0]
        self.time = out[0]
        self.n_r_max, self.l_max, self.minc = out[1]
        self.rad = file.fort_read(precision, shape=(self.n_r_max))
        self.e_pol_l = file.fort_read(precision, shape=(self.l_max, self.n_r_max))
        self.e_pol_m = file.fort_read(precision, shape=(self.l_max+1, self.n_r_max))
        self.e_tor_l = file.fort_read(precision, shape=(self.l_max, self.n_r_max))
        self.e_tor_m = file.fort_read(precision, shape=(self.l_max+1, self.n_r_max))

        self.ell = np.arange(self.l_max+1)
        file.close()

        if iplot:
            self.plot(levels, cm)


    def plot(self, levels, cm):
        """
        Plotting function

        :param levels: number of contour levels
        :type levels: int
        :param cm: name of the colormap
        :type cm: str
        """
        fig0 = plt.figure()
        ax0 = fig0.add_subplot(111)
        vmax = np.log10(self.e_pol_l).max()
        vmin = vmax-14
        levs = np.linspace(vmin, vmax, levels)
        im = ax0.contourf(self.ell[1:], self.rad, np.log10(self.e_pol_l.T), 
                          levs, cmap=plt.get_cmap(cm), extend='both')
        fig0.colorbar(im)
        if labTex:
            ax0.set_xlabel('Degree $\ell$')
            ax0.set_ylabel('Radius $r$')
        else:
            ax0.set_xlabel('Degree l')
            ax0.set_ylabel('Radius')
        ax0.set_xscale('log')
        ax0.set_title('E pol')

        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        vmax = np.log10(self.e_tor_l).max()
        vmin = vmax-14
        levs = np.linspace(vmin, vmax, levels)
        im = ax1.contourf(self.ell[1:], self.rad, np.log10(self.e_tor_l.T), 
                          levs, cmap=plt.get_cmap(cm), extend='both')
        fig1.colorbar(im)
        if labTex:
            ax1.set_xlabel('Degree $\ell$')
            ax1.set_ylabel('Radius $r$')
        else:
            ax1.set_xlabel('Degree l')
            ax1.set_ylabel('Radius')
        ax1.set_xscale('log')
        ax1.set_title('E tor')

        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        vmax = np.log10(self.e_pol_m).max()
        vmin = vmax-14
        levs = np.linspace(vmin, vmax, levels)
        im = ax2.contourf(self.ell[::self.minc]+1, self.rad, 
                          np.log10(self.e_pol_m[::self.minc,:].T), 
                          levs, cmap=plt.get_cmap(cm), extend='both')
        fig2.colorbar(im)
        if labTex:
            ax2.set_xlabel('Order $m+1$')
            ax2.set_ylabel('Radius $r$')
        else:
            ax2.set_xlabel('Order m+1')
            ax2.set_ylabel('Radius')
        ax2.set_xscale('log')
        ax2.set_title('E pol')

        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        vmax = np.log10(self.e_tor_m).max()
        vmin = vmax-14
        levs = np.linspace(vmin, vmax, levels)
        im = ax3.contourf(self.ell[::self.minc]+1, self.rad, 
                          np.log10(self.e_tor_m[::self.minc,:].T), 
                          levs, cmap=plt.get_cmap(cm), extend='both')
        fig3.colorbar(im)
        if labTex:
            ax3.set_xlabel('Order $m+1$')
            ax3.set_ylabel('Radius $r$')
        else:
            ax3.set_xlabel('Order m+1')
            ax3.set_ylabel('Radius')
        ax3.set_xscale('log')
        ax3.set_title('E tor')

