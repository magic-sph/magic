# -*- coding: utf-8 -*-
import os, re
import matplotlib.pyplot as P
import numpy as N
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

        if field in ('eKin', 'ekin', 'e_kin', 'Ekin', 'E_kin', 'eKinR'):
            if ave:
                self.name = 'kin_spec_ave'
            else:
                self.name = 'kin_spec_'
        elif field in ('u2', 'usquare', 'u_square', 'uSquare', 'U2'):
            if ave:
                self.name = 'u2_spec_ave'
            else:
                self.name = 'u2_spec_'
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
                    print('No such tag... try again')
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
            print('No such file')
            return

        if ave is False:
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

        if iplot:
            self.plot()

    def plot(self):
        """
        Plotting function
        """
        if self.name == 'kin_spec_ave' or self.name == 'kin_spec_':
            if self.gather:
                fig = P.figure()
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
                fig = P.figure()
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

                fig = P.figure()
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
                fig = P.figure()
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
                fig = P.figure()
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

                fig = P.figure()
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
                fig = P.figure()
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
                fig = P.figure()
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

                fig = P.figure()
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
                files = scanDir('%s*%s' % (self.name, tag))
                if len(files) != 0:
                    filename = os.path.join(datadir, files[-1])
                else:
                    print('No such tag... try again')
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

        self.ell = N.arange(self.l_max+1)
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
        fig0 = P.figure()
        ax0 = fig0.add_subplot(111)
        vmax = N.log10(self.e_pol_l).max()
        vmin = vmax-14
        levs = N.linspace(vmin, vmax, levels)
        im = ax0.contourf(self.ell[1:], self.rad, N.log10(self.e_pol_l.T), 
                          levs, cmap=P.get_cmap(cm), extend='both')
        fig0.colorbar(im)
        if labTex:
            ax0.set_xlabel('Degree $\ell$')
            ax0.set_ylabel('Radius $r$')
        else:
            ax0.set_xlabel('Degree l')
            ax0.set_ylabel('Radius')
        ax0.set_xscale('log')
        ax0.set_title('E pol')

        fig1 = P.figure()
        ax1 = fig1.add_subplot(111)
        vmax = N.log10(self.e_tor_l).max()
        vmin = vmax-14
        levs = N.linspace(vmin, vmax, levels)
        im = ax1.contourf(self.ell[1:], self.rad, N.log10(self.e_tor_l.T), 
                          levs, cmap=P.get_cmap(cm), extend='both')
        fig1.colorbar(im)
        if labTex:
            ax1.set_xlabel('Degree $\ell$')
            ax1.set_ylabel('Radius $r$')
        else:
            ax1.set_xlabel('Degree l')
            ax1.set_ylabel('Radius')
        ax1.set_xscale('log')
        ax1.set_title('E tor')

        fig2 = P.figure()
        ax2 = fig2.add_subplot(111)
        vmax = N.log10(self.e_pol_m).max()
        vmin = vmax-14
        levs = N.linspace(vmin, vmax, levels)
        im = ax2.contourf(self.ell[::self.minc]+1, self.rad, 
                          N.log10(self.e_pol_m[::self.minc,:].T), 
                          levs, cmap=P.get_cmap(cm), extend='both')
        fig2.colorbar(im)
        if labTex:
            ax2.set_xlabel('Order $m+1$')
            ax2.set_ylabel('Radius $r$')
        else:
            ax2.set_xlabel('Order m+1')
            ax2.set_ylabel('Radius')
        ax2.set_xscale('log')
        ax2.set_title('E pol')

        fig3 = P.figure()
        ax3 = fig3.add_subplot(111)
        vmax = N.log10(self.e_tor_m).max()
        vmin = vmax-14
        levs = N.linspace(vmin, vmax, levels)
        im = ax3.contourf(self.ell[::self.minc]+1, self.rad, 
                          N.log10(self.e_tor_m[::self.minc,:].T), 
                          levs, cmap=P.get_cmap(cm), extend='both')
        fig3.colorbar(im)
        if labTex:
            ax3.set_xlabel('Order $m+1$')
            ax3.set_ylabel('Radius $r$')
        else:
            ax3.set_xlabel('Order m+1')
            ax3.set_ylabel('Radius')
        ax3.set_xscale('log')
        ax3.set_title('E tor')

