# -*- coding: utf-8 -*-
import os
import re
import copy
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
    >>> sp = MagicSpectrum(field='e_mag', tag='test', ave=True)
    """

    def __init__(self, datadir='.', field='e_kin', iplot=True, ispec=None,
                 ave=False, normalize=False, tag=None, quiet=False):
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
        :param datadir: current working directory
        :type datadir: str
        :param quiet: when set to True, makes the output silent (default False)
        :type quiet: bool
        """
        self.normalize = normalize
        self.ave = ave
        if field in ('eKin', 'ekin', 'e_kin', 'Ekin', 'E_kin', 'eKinR', 'kin'):
            if self.ave:
                self.name = 'kin_spec_ave'
            else:
                self.name = 'kin_spec_'
        elif field in ('u2', 'usquare', 'u_square', 'uSquare', 'U2'):
            if self.ave:
                self.name = 'u2_spec_ave'
            else:
                self.name = 'u2_spec_'
        elif field in ('eMag', 'emag', 'e_mag', 'Emag', 'E_mag', 'eMagR', 'mag'):
            if self.ave:
                self.name = 'mag_spec_ave'
            else:
                self.name = 'mag_spec_'
        elif field in ('dtVrms'):
            self.name = 'dtVrms_spec'
            self.ave = True

        elif field in ('T','temperature','S','entropy'):
            self.name = 'T_spec_'

        if self.ave: # Time-averaged spectra

            if tag is not None:
                pattern = os.path.join(datadir, '{}.{}'.format(self.name, tag))
                files = scanDir(pattern)
                # Either the log.tag directly exists and the setup is easy to
                # obtain
                if os.path.exists(os.path.join(datadir, 'log.{}'.format(tag))):
                    MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                        nml='log.{}'.format(tag))
                # Or the tag is a bit more complicated and we need to find
                # the corresponding log file
                else:
                    mask = re.compile(r'{}\/{}\.(.*)'.format(datadir, self.name))
                    if mask.match(files[-1]):
                        ending = mask.search(files[-1]).groups(0)[0]
                        pattern = os.path.join(datadir, 'log.{}'.format(ending))
                        if os.path.exists(pattern):
                            MagicSetup.__init__(self, datadir=datadir,
                                                quiet=True, nml='log.{}'.format(ending))

                # Sum the files that correspond to the tag
                mask = re.compile(r'{}\.(.*)'.format(self.name))
                for k, file in enumerate(files):
                    if not quiet:
                        print('reading {}'.format(file))

                    tag = mask.search(file).groups(0)[0]
                    nml = MagicSetup(nml='log.{}'.format(tag), datadir=datadir,
                                     quiet=True)
                    filename = file
                    data = fast_read(filename)

                    if k == 0:
                        speclut = SpecLookUpTable(data, self.name, nml.start_time,
                                                  nml.stop_time)
                    else:
                        speclut += SpecLookUpTable(data, self.name, nml.start_time,
                                                    nml.stop_time)

            else: # Tag is None: take the most recent one
                pattern = os.path.join(datadir, '{}.*'.format(self.name))
                files = scanDir(pattern)
                filename = files[-1]
                if not quiet:
                    print('reading {}'.format(filename))
                # Determine the setup
                mask = re.compile(r'{}\.(.*)'.format(self.name))
                ending = mask.search(files[-1]).groups(0)[0]
                if os.path.exists('log.{}'.format(ending)):
                    try:
                        MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml='log.{}'.format(ending))
                    except AttributeError:
                        self.start_time = None
                        self.stop_time = None
                        pass

                if not hasattr(self, 'stop_time'):
                    self.stop_time = None
                data = fast_read(filename)
                speclut = SpecLookUpTable(data, self.name, self.start_time,
                                          self.stop_time)

        else: # Snapshot spectra

            if tag is not None:
                if ispec is not None:
                    file = '{}{}.{}'.format(self.name, ispec, tag)
                    filename = os.path.join(datadir, file)
                else:
                    pattern = os.path.join(datadir, '{}*{}'.format(self.name, tag))
                    files = scanDir(pattern)
                    if len(files) != 0:
                        filename = files[-1]
                    else:
                        print('No such tag... try again')
                        return

                if os.path.exists(os.path.join(datadir, 'log.{}'.format(tag))):
                    try:
                        MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml='log.{}'.format(tag))
                    except AttributeError:
                        self.start_time = None
                        self.stop_time = None
            else:
                if ispec is not None:
                    pattern = os.path.join(datadir, '{}{}*'.format(self.name, ispec))
                    files = scanDir(pattern)
                    filename = files[-1]
                else:
                    pattern = os.path.join(datadir, '{}*'.format(self.name))
                    files = scanDir(pattern)
                    filename = files[-1]

                # Determine the setup
                mask = re.compile(r'.*\.(.*)')
                ending = mask.search(files[-1]).groups(0)[0]
                if os.path.exists(os.path.join(datadir, 'log.{}'.format(ending))):
                    try:
                        MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml='log.{}'.format(ending))
                    except AttributeError:
                        self.start_time = None
                        self.stop_time = None
                        pass

            if not quiet:
                print('reading {}'.format(filename))
            if not hasattr(self, 'stop_time'):
                self.stop_time = None
            data = fast_read(filename, skiplines=1)
            speclut = SpecLookUpTable(data, self.name, self.start_time,
                                      self.stop_time)

        # Copy look-up table arguments into MagicSpectrum object
        for attr in speclut.__dict__:
            setattr(self, attr, speclut.__dict__[attr])

        if iplot:
            self.plot()

    def plot(self):
        """
        Plotting function
        """
        if self.name == 'kin_spec_ave' or self.name == 'kin_spec_':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            if self.normalize:
                y = self.ekin_poll+self.ekin_torl
                ax.loglog(self.index[1:], y[1:]/y[1:].max(),)
            else:
                ax.loglog(self.index[1:], self.ekin_poll[1:], label='poloidal')
                if self.ave:
                    ax.fill_between(self.index[1:], self.ekin_poll[1:]-\
                                    self.ekin_poll_SD[1:], self.ekin_poll[1:]+\
                                    self.ekin_poll_SD[1:], alpha=0.2)
                ax.loglog(self.index[1:], self.ekin_torl[1:], label='toroidal')
                if self.ave:
                    ax.fill_between(self.index[1:], self.ekin_torl[1:]-\
                                    self.ekin_torl_SD[1:], self.ekin_torl[1:]+\
                                    self.ekin_torl_SD[1:], alpha=0.2)
            if labTex:
                ax.set_xlabel('Degree $\ell$')
            else:
                ax.set_xlabel('Degree l')
            ax.set_ylabel('Kinetic energy')
            ax.set_xlim(1, self.index[-1])
            ax.legend(loc='upper right', frameon=False)
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            if self.normalize:
                y = self.ekin_polm[::self.minc]+self.ekin_torm[::self.minc]
                ax.loglog(self.index[::self.minc]+1, y/y.max())
            else:
                ax.loglog(self.index[::self.minc]+1, self.ekin_polm[::self.minc],
                          label='poloidal')
                if self.ave:
                    ax.fill_between(self.index[::self.minc]+1,
                                    self.ekin_polm[::self.minc]-\
                                    self.ekin_polm_SD[::self.minc],
                                    self.ekin_polm[::self.minc]+\
                                    self.ekin_polm_SD[::self.minc], alpha=0.2)
                ax.loglog(self.index[::self.minc]+1, self.ekin_torm[::self.minc],
                          label='toroidal')
                if self.ave:
                    ax.fill_between(self.index[::self.minc]+1,
                                    self.ekin_torm[::self.minc]-\
                                    self.ekin_torm_SD[::self.minc],
                                    self.ekin_torm[::self.minc]+\
                                    self.ekin_torm_SD[::self.minc], alpha=0.2)
            if labTex:
                ax.set_xlabel('$m$ + 1')
            else:
                ax.set_xlabel('m + 1')
            ax.set_ylabel('Kinetic energy')
            ax.set_xlim(1, self.index[-1]+1)
            ax.legend(loc='upper right', frameon=False)
            fig.tight_layout()
        elif self.name == 'u2_spec_ave' or self.name == 'u2_spec_':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.loglog(self.index[1:], self.ekin_poll[1:], label='poloidal')
            if self.ave:
                ax.fill_between(self.index[1:], self.ekin_poll[1:]-\
                                self.ekin_poll_SD[1:], self.ekin_poll[1:]+\
                                self.ekin_poll_SD[1:], alpha=0.2)
            ax.loglog(self.index[1:], self.ekin_torl[1:], label='toroidal')
            if self.ave:
                ax.fill_between(self.index[1:], self.ekin_torl[1:]-\
                                self.ekin_torl_SD[1:], self.ekin_torl[1:]+\
                                self.ekin_torl_SD[1:], alpha=0.2)
            if labTex:
                ax.set_xlabel('Degree $\ell$')
                ax.set_ylabel(r'${\cal U}^2$')
            else:
                ax.set_xlabel('Degree l')
                ax.set_ylabel('Velocity square')
            ax.set_xlim(1, self.index[-1])
            ax.legend(loc='upper right', frameon=False)
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.loglog(self.index[::self.minc]+1, self.ekin_polm[::self.minc],
                      label='poloidal')
            if self.ave:
                ax.fill_between(self.index[::self.minc]+1,
                                self.ekin_polm[::self.minc]-\
                                self.ekin_polm_SD[::self.minc],
                                self.ekin_polm[::self.minc]+\
                                self.ekin_polm_SD[::self.minc], alpha=0.2)
            ax.loglog(self.index[::self.minc]+1, self.ekin_torm[::self.minc],
                      label='toroidal')
            if self.ave:
                ax.fill_between(self.index[::self.minc]+1,
                                self.ekin_torm[::self.minc]-\
                                self.ekin_torm_SD[::self.minc],
                                self.ekin_torm[::self.minc]+\
                                self.ekin_torm_SD[::self.minc], alpha=0.2)
            if labTex:
                ax.set_xlabel(r'$m + 1$')
                ax.set_ylabel(r'${\cal U}^2$')
            else:
                ax.set_xlabel('m + 1')
                ax.set_ylabel('Velocity square')
            ax.set_xlim(1, self.index.max[-1]+1)
            ax.legend(loc='upper right', frameon=False)
            fig.tight_layout()
        elif self.name == 'mag_spec_ave' or self.name == 'mag_spec_':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.loglog(self.index[1:], self.emag_poll[1:], label='poloidal')
            if self.ave:
                ax.fill_between(self.index[1:], self.emag_poll[1:]-\
                                self.emag_poll_SD[1:], self.emag_poll[1:]+\
                                self.emag_poll_SD[1:], alpha=0.2)
            ax.loglog(self.index[1:], self.emag_torl[1:], label='toroidal')
            if self.ave:
                ax.fill_between(self.index[1:], self.emag_torl[1:]-\
                                self.emag_torl_SD[1:], self.emag_torl[1:]+\
                                self.emag_torl_SD[1:], alpha=0.2)
            ax.loglog(self.index[1:], self.emagcmb_l[1:], label='cmb')
            if self.ave:
                ax.fill_between(self.index[1:], self.emagcmb_l[1:]-\
                                self.emagcmb_l_SD[1:], self.emagcmb_l[1:]+\
                                self.emagcmb_l_SD[1:], alpha=0.2)
            if labTex:
                ax.set_xlabel('Degree $\ell$')
            else:
                ax.set_xlabel('Degree l')
            ax.set_ylabel('Magnetic Energy')
            ax.set_xlim(1, self.index[-1])
            ax.legend(loc='upper right', frameon=False)
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.loglog(self.index[::self.minc]+1, self.emag_polm[::self.minc],
                      label='poloidal')
            if self.ave:
                ax.fill_between(self.index[::self.minc]+1,
                                self.emag_polm[::self.minc]-\
                                self.emag_polm_SD[::self.minc],
                                self.emag_polm[::self.minc]+\
                                self.emag_polm_SD[::self.minc], alpha=0.2)
            ax.loglog(self.index[::self.minc]+1, self.emag_torm[::self.minc],
                      label='toroidal')
            if self.ave:
                ax.fill_between(self.index[::self.minc]+1,
                                self.emag_torm[::self.minc]-\
                                self.emag_torm_SD[::self.minc],
                                self.emag_torm[::self.minc]+\
                                self.emag_torm_SD[::self.minc], alpha=0.2)
            ax.loglog(self.index[::self.minc]+1, self.emagcmb_m[::self.minc],
                      label='cmb')
            if self.ave:
                ax.fill_between(self.index[::self.minc]+1,
                                self.emagcmb_m[::self.minc]-\
                                self.emagcmb_m_SD[::self.minc],
                                self.emagcmb_m[::self.minc]+\
                                self.emagcmb_m_SD[::self.minc], alpha=0.2)
            if labTex:
                ax.set_xlabel('$m$+1')
            else:
                ax.set_xlabel('m+1')
            ax.set_ylabel('Magnetic energy')
            ax.set_xlim(1, self.index[-1]+1)
            ax.legend(loc='upper right', frameon=False)
            fig.tight_layout()

        elif self.name == 'dtVrms_spec':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.loglog(self.index[1:], self.CorRms[1:], label='Coriolis')
            ax.fill_between(self.index[1:], self.CorRms[1:]-self.CorRms_SD[1:,],
                            self.CorRms[1:]+self.CorRms_SD[1:,], alpha=0.2)
            ax.loglog(self.index[1:], self.PreRms[1:], label='Pressure')
            ax.fill_between(self.index[1:], self.PreRms[1:]-self.PreRms_SD[1:,],
                            self.PreRms[1:]+self.PreRms_SD[1:,], alpha=0.2)
            if self.LFRms.max() > 1e-10:
                ax.loglog(self.index[1:], self.LFRms[1:], label='Lorentz')
                ax.fill_between(self.index[1:], self.LFRms[1:]-self.LFRms_SD[1:,],
                                self.LFRms[1:]+self.LFRms_SD[1:,], alpha=0.2)
            if self.BuoRms.max() > 1e-10:
                ax.loglog(self.index[1:], self.BuoRms[1:], label='Buoyancy')
                ax.fill_between(self.index[1:], self.BuoRms[1:]-self.BuoRms_SD[1:,],
                                self.BuoRms[1:]+self.BuoRms_SD[1:,], alpha=0.2)
            if self.ChemRms.max() > 1e-10:
                ax.loglog(self.index[1:], self.ChemRms[1:], label='Chem. Buoyancy')
                ax.fill_between(self.index[1:], self.ChemRms[1:]-self.ChemRms_SD[1:,],
                                self.ChemRms[1:]+self.ChemRms_SD[1:,], alpha=0.2)
            ax.loglog(self.index[1:], self.InerRms[1:], label='Inertia')
            ax.fill_between(self.index[1:], self.InerRms[1:]-self.InerRms_SD[1:,],
                            self.InerRms[1:]+self.InerRms_SD[1:,], alpha=0.2)
            ax.loglog(self.index[1:], self.DifRms[1:], label='Viscosity')
            ax.fill_between(self.index[1:], self.DifRms[1:]-self.DifRms_SD[1:,],
                            self.DifRms[1:]+self.DifRms_SD[1:,], alpha=0.2)
            ax.loglog(self.index[1:], self.geos[1:], label='Coriolis-Pressure')
            ax.fill_between(self.index[1:], self.geos[1:]-self.geos_SD[1:,],
                            self.geos[1:]+self.geos_SD[1:,], alpha=0.2)
            #ax.loglog(self.index[1:], self.arcMag[1:], ls='--',
            #          label='Coriolis-Pressure-Buoyancy-Lorentz')

            if labTex:
                ax.set_xlabel('$\ell$')
            else:
                ax.set_xlabel('l')
            ax.set_ylabel('RMS forces')
            ax.set_xlim(1, self.index[-1])
            ax.legend(loc='lower right', frameon=False, ncol=2)
            fig.tight_layout()

        elif self.name == 'T_spec_':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.loglog(self.index[1:], self.T_l[1:])
            ax.loglog(self.index[1:], self.T_icb_l[1:])
            if labTex:
                ax.set_xlabel('$\ell$')
            else:
                ax.set_xlabel('l')
            ax.set_xlim(1, self.index[-1])
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.loglog(self.index[::self.minc]+1,
                      self.T_m[::self.minc],color='g')
            if labTex:
                ax.set_xlabel('Order $m+1$')
            else:
                ax.set_xlabel('m+1')
            ax.set_xlim(1, self.index[-1]+1)
            fig.tight_layout()

class SpecLookUpTable:
    """
    The purpose of this class is to create a lookup table between the numpy
    array that comes from the reading of the spec files and the corresponding
    columns.
    """

    def __init__(self, data, name, tstart=None, tstop=None):
        """
        :param data: numpy array that contains the data
        :type data: numpy.ndarray
        :param name: name of the field (i.e. 'eKinR', 'eMagR', 'powerR', ...)
        :type name: str
        :param tstart: starting time that was used to compute the time average
        :type tstart: float
        :param tstop: stop time that was used to compute the time average
        :type tstop: float
        """

        self.name = name
        self.start_time = tstart
        self.stop_time = tstop

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

            self.index = self.index-1
        elif self.name == 'T_spec_':
            self.T_l = data[:,1]
            self.T_m = data[:,2]
            self.T_icb_l = data[:,3]
            self.T_icb_m = data[:,4]
            self.dT_icb_l = data[:,5]
            self.dT_icb_m = data[:,6]

    def __add__(self, new):
        """
        This is a python built-in method to stack two look-up tables.
        """

        out = copy.deepcopy(new)
        if self.start_time is not None:
            fac_old = self.stop_time-self.start_time
            out.start_time = self.start_time
        else:
            fac_old = 0.
        if new.stop_time is not None:
            fac_new = new.stop_time-new.start_time
            out.stop_time = new.stop_time
        else:
            fac_new = 0.
        if fac_old != 0 or fac_new != 0:
            fac_tot = fac_new+fac_old
        else:
            fac_tot = 1.

        idx_old_max = len(self.index)
        idx_new_max = len(new.index)

        if idx_old_max == idx_new_max:
            for attr in new.__dict__.keys():
                if attr not in ['index', 'name', 'start_time', 'stop_time']:
                    # Only stack if both new and old have the attribute available
                    if attr in self.__dict__:
                        # Standard deviation
                        if attr.endswith('SD'):
                            if abs(self.__dict__[attr]).max() > 0.:
                                out.__dict__[attr] = np.sqrt(( \
                                                  fac_new*new.__dict__[attr]**2 + \
                                                  fac_old*self.__dict__[attr]**2) / \
                                                  fac_tot)
                            else:
                                out.__dict__[attr] = new.__dict__[attr]
                        # Regular field
                        else:
                            if abs(self.__dict__[attr]).max() > 0.:
                                out.__dict__[attr] = (fac_new*new.__dict__[attr] + \
                                                      fac_old*self.__dict__[attr]) / \
                                                      fac_tot
                            else:
                                out.__dict__[attr] = new.__dict__[attr]
        else: # Different truncations
            idx_min = min(idx_old_max, idx_new_max)
            for attr in new.__dict__.keys():
                if attr not in ['index', 'name', 'start_time', 'stop_time']:
                    # Only stack if both new and old have the attribute available
                    if attr in self.__dict__:
                        # Standard deviation
                        if attr.endswith('SD'):
                            if abs(self.__dict__[attr]).max() > 0.:
                                out.__dict__[attr][:idx_min] = \
                                   np.sqrt((fac_new*new.__dict__[attr][:idx_min]**2+\
                                            fac_old*self.__dict__[attr][:idx_min]**2)\
                                            /fac_tot)
                            else:
                                out.__dict__[attr] = new.__dict__[attr]
                        # Regular field
                        else:
                            if abs(self.__dict__[attr]).max() > 0.:
                                out.__dict__[attr][:idx_min] = \
                                   (fac_new*new.__dict__[attr][:idx_min] + \
                                    fac_old*self.__dict__[attr][:idx_min]) / fac_tot
                            else:
                                out.__dict__[attr] = new.__dict__[attr]

        return out


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
                file = '{}{}.{}'.format(self.name, ispec, tag)
                filename = os.path.join(datadir, file)
            else:
                pattern = os.path.join(datadir, '{}*{}'.format(self.name, tag))
                files = scanDir(pattern)
                if len(files) != 0:
                    filename = files[-1]
                else:
                    print('No such tag... try again')
                    return

            if os.path.exists(os.path.join(datadir, 'log.{}'.format(tag))):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.{}'.format(tag))
        else:
            if ispec is not None:
                pattern = os.path.join(datadir, '{}{}*'.format(self.name, ispec))
                files = scanDir(pattern)
                filename = files[-1]
            else:
                pattern = os.path.join(datadir, '{}*'.format(self.name))
                files = scanDir(pattern)
                filename = files[-1]
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists(os.path.join(datadir, 'log.{}'.format(ending))):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.{}'.format(ending))

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
