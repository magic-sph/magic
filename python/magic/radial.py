# -*- coding: utf-8 -*-
import os
import re
import copy
import matplotlib.pyplot as plt
import numpy as np
from .log import MagicSetup
from .libmagic import fast_read,scanDir
import scipy.interpolate as sint


class MagicRadial(MagicSetup):
    """
    This class can be used to read and display the time and
    horizontally averaged files:

        * Kinetic energy: :ref:`eKinR.TAG <secEkinRFile>`
        * Magnetic energy: :ref:`eMagR.TAG <secEmagRfile>`
        * Anelastic reference state: :ref:`anel.TAG <secAnelFile>`
        * Variable electrical conductivity: :ref:`varCond.TAG <secVarCondFile>`
        * Variable thermal diffusivity: :ref:`varDiff.TAG <secVarDiffFile>`
        * Variable kinematic viscosity: :ref:`varVisc.TAG <secVarViscFile>`
        * Diagnostic parameters: :ref:`parR.TAG <secPaRfile>`
        * Power budget: :ref:`powerR.TAG <secPowerRfile>`
        * Heat fluxes: :ref:`fluxesR.TAG <secFluxesRfile>`
        * Mean entropy, temperature and pressure: :ref:`heatR.TAG <secHeatRfile>`
        * Radial profiles used for boundary layers: :ref:`bLayersR.TAG <secBLayersRfile>`
        * Parallel/perpendicular decomposition: :ref:`perpParR.TAG <secPerpParRfile>`

    >>> rad = MagicRadial(field='eKinR') # display the content of eKinR.tag
    >>> print(rad.radius, rad.ekin_pol_axi) # print radius and poloidal energy
    """

    def __init__(self, datadir='.', field='eKin', iplot=True, tag=None, tags=None,
                 normalize_radius=False, quiet=False):
        """
        :param datadir: working directory
        :type datadir: str
        :param field: the field you want to plot
        :type field: str
        :param iplot: to plot the output, default is True
        :type iplot: bool
        :param tag: a specific tag, default is None
        :type tag: str
        :param tags: a list that contains multiple tags: useful to sum
                     several radial files
        :type tags: list
        :param quiet: when set to True, makes the output silent (default False)
        :type quiet: bool
        """

        if field in ('eKin', 'ekin', 'e_kin', 'Ekin', 'E_kin', 'eKinR'):
            self.name ='eKinR'
        elif field in ('eMag', 'emag', 'e_mag', 'Emag', 'E_mag', 'eMagR'):
            self.name = 'eMagR'
        elif field in ('anel', 'ref', 'rho', 'density', 'background'):
            self.name = 'anel'
        elif field in ('varDiff', 'vardiff'):
            self.name = 'varDiff'
        elif field in ('varVisc', 'varvisc'):
            self.name = 'varVisc'
        elif field in ('varCond', 'varcond'):
            self.name = 'varCond'
        elif field in ('power', 'powerR'):
            self.name = 'powerR'
        elif field in ('parrad'):
            self.name = 'parrad'
        elif field in ('bLayersR'):
            self.name = 'bLayersR'
        elif field in ('parR'):
            self.name = 'parR'
        elif field in ('fluxesR'):
            self.name = 'fluxesR'
        elif field in ('heatR'):
            self.name = 'heatR'
        elif field in ('perpParR'):
            self.name = 'perpParR'
        else:
            print('No corresponding radial profiles... Try again')

        self.normalize_radius = normalize_radius

        if tags is None:
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
                    if self.name == 'varCond' or self.name == 'varVisc' or \
                       self.name == 'varDiff' or self.name == 'anel':
                        data = fast_read(filename, skiplines=1)
                    else:
                        data = fast_read(filename, skiplines=0)

                    if k == 0:
                        radlut = RadLookUpTable(data, self.name, nml.start_time,
                                                nml.stop_time)
                    else:
                        radlut += RadLookUpTable(data, self.name, nml.start_time,
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

                if self.name == 'varCond' or self.name == 'varVisc' or \
                   self.name == 'varDiff' or self.name == 'anel':
                    data = fast_read(filename, skiplines=1)
                else:
                    data = fast_read(filename, skiplines=0)

                radlut = RadLookUpTable(data, self.name, self.start_time,
                                        self.stop_time)
        else:
            if os.path.exists(os.path.join(datadir, 'log.{}'.format(tags[-1]))):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.{}'.format(tags[-1]))
            else:
                self.start_time = None
                self.stop_time = None
            for k, tagg in enumerate(tags):
                nml = MagicSetup(nml='log.{}'.format(tagg), datadir=datadir,
                                 quiet=True)
                file = '{}.{}'.format(self.name, tagg)
                filename = os.path.join(datadir, file)
                if not quiet:
                    print('reading {}'.format(filename))
                if self.name == 'varCond' or self.name == 'varVisc' or \
                   self.name == 'varDiff' or self.name == 'anel':
                    data = fast_read(filename, skiplines=1)
                else:
                    data = fast_read(filename, skiplines=0)

                if k == 0:
                    radlut = RadLookUpTable(data, self.name, nml.start_time,
                                            nml.stop_time)
                else:
                    radlut += RadLookUpTable(data, self.name, nml.start_time,
                                             nml.stop_time)

        # Copy look-up table arguments into MagicRadial object
        for attr in radlut.__dict__:
            setattr(self, attr, radlut.__dict__[attr])

        # Plotting function
        if iplot:
            self.plot()

    def plot(self):
        """
        Display the result when ``iplot=True``
        """
        if self.normalize_radius:
            x_axis = self.radius/self.radius[0]
        else:
            x_axis = self.radius

        if self.name == 'eKinR':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_axis, self.ekin_pol, ls='-', c='#30a2da',
                    label='ekin pol')
            ax.plot(x_axis, self.ekin_tor, ls='-', c='#fc4f30',
                    label='ekin tor')
            ax.plot(x_axis, self.ekin_pol_axi, ls='--', c='#30a2da',
                    label='ekin pol axi')
            ax.plot(x_axis, self.ekin_tor_axi, ls='--', c='#fc4f30',
                    label='ekin tor axi')
            ax.plot(x_axis, self.ekin_pol+self.ekin_tor, ls='-', c='#31363B',
                    label='ekin tot')
            ax.set_xlabel('Radius')
            ax.set_ylabel('Kinetic energy')
            ax.set_xlim(x_axis.min(), x_axis.max())
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()
        elif self.name == 'eMagR':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_axis, self.emag_pol, ls='-', c='#30a2da',
                    label='emag pol')
            ax.plot(x_axis, self.emag_tor, ls='-', c='#fc4f30',
                    label='emag tor')
            ax.plot(x_axis, self.emag_pol_axi, ls='--', c='#30a2da',
                    label='emag pol axi')
            ax.plot(x_axis, self.emag_tor_axi, ls='--', c='#fc4f30',
                    label='emag tor axi')
            ax.plot(x_axis, self.emag_pol+self.emag_tor, ls='-', c='#31363B',
                    label='emag tot')

            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Radius')
            ax.set_ylabel('Magnetic energy')
            ax.set_xlim(x_axis.min(), x_axis.max())
            if hasattr(self, 'con_radratio'):
                if self.nVarCond == 2:
                    ax.axvline(self.con_radratio*self.radius[0], color='k',
                              linestyle='--')
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_axis, self.dip_ratio)
            ax.set_xlabel('Radius')
            ax.set_ylabel('Dipolarity')
            ax.set_xlim(x_axis.min(), x_axis.max())
            fig.tight_layout()
        elif self.name == 'anel':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(x_axis, self.temp0, label='Temperature')
            ax.semilogy(x_axis, self.rho0, label='Density')
            ax.set_ylabel('Reference state')
            ax.set_xlabel('Radius')
            ax.set_xlim(x_axis.min(), x_axis.max())
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_axis, self.beta, label='beta')
            ax.plot(x_axis, self.dbeta, label='dbeta/dr')
            ax.set_xlabel('Radius')
            ax.set_ylabel('Derivatives of rho')
            ax.set_xlim(x_axis.min(), x_axis.max())
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_axis, self.grav)
            ax.set_xlabel('Radius')
            ax.set_ylabel('Gravity')
            ax.set_xlim(x_axis.min(), x_axis.max())
            fig.tight_layout()
        elif self.name == 'varDiff':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(x_axis, self.conduc, label='conductivity')
            ax.semilogy(x_axis, self.kappa, label='diffusivity')
            ax.semilogy(x_axis, self.prandtl, label='Prandtl')
            ax.set_ylabel('Thermal properties')
            ax.set_xlabel('Radius')
            ax.set_xlim(x_axis.min(), x_axis.max())
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_axis, self.dLkappa, label='dLkappa')
            ax.set_xlabel('Radius')
            ax.set_ylabel('$d\ln\kappa / dr$')
            ax.set_xlim(x_axis.min(), x_axis.max())
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()
        elif self.name == 'varVisc':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(x_axis, self.dynVisc, label='dyn. visc')
            ax.semilogy(x_axis, self.kinVisc, label='kin. visc')
            ax.semilogy(x_axis, self.prandtl, label='Prandtl')
            if self.mode not in (1,5):
                ax.semilogy(x_axis, self.prandtlmag, label='Pm')
            ax.semilogy(x_axis, self.ekman, label='Ekman')
            ax.set_ylabel('Thermal properties')
            ax.set_xlabel('Radius')
            ax.set_xlim(x_axis.min(), x_axis.max())
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_axis, self.dLvisc, label='dLvisc')
            ax.set_xlabel('Radius')
            ax.set_ylabel(r'$d\ln\nu / dr$')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()
        elif self.name == 'varCond':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(x_axis, self.conduc, label='conductivity')
            ax.set_ylabel('Electrical conductivity')
            ax.set_xlabel('Radius')
            ax.set_xlim(x_axis.min(), x_axis.max())
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()
        elif self.name == 'powerR':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.buoPower, label='Power therm')
            if hasattr(self, 'buoPower_SD'):
                ax.fill_between(x_axis, self.buoPower-self.buoPower_SD,
                                self.buoPower+self.buoPower_SD,
                                color='C0', alpha=0.3)
            ax.plot(self.radius, self.viscDiss, label='visc diss')
            if hasattr(self, 'viscDiss_SD'):
                ax.fill_between(x_axis, self.viscDiss-self.viscDiss_SD,
                                self.viscDiss+self.viscDiss_SD,
                                color='C1', alpha=0.3)
            if abs(self.buoPower_chem).max() > 0:
                ax.plot(self.radius, self.buoPower_chem, color='C3',
                        label='Power chem')
                if hasattr(self, 'buoPower_chem_SD'):
                    ax.fill_between(x_axis, self.buoPower_chem-self.buoPower_chem_SD,
                                    self.buoPower_chem+self.buoPower_chem_SD,
                                    color='C3', alpha=0.3)
            if self.ohmDiss.max() != 0.:
                ax.plot(x_axis, self.ohmDiss, label='ohm diss', color='C2')
                if hasattr(self, 'ohmDiss_SD'):
                    ax.fill_between(x_axis, self.ohmDiss-self.ohmDiss_SD,
                                    self.ohmDiss+self.ohmDiss_SD,
                                    color='C2', alpha=0.3)
            ax.set_xlabel('Radius')
            ax.set_xlim(x_axis.min(), x_axis.max())
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()
        elif self.name == 'parrad' or self.name == 'bLayersR':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_axis, self.entropy, label='Entropy')
            ax.twinx()
            ax.plot(x_axis, self.entropy_SD/self.entropy_SD.max(),
                    label='Standard dev. of Entropy')
            ax.set_xlabel('Radius')
            ax.set_xlim(x_axis.min(), x_axis.max())
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_axis, self.uh, label='uh')
            if hasattr(self, 'uh_SD'):
                ax.fill_between(x_axis, self.uh-self.uh_SD, self.uh+self.uh_SD,
                                color='C0', alpha=0.3)
            ax.plot(x_axis, self.duhdr, label='duhdr')
            if hasattr(self, 'duhdr_SD'):
                ax.fill_between(x_axis, self.duhdr-self.duhdr_SD,
                                self.duhdr+self.duhdr_SD,
                                color='C1', alpha=0.3)
            ax.set_xlabel('Radius')
            ax.set_xlim(x_axis.min(), x_axis.max())
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()
        elif self.name == 'parR':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_axis, self.rm)
            if hasattr(self, 'rm_SD'):
                ax.fill_between(x_axis, self.rm-self.rm_SD, self.rm+self.rm_SD,
                                color='C0', alpha=0.3)
            ax.set_xlim(x_axis.min(), x_axis.max())
            ax.set_xlabel('Radius')
            ax.set_ylabel('Rm')
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_axis, self.rol, label='Rol')
            if hasattr(self, 'rol_SD'):
                ax.fill_between(x_axis, self.rol-self.rol_SD, self.rol+self.rol_SD,
                                color='C1', alpha=0.3)
            ax.plot(x_axis, self.urol, label='u Rol')
            if hasattr(self, 'urol_SD'):
                ax.fill_between(x_axis, self.urol-self.urol_SD,
                                self.urol+self.urol_SD,
                                color='C2', alpha=0.3)
            ax.set_xlabel('Radius')
            ax.set_ylabel('Rol')
            ax.set_xlim(x_axis.min(), x_axis.max())
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_axis[1:-1], self.dlV[1:-1], label='dlV')
            if hasattr(self, 'dlV_SD'):
                ax.fill_between(x_axis[1:-1], self.dlV[1:-1]-self.dlV_SD[1:-1],
                                self.dlV[1:-1]+self.dlV_SD[1:-1],
                                color='C0', alpha=0.3)
            ax.plot(x_axis[1:-1], self.dlVc[1:-1], label='dlVc')
            if hasattr(self, 'dlVc_SD'):
                ax.fill_between(x_axis[1:-1], self.dlVc[1:-1]-self.dlVc_SD[1:-1],
                                self.dlVc[1:-1]+self.dlVc_SD[1:-1],
                                color='C1', alpha=0.3)
            if hasattr(self, 'dlPolPeak'):
                ax.plot(x_axis[1:-1], self.dlPolPeak[1:-1], label='dl pol. peak')
                ax.fill_between(x_axis[1:-1],
                                self.dlPolPeak[1:-1]-self.dlPolPeak_SD[1:-1],
                                self.dlPolPeak[1:-1]+self.dlPolPeak_SD[1:-1],
                                color='C2', alpha=0.3)
            ax.set_yscale('log')
            ax.set_xlabel('Radius')
            ax.set_ylabel('Lengthscales')
            ax.set_xlim(x_axis.min(), x_axis.max())
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()
        elif self.name == 'fluxesR':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_axis, self.fcond, label='Fcond')
            if hasattr(self, 'fcond_SD'):
                ax.fill_between(x_axis, self.fcond-self.fcond_SD,
                                self.fcond+self.fcond_SD,
                                color='C0', alpha=0.3)
            ax.plot(x_axis, self.fconv, label='Fconv')
            if hasattr(self, 'fconv_SD'):
                ax.fill_between(x_axis, self.fconv-self.fconv_SD,
                                self.fconv+self.fconv_SD,
                                color='C1', alpha=0.3)
            ax.plot(x_axis, self.fkin, label='Fkin')
            if hasattr(self, 'fkin_SD'):
                ax.fill_between(x_axis, self.fkin-self.fkin_SD,
                                self.fkin+self.fkin_SD,
                                color='C2', alpha=0.3)
            ax.plot(x_axis, self.fvisc, label='Fvisc')
            if hasattr(self, 'fvisc_SD'):
                ax.fill_between(x_axis, self.fvisc-self.fvisc_SD,
                                self.fvisc+self.fvisc_SD,
                                color='C3', alpha=0.3)
            if self.prmag != 0:
                ax.plot(x_axis, self.fpoyn, label='Fpoyn')
                if hasattr(self, 'fpoyn_SD'):
                    ax.fill_between(x_axis, self.fpoyn-self.fpoyn_SD,
                                    self.fpoyn+self.fpoyn_SD,
                                    color='C4', alpha=0.3)
                ax.plot(x_axis, self.fres/self.fcond[0], label='Fres')
                if hasattr(self, 'fres_SD'):
                    ax.fill_between(x_axis, self.fres-self.fres_SD,
                                    self.fres+self.fres_SD,
                                    color='C5', alpha=0.3)

            ax.set_xlabel('Radius')
            ax.set_ylabel('Fluxes')
            ax.plot(x_axis, self.ftot, label='Ftot')
            ax.legend(loc='best', frameon=False)
            ax.set_xlim(x_axis[-1], x_axis[0])
            fig.tight_layout()
        elif self.name == 'heatR':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_axis, self.entropy, label='s', color='C0')
            if hasattr(self, 'entropy_SD'):
                ax.fill_between(x_axis, self.entropy-self.entropy_SD,
                                self.entropy+self.entropy_SD,
                                color='C0', alpha=0.3)
            if self.DissNb > 0:
                ax.plot(x_axis, self.temperature, label='T', color='C1')
                if hasattr(self, 'temperature_SD'):
                    ax.fill_between(x_axis, self.temperature-self.temperature_SD,
                                    self.temperature+self.temperature_SD,
                                    color='C1', alpha=0.3)
            if self.xi.max() > 1e-10:
                ax.plot(x_axis, self.xi, label='xi', color='C2')
                if hasattr(self, 'xi_SD'):
                    ax.fill_between(x_axis, self.xi-self.xi_SD,
                                    self.xi+self.xi_SD,
                                    color='C2', alpha=0.3)
            ax.set_xlabel('Radius')
            ax.set_ylabel('Temperature, Entropy, Composition')
            ax.legend(loc='best', frameon=False)
            ax.set_xlim(x_axis[-1], x_axis[0])
            fig.tight_layout()

            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            ax1.plot(x_axis, self.pressure, label='p')
            ax1.set_xlabel('Radius')
            ax1.set_ylabel('Pressure')
            ax1.set_xlim(x_axis[-1], x_axis[0])
            fig.tight_layout()
        elif self.name == 'perpParR':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_axis, self.Eperp, ls='-', color='C0', label='E perp.')
            if hasattr(self, 'Eperp_SD'):
                ax.fill_between(x_axis, self.Eperp-self.Eperp_SD,
                                self.Eperp+self.Eperp_SD,
                                color='C0', alpha=0.3)
            ax.plot(x_axis, self.Epar, ls='-', color='C1', label='E par.')
            if hasattr(self, 'Epar_SD'):
                ax.fill_between(x_axis, self.Epar-self.Epar_SD,
                                self.Epar+self.Epar_SD,
                                color='C1', alpha=0.3)
            ax.plot(x_axis, self.Eperp_axi, ls='--', color='C0', label='E eperp. ax.')
            if hasattr(self, 'Eperp_axi_SD'):
                ax.fill_between(x_axis, self.Eperp_axi-self.Eperp_axi_SD,
                                self.Eperp_axi+self.Eperp_axi_SD,
                                color='C0', alpha=0.3)
            ax.plot(x_axis, self.Epar_axi, ls='--', color='C1', label='E par. ax.')
            if hasattr(self, 'Epar_axi_SD'):
                ax.fill_between(x_axis, self.Epar_axi-self.Epar_axi_SD,
                                self.Epar_axi+self.Epar_axi_SD,
                                color='C1', alpha=0.3)
            ax.set_xlabel('Radius')
            ax.set_ylabel('Kinetic energy')
            ax.set_xlim(x_axis.min(), x_axis.max())
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()

        if hasattr(self, 'con_radratio'):
            if self.nVarCond == 2:
                ax.axvline(self.con_radratio*self.radius[0], color='k',
                          linestyle='--')


class RadLookUpTable:
    """
    The purpose of this class is to create a lookup table between the numpy
    array that comes from the reading of the radial file and the corresponding
    column.
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

        if self.name == 'eKinR':
            self.radius = data[:, 0]
            self.ekin_pol = data[:, 1]
            self.ekin_pol_axi = data[:, 2]
            self.ekin_tor = data[:, 3]
            self.ekin_tor_axi = data[:, 4]
            self.ekin_pol_surf = data[:, 5]
            self.ekin_pol_axi_surf = data[:, 6]
            self.ekin_tor_surf = data[:, 7]
            self.ekin_tor_axi_surf = data[:, 8]
            self.ekin_tot_r = self.ekin_pol + self.ekin_tor
        elif self.name == 'eMagR':
            self.radius = data[:, 0]
            self.emag_pol = data[:, 1]
            self.emag_pol_axi = data[:, 2]
            self.emag_tor = data[:, 3]
            self.emag_tor_axi = data[:, 4]
            self.emag_pol_surf = data[:, 5]
            self.emag_pol_axi_surf = data[:, 6]
            self.emag_tor_surf = data[:, 7]
            self.emag_tor_axi_surf = data[:, 8]
            self.dip_ratio = data[:, 9]
            self.emag_tot_r = self.emag_pol+self.emag_tor
            self.ecmb = self.emag_tot_r[0]
        elif self.name == 'anel':
            self.radius = data[:, 0]
            self.temp0 = data[:, 1]
            self.rho0 = data[:, 2]
            self.beta = data[:, 3]
            self.dbeta = data[:, 4]
            self.grav = data[:, 5]
            try:
                self.dsdr = data[:, 6]
            except IndexError:
                self.dsdr = np.zeros_like(self.radius)
            try:
                self.divkgradT = data[:, 7]
            except IndexError:
                self.divkgradT = np.zeros_like(self.radius)
            try:
                self.alpha0 = data[:, 8]
            except IndexError:
                self.alpha0 = np.zeros_like(self.radius)
            try:
                self.ogrun = data[:, 9]
            except IndexError:
                self.ogrun = np.zeros_like(self.radius)
            try:
                self.dLtemp0 = data[:, 10]
            except IndexError:
                self.dLtemp0 = np.zeros_like(self.radius)
        elif self.name == 'varDiff':
            self.radius = data[:, 0]
            self.conduc = data[:, 1]
            self.kappa = data[:, 2]
            self.dLkappa = data[:, 3]
            self.prandtl = data[:, 4]
        elif self.name == 'varVisc':
            self.radius = data[:, 0]
            self.dynVisc = data[:, 1]
            self.kinVisc = data[:, 2]
            self.dLvisc = data[:, 3]
            self.ekman = data[:, 4]
            self.prandtl = data[:, 5]
            self.prandtlmag = data[:, 6]
        elif self.name == 'varCond':
            self.radius = data[:, 0]
            self.conduc = data[:, 1]
            self.lmbda = data[:, 2]
        elif self.name == 'powerR':
            self.radius = data[:, 0]
            self.buoPower = data[:, 1]
            if data.shape[1] == 9 or data.shape[1] == 5:
                self.buoPower_chem = data[:, 2]
                self.viscDiss = data[:, 3]
                self.ohmDiss = data[:, 4]
                if data.shape[1] == 9:
                    self.buoPower_SD = data[:, 5]
                    self.buoPower_chem_SD = data[:, 6]
                    self.viscDiss_SD = data[:, 7]
                    self.ohmDiss_SD = data[:, 8]
            elif data.shape[1] == 7 or data.shape[1] == 4:
                self.viscDiss = data[:, 2]
                self.ohmDiss = data[:, 3]
                if data.shape[0] == 7:
                    self.buoPower_SD = data[:, 4]
                    self.viscDiss_SD = data[:, 5]
                    self.ohmDiss_SD = data[:, 6]
        elif self.name == 'parrad':
            self.radius = data[:, 0]
            self.rm = data[:, 1]
            self.rol = data[:, 2]
            self.urol = data[:, 3]
            self.dlV = data[:, 4]
            self.udlV = data[:, 5]
            self.udlVc = data[:, 8]
            self.entropy = data[:, 9]
            self.entropy_SD = np.sqrt(abs(data[:, 10]))
            self.uh = data[:, 11]
            self.duhdr = data[:, 12]
        elif self.name == 'parR':
            self.radius = data[:, 0]
            self.rm = data[:, 1]
            self.rol = data[:, 2]
            self.urol = data[:, 3]
            self.dlV = data[:, 4]
            self.dlVc = data[:, 5]
            if data.shape[1] == 8:
                self.udlV = data[:, 6]
                self.udlVc = data[:, 7]
            elif data.shape[1] == 13:
                self.dlPolPeak = data[:, 6]
                self.rm_SD = data[:, 7]
                self.rol_SD = data[:, 8]
                self.urol_SD = data[:, 9]
                self.dlV_SD = data[:, 10]
                self.dlVc_SD = data[:, 11]
                self.dlPolPeak_SD = data[:, 12]
        elif self.name == 'bLayersR':
            self.radius = data[:, 0]
            self.entropy = data[:, 1]
            if data.shape[1] == 5:
                self.entropy_SD = np.sqrt(abs(data[:, 2]))
                self.uh = data[:, 3]
                self.duhdr = data[:, 4]
            elif data.shape[1] == 6:
                self.entropy_SD = np.sqrt(abs(data[:, 2]))
                self.uh = data[:, 3]
                self.duhdr = data[:, 4]
                self.dissS = data[:, 5]
            elif data.shape[1] == 9:
                self.uh =data[:, 2]
                self.duhdr =data[:, 3]
                self.dissS =data[:, 4]
                self.entropy_SD =data[:, 5]
                self.uh_SD =data[:, 6]
                self.duhdr_SD =data[:, 7]
                self.dissS_SD =data[:, 8]
        elif self.name == 'fluxesR':
            self.radius = data[:, 0]
            self.fcond = data[:, 1]
            self.fconv = data[:, 2]
            self.fkin = data[:, 3]
            self.fvisc = data[:, 4]
            self.fpoyn = data[:, 5]
            self.fres = data[:, 6]
            self.ftot = self.fcond+self.fconv+self.fkin+self.fvisc+\
                        self.fpoyn+self.fres
            if data.shape[1] == 13:
                self.fcond_SD = data[:, 7]
                self.fconv_SD = data[:, 8]
                self.fkin_SD = data[:, 9]
                self.fvisc_SD = data[:, 10]
                self.fpoyn_SD = data[:, 11]
                self.fres_SD = data[:, 12]
        elif self.name == 'heatR':
            self.radius = data[:, 0]
            self.entropy = data[:, 1]
            self.temperature = data[:, 2]
            self.pressure = data[:, 3]
            self.density = data[:, 4]
            self.xi = data[:, 5]
            if data.shape[1] == 11:
                self.entropy_SD = data[:, 6]
                self.temperature_SD = data[:, 7]
                self.pressure_SD = data[:, 8]
                self.density_SD = data[:, 9]
                self.xi_SD = data[:, 10]
        elif self.name == 'perpParR':
            self.radius = data[:, 0]
            self.Eperp = data[:, 1]
            self.Epar = data[:, 2]
            self.Eperp_axi = data[:, 3]
            self.Epar_axi = data[:, 4]
            if data.shape[1] == 9:
                self.Eperp_SD = data[:, 5]
                self.Epar_SD = data[:, 6]
                self.Eperp_axi_SD = data[:, 7]
                self.Epar_axi_SD = data[:, 8]

    def __add__(self, new):
        """
        This method allows to sum two look up tables together. It is also
        working if the number of radial grid points has changed. In that
        case a spline interpolation is done to match the newest grid
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

        n_r_new = len(new.radius)
        n_r_old = len(self.radius)

        if n_r_new == n_r_old:
            for attr in new.__dict__.keys():
                if attr not in ['radius', 'name', 'start_time', 'stop_time']:
                    # Only stack if both new and old have the attribute available
                    if attr in self.__dict__:
                        # Standard deviation
                        if attr.endswith('SD'):
                            out.__dict__[attr] = np.sqrt(( \
                                                  fac_new*new.__dict__[attr]**2 + \
                                                  fac_old*self.__dict__[attr]**2) / \
                                                  fac_tot)
                        # Regular field
                        else:
                            out.__dict__[attr] = (fac_new*new.__dict__[attr] + \
                                                  fac_old*self.__dict__[attr]) / \
                                                  fac_tot
        else: # Different radial grid then interpolate on the new grid using splines
            rold = self.radius[::-1] # Splines need r increasing
            rnew = new.radius[::-1]
            for attr in new.__dict__.keys():
                if attr not in ['radius', 'name', 'start_time', 'stop_time']:
                    # Only stack if both new and old have the attribute available
                    if attr in self.__dict__:
                        datOldGrid = self.__dict__[attr]
                        tckp = sint.splrep(rold, datOldGrid[::-1])
                        datNewGrid = sint.splev(rnew, tckp)
                        self.__dict__[attr] = datNewGrid[::-1]
                        # Standard deviation
                        if attr.endswith('SD'):
                            out.__dict__[attr] = np.sqrt(( \
                                                  fac_new*new.__dict__[attr]**2 + \
                                                  fac_old*self.__dict__[attr]**2) / \
                                                  fac_tot)
                        # Regular field
                        else:
                            out.__dict__[attr] = (fac_new*new.__dict__[attr] + \
                                                  fac_old*self.__dict__[attr]) / \
                                                  fac_tot

        return out
