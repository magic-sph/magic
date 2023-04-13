#-*- coding: utf-8 -*-
import os
import re
import copy
import matplotlib.pyplot as plt
import numpy as np
from .log import MagicSetup
import glob
from .libmagic import fast_read, scanDir


class MagicTs(MagicSetup):
    """
    This python class is used to read and plot the different time series
    written by the code:

       * Kinetic energy: :ref:`e_kin.TAG <secEkinFile>`
       * Magnetic energy of the outer core: :ref:`e_mag_oc.TAG <secEmagocFile>`
       * Magnetic energy of the inner core: :ref:`e_mag_ic.TAG <secEmagicFile>`
       * Dipole information: :ref:`dipole.TAG <secDipoleFile>`
       * Rotation: :ref:`rot.TAG <secRotFile>`
       * Diagnostic parameters: :ref:`par.TAG <secParFile>`
       * Geostrophy: :ref:`geos.TAG <secGeosFile>`
       * Taylorization measures: :ref:`Tay.TAG <secTayFile>`
       * Heat transfer: :ref:`heat.TAG <secHeatFile>`
       * Helicity: :ref:`helicity.TAG <secHelicityFile>`
       * Velocity square: :ref:`u_square.TAG <secu_squareFile>`
       * Angular momentum: :ref:`AM.TAG <secAMFile>`
       * Power budget: :ref:`power.TAG <secpowerFile>`
       * Earth-likeness of the CMB field: :ref:`earth_like.TAG <secEarthLikeFile>`
       * Parallel and perpendicular decomposition: :ref:`perpPar.TAG <secperpParFile>`
       * Phase field: :ref:`phase.TAG <secphaseFile>`
       * Hemisphericity: :ref:`hemi.TAG <secHemiFile>`
       * RMS force balance: :ref:`dtVrms.TAG <secdtVrmsFile>`
       * RMS induction terms: :ref:`dtBrms.TAG <secdtBrmsFile>`
       * Time-evolution of m-spectra: :ref:`am_[kin|mag]_[pol|tor].TAG <secTimeSpectraFiles>`

    Here are a couple of examples of how to use this function.

    >>> # plot the most recent e_kin.TAG file found in the directoy
    >>> MagicTs(field='e_kin')
    >>>
    >>> # stack **all** the power.TAG file found in the directory
    >>> ts = MagicTs(field='power', all=True)
    >>> print(ts.time, ts.buoPower) # print time and buoyancy power
    >>>
    >>> # If you only want to read the file ``heat.N0m2z``
    >>> ts = MagicTs(field='heat', tag='N0m2z', iplot=False)
    """

    def __init__(self, datadir='.', field='e_kin', iplot=True, all=False, tag=None):
        """
        :param datadir: working directory
        :type datadir: str
        :param field: the file you want to plot
        :type field: str
        :param iplot: when set to True, display the plots (default True)
        :type iplot: bool
        :param all: when set to True, the complete time series is reconstructed by
                    stacking all the corresponding files from the working directory
                    (default False)
        :type all: bool
        :param tag: read the time series that exactly corresponds to the specified tag
        :type tag: str
        """
        self.field = field
        pattern = os.path.join(datadir, 'log.*')
        logFiles = scanDir(pattern)

        if self.field in ('am_mag_pol','am_mag_tor','am_kin_pol','am_kin_tor'):
            binary = True
        else:
            binary = False

        if tag is not None:
            pattern = os.path.join(datadir, '{}.{}'.format(self.field, tag))
            files = scanDir(pattern)

            # Either the log.tag directly exists and the setup is easy to obtain
            if os.path.exists(os.path.join(datadir, 'log.{}'.format(tag))):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.{}'.format(tag))
            # Or the tag is a bit more complicated and we need to find
            # the corresponding log file
            else:
                st = os.path.join(datadir, '{}\.(.*)'.format(self.field))
                mask = re.compile(st)
                if mask.match(files[-1]):
                    ending = mask.search(files[-1]).groups(0)[0]
                    pattern = os.path.join(datadir, 'log.{}'.format(ending))
                    if logFiles.__contains__(pattern):
                        MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml='log.{}'.format(ending))

            # Concatenate the files that correspond to the tag
            for k, file in enumerate(files):
                filename = file
                data = fast_read(filename, binary=binary)
                if k == 0:
                    tslut = TsLookUpTable(data, self.field)
                else:
                    tslut += TsLookUpTable(data, self.field)

        # If no tag is specified, the most recent is plotted
        elif not all:
            if len(logFiles) != 0:
                MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])
                name = '{}.{}'.format(self.field, self.tag)
                filename = os.path.join(datadir, name)
                data = fast_read(filename, binary=binary)
            else:
                mot = '{}.*'.format(self.field)
                dat = [(os.stat(i).st_mtime, i) for i in glob.glob(mot)]
                dat.sort()
                filename = dat[-1][1]
                data = fast_read(filename, binary=binary)
            tslut = TsLookUpTable(data, self.field)

        # If no tag is specified but all=True, all the directory is plotted
        else:
            if len(logFiles) != 0:
                MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])
            pattern = os.path.join(datadir, '{}.*'.format(self.field))
            files = scanDir(pattern)
            for k, file in enumerate(files):
                filename = file
                data = fast_read(filename, binary=binary)
                if len(data) > 0: # File is not empty
                    if k == 0:
                        tslut = TsLookUpTable(data, self.field)
                    else:
                        tslut += TsLookUpTable(data, self.field)

        try:
            # Copy look-up table arguments into MagicTs object
            for attr in tslut.__dict__:
                setattr(self, attr, tslut.__dict__[attr])

            if iplot:
                self.plot()
        except NameError: # In case tslut in not Defined
            print('No file correponding to field "{}" has been found'.format(self.field))

    def plot(self):
        """
        Plotting subroutines. Only called if 'iplot=True'
        """
        if self.field == 'e_kin':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.ekin_pol, ls='-', c='C0', label='ekin pol')
            ax.plot(self.time, self.ekin_tor, ls='-', c='C1', label='ekin tor')
            ax.plot(self.time, self.ekin_pol_axi, ls='--', c='C0',
                    label='ekin pol axi')
            ax.plot(self.time, self.ekin_tor_axi, ls='--', c='C1',
                    label='ekin tor axi')
            ax.plot(self.time, self.ekin_tot, ls='-', c='0.25', label='ekin tot')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Ekin')
            ax.set_yscale('log')
            ax.set_xlim(self.time[0], self.time[-1])
            fig.tight_layout()
        elif self.field == 'e_mag_oc':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.emagoc_pol, ls='-', c='C0', label='emag pol')
            ax.plot(self.time, self.emagoc_tor, ls='-', c='C1', label='emag tor')
            ax.plot(self.time, self.emagoc_pol_axi, ls='--', c='C0',
                    label='emag pol axi')
            ax.plot(self.time, self.emagoc_tor_axi, ls='--', c='C1',
                    label='emag tor axi')
            ax.plot(self.time, self.emag_tot, ls='-', c='0.25', label='emag tot')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Emag')
            ax.set_yscale('log')
            ax.set_xlim(self.time[0], self.time[-1])
            fig.tight_layout()

            # fig,ax = plt.subplots(1)
            # ax.plot(self.time, self.emag_es, ls='-',
            #         label=r'${E_B}^S$')
            # ax.plot(self.time, self.emag_eas, ls='-',
            #         label=r'${E_B}^A$')
            # ax.legend(loc='best', frameon=False)
            # ax.set_xlabel('Time')
            # ax.set_ylabel('Emag')

        elif self.field == 'e_mag_ic':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.emagic_pol, ls='-', c='C0', label='emagic pol')
            ax.plot(self.time, self.emagic_tor, ls='-', c='C1', label='emagic tor')
            ax.plot(self.time, self.emagic_pol_axi, ls='--', c='C0',
                    label='emagic pol axi')
            ax.plot(self.time, self.emagic_tor_axi, ls='--', c='C1',
                    label='emagic tor axi')
            ax.plot(self.time, self.emagic_pol+self.emagic_tor, ls='-', c='0.25',
                    label='emagic tot')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Emag inner core')
            ax.set_yscale('log')
            ax.set_xlim(self.time[0], self.time[-1])
            fig.tight_layout()
        elif self.field == 'rot':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.omega_ic, ls='-', label='Omega IC')
            ax.set_xlabel('Time')
            ax.set_ylabel('Rotation inner core')
            ax.legend(loc='best', frameon=False)
            ax.set_xlim(self.time[0], self.time[-1])
            fig.tight_layout()
            
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            ax1.plot(self.time, self.lorentz_torque_ic, ls='-', c='C0',
                     label='Lorentz torque on IC')
            ax1.plot(self.time,self.viscous_torque_ic, ls='-', c='C1',
                     label='Viscous torque on IC')
            ax1.legend(loc='best', frameon=False)
            ax1.set_xlabel('Time')
            ax1.set_ylabel('Torque on IC')
            ax1.set_xlim(self.time[0], self.time[-1])
            ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            ax1.axhline(0., color='0.5', ls='--', lw=1)
            fig1.tight_layout()
        elif self.field == 'timestep':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.step(self.time, self.dt)
            ax.set_yscale('log')
            ax.set_xlabel('Time')
            ax.set_ylabel('Time step size')
            fig.tight_layout()
        elif self.field == 'dipole':
            if self.ktopb != 2:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(self.time, self.theta_dip, label='theta_dip')
                #ax.plot(self.time, self.phi_dip, 'r-', label='phi_dip')
                ax.set_ylabel('Dipole angle')
                ax.set_xlabel('Time')
                ax.set_ylim(-1., 181)
                ax.set_xlim(self.time[0], self.time[-1])
                fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.dipTot, label='Total dipolarity')
            ax.plot(self.time, self.dipolarity, ls='--', label='Axisym dipolarity')
            ax.plot(self.time, self.dipTot_cmb, ls='-', c='C2',
                    label='Total dipolarity CMB')
            ax.plot(self.time, self.dip_cmb, ls='--', c='C2',
                    label='Axisym dipolarity')
            if hasattr(self, 'l_geo'):
                lcut = self.l_geo
            else:
                lcut = 11
            ax.plot(self.time, self.dip_l11, ls='-', c='C3',
                    label='Axisym dip l={:d}'.format(lcut))
            ax.plot(self.time, self.dipTot_l11, ls='--', c='C3',
                    label='Total dip l={:d}'.format(lcut))
            # ax.plot(self.time, self.dip3, ls='-', c='#e5ae38',
            #         label='Epol axi/Ecmb')
            ax.legend(loc='best', frameon=False)
            ax.set_ylabel('Dipolarity')
            ax.set_xlabel('Time')
            ax.set_ylim(0, 1)
            ax.set_xlim(self.time[0], self.time[-1])
            fig.tight_layout()
        elif self.field == 'AM':
            fig = plt.figure()
            ax = fig.add_subplot(211)
            ax.plot(self.time, self.am_oc_z, label='Outer core')
            ax.plot(self.time, self.am_ic, label='Inner core')
            ax.plot(self.time, self.amz, ls='-', c='0.25', label='Total')
            ax.set_xlim(self.time[0], self.time[-1])
            ax.legend(loc='best', frameon=False)
            ax.set_ylabel('AM')
            ax = fig.add_subplot(212)
            ax.semilogy(self.time[1:], np.abs(self.damzdt[1:]))
            ax.set_xlabel('Time')
            ax.set_ylabel('dAmz / dt')
            ax.set_xlim(self.time[1], self.time[-1])
            fig.tight_layout()
        elif self.field == 'par':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            if self.mode == 1 or self.prmag == 0.:
                ax.semilogy(self.time, self.rm, label='Reynolds')
            else:
                ax.semilogy(self.time, self.rm, label='Magnetic Reynolds')
            if self.elsasser.max() > 0.:
                ax.semilogy(self.time, self.elsasser, label='Elsasser')
                ax.semilogy(self.time, self.els_cmb, label='Elsasser CMB')
            ax.semilogy(self.time, self.rossby_l, label='Rossby l')
            if hasattr(self, 'rolc'):
                ax.semilogy(self.time, self.rolc, label='Roc l')
            ax.legend(loc='lower right', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Params')
            ax.set_xlim(self.time[0], self.time[-1])
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(self.time, self.dlV, label='Integral (ell)')
            ax.semilogy(self.time, self.dlVc, label='Integral (ell c)')
            ax.semilogy(self.time, self.dmV, label='Integral (m)')
            ax.semilogy(self.time, self.dlPolPeak, label='Peak (pol)')
            if abs(self.lbDiss).max() > 0.:
                ax.semilogy(self.time, self.lbDiss, label='Magnetic dissipation')
            if abs(self.lvDiss).max() > 0.:
                ax.semilogy(self.time, self.lvDiss, label='Viscous dissipation')
            ax.set_xlabel('Time')
            ax.set_ylabel('Lengthscales')
            ax.set_xlim(self.time[0], self.time[-1])
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()

            if self.dipolarity.max() > 0.:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(self.time, self.dipolarity, label='Dipolarity')
                ax.plot(self.time, self.dip_cmb, label='Dipolarity CMB')
                ax.legend(loc='upper right', frameon=False)
                ax.set_xlim(self.time[0], self.time[-1])
                ax.set_xlabel('Time')
                ax.set_ylabel('Dipolarity')
                ax.set_ylim(0, 1)
                fig.tight_layout()
        elif self.field == 'geos':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.geos, label='Total')
            ax.plot(self.time, self.geosM, label='Meridional')
            ax.plot(self.time, self.geosZ, label='Zonal')
            ax.plot(self.time, self.geosNAP, label='Non-axi perp')
            ax.set_xlabel('Time')
            ax.set_ylabel('Geostrophy')
            ax.set_xlim(self.time[0], self.time[-1])
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.corr_vz_otc, label='uz')
            ax.plot(self.time, self.corr_vortz_otc, label='z vorticity')
            ax.plot(self.time, self.corr_hel_otc, label='Helicity')
            ax.set_xlabel('Time')
            ax.set_ylabel('z correlations')
            ax.set_xlim(self.time[0], self.time[-1])
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()
        elif self.field == 'phase':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax1 = ax.twinx()
            ax.plot(self.time, self.rmelt, label='r melt', color='C0')
            ax1.plot(self.time, self.trmelt, label='T(r melt)', color='C1')
            ax.set_xlim(self.time[0], self.time[-1])
            ax.set_xlabel('Time')
            ax.set_ylabel('r melt')
            ax1.set_ylabel('T(r melt)')
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(self.time, self.ekinS/self.ekinL)
            ax.set_xlim(self.time[0], self.time[-1])
            ax.set_xlabel('Time')
            ax.set_ylabel('Relative energy fraction in solidus')
            fig.tight_layout()
        elif self.field == 'hemi':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.hemi_emag, label='Emag')
            ax.plot(self.time, self.hemi_br, label='|Br| volume')
            ax.plot(self.time, self.hemi_cmb, label='|Br| CMB')
            ax.set_xlabel('Time')
            ax.set_ylabel('Hemisphericity')
            ax.set_xlim(self.time[0], self.time[-1])
            ax.set_ylim(0., 1.)
            #ax.set_xlim(self.time[0], self.time[-1])
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()
        elif self.field == 'earth_like':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.axial_dipole, label='AD/NAD')
            ax.plot(self.time, self.symmetry, label='O/E')
            ax.plot(self.time, self.zonality, label='Z/NZ')
            ax.plot(self.time, self.flux_concentration, label='FCF')
            ax.set_xlim(self.time[0], self.time[-1])
            ax.set_xlabel('Time')
            ax.set_ylabel('Rating parameters')
            ax.legend(loc='upper right', frameon=False)
            fig.tight_layout()

            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            ax1.plot(self.time, self.chi_square)
            ax1.set_xlim(self.time[0], self.time[-1])
            ax1.set_xlabel('Time')
            ax1.set_ylabel('Chi square')
            fig1.tight_layout()
        elif self.field == 'misc':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.topnuss, label='Top Nusselt')
            ax.plot(self.time, self.botnuss, label='Bottom Nusselt')
            ax.legend(loc='lower right', frameon=False)

            ax.set_xlim(self.time[0], self.time[-1])
            ax.set_xlabel('Time')
            ax.set_ylabel('Nusselt number')
            fig.tight_layout()
            if self.helrms.max() != 0.:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(self.time, self.helrms)
                ax.set_xlim(self.time[0], self.time[-1])
                ax.set_xlabel('Time')
                ax.set_ylabel('Helicity')
                fig.tight_layout()
        elif self.field == 'heat':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            if self.kbots == 2 and self.ktops == 2:
                ax.plot(self.time, self.deltaTnuss, label=r'$Nu_{\Delta T}$')
            else:
                ax.plot(self.time, self.topnuss, label='Top Nusselt')
                ax.plot(self.time, self.botnuss, label='Bottom Nusselt')
            ax.set_xlim(self.time[0], self.time[-1])
            ax.legend(loc='lower right', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Nusselt number')
            fig.tight_layout()

            if self.topsherwood.max() != 1.0 or self.deltasherwood.max() != 1.0:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                if self.kbotxi == 2 and self.ktopxi == 2:
                    ax.plot(self.time, self.deltasherwood, label=r'$Sh_{\Delta \xi}$')
                else:
                    ax.plot(self.time, self.topsherwood, label='Top Sherwood')
                    ax.plot(self.time, self.botsherwood, label='Bottom Sherwood')
                ax.legend(loc='lower right', frameon=False)
                ax.set_xlim(self.time[0], self.time[-1])
                ax.set_xlabel('Time')
                ax.set_ylabel('Sherwood number')
                fig.tight_layout()
        elif self.field == 'helicity':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.helRMSN, label='Northern Hemisphere')
            ax.plot(self.time, self.helRMSS, label='Southern Hemisphere')
            ax.set_xlim(self.time[0], self.time[-1])
            ax.legend(loc='lower right', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Helicity')
            fig.tight_layout()
        elif self.field == 'u_square':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.ekin_pol, ls='-', c='C0', label='ekin pol')
            ax.plot(self.time, self.ekin_tor, ls='-', c='C1', label='ekin tor')
            ax.plot(self.time, self.ekin_pol_axi, ls='--', c='C0',
                    label='ekin pol axi')
            ax.plot(self.time, self.ekin_tor_axi, ls='--', c='C1',
                    label='ekin tor axi')
            ax.plot(self.time, self.ekin_tot, ls='-', c='0.25', label='ekin tot')
            ax.set_xlim(self.time[0], self.time[-1])
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('u**2')
            ax.set_yscale('log')
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            if self.mode == 1 or self.prmag == 0.:
                ax.semilogy(self.time, self.rm, label='Reynolds')
            else:
                ax.semilogy(self.time, self.rm, label='Magnetic Reynolds')
            ax.semilogy(self.time, self.ro, label='Rossby')
            ax.semilogy(self.time, self.rossby_l, label='Rossby l')
            ax.semilogy(self.time, self.dl, label='l')
            ax.set_xlim(self.time[0], self.time[-1])
            ax.legend(loc='lower right', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Params')
            fig.tight_layout()
        elif self.field in ('dtVrms'):
            fig = plt.figure() # Poloidal forces
            ax = fig.add_subplot(111)
            ax.semilogy(self.time, self.CorRms, label='Coriolis')
            ax.semilogy(self.time, self.PreRms, label='Pressure')
            ax.semilogy(self.time, self.LFRms, label='Lorentz')
            ax.semilogy(self.time, self.BuoRms, label='Thermal Buoyancy')
            if abs(self.ChemRms).max() > 0:
                ax.semilogy(self.time, self.ChemRms, label='Chemical Buoyancy')
            ax.semilogy(self.time, self.InerRms, label='Inertia')
            ax.semilogy(self.time, self.DifRms, label='Diffusion')

            ax.set_xlim(self.time[0], self.time[-1])
            ax.legend(loc='best', frameon=False, ncol=2)
            ax.set_xlabel('Time')
            ax.set_ylabel('RMS forces')
            fig.tight_layout()

            fig = plt.figure() # Toroidal forces
            ax = fig.add_subplot(111)
            ax.semilogy(self.time, self.geos, label='Geostrophic balance')
            ax.semilogy(self.time, self.mageos, label='Magnetostrophic')
            ax.semilogy(self.time, self.arc, label='Archimedean')
            ax.semilogy(self.time, self.arcMag, label='Archimedean+Lorentz')
            ax.semilogy(self.time, self.corLor, label='Coriolis/Lorentz')
            ax.semilogy(self.time, self.preLor, label='Pressure/Lorentz')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('RMS balances')
            ax.set_xlim(self.time[0], self.time[-1])
            fig.tight_layout()
        elif self.field == 'perpPar':
            fig = plt.figure()
            ax= fig.add_subplot(111)
            ax.plot(self.time, self.eperp, ls='-', c='C0', label='ekin perp')
            ax.plot(self.time, self.epar, ls='-', c='C1', label='ekin par')
            ax.plot(self.time, self.eperp_axi, ls='--', c='C0',
                    label='ekin perp axi')
            ax.plot(self.time, self.epar_axi, ls='--', c='C1',
                    label='ekin par axi')
            ax.plot(self.time, self.ekin_tot, ls='-', c='0.25', label='ekin tot')
            ax.plot(self.time, self.ekin_tot, 'k-')
            ax.set_xlim(self.time[0], self.time[-1])
            ax.set_yscale('log')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Kinetic energy')

            ax.set_xlim(self.time[0], self.time[-1])
            fig.tight_layout()
        elif self.field in ('power'):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            if self.buoPower.max() != 0.:
               ax.semilogy(self.time, self.buoPower, label='Thermal buoyancy')
            if self.buoPower_chem.max() != 0.:
                ax.semilogy(self.time, self.buoPower_chem,
                            label='Chemical buoyancy')
            if self.ohmDiss.max() != 0.:
                ax.semilogy(self.time, -self.ohmDiss, label='Ohmic diss.')
            ax.semilogy(self.time, -self.viscDiss, label='Viscous diss.')
            ax.set_xlim(self.time[0], self.time[-1])
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Power')
            fig.tight_layout()

            if hasattr(self,'fohm'):
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(self.time, self.fohm)
                ax.set_xlim(self.time[0], self.time[-1])
                ax.set_ylim(0., 1.)
                ax.set_xlabel('Time')
                ax.set_ylabel('fohm')
                fig.tight_layout()
        elif self.field in ('dtBrms'):
            fig = plt.figure() # Poloidal
            ax = fig.add_subplot(111)
            ax.semilogy(self.time, self.DynPolRms, label='Induction')
            ax.semilogy(self.time, self.DifPolRms, label='Diffusion')
            ax.semilogy(self.time, self.dtBpolRms, label='Time derivative')

            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Poloidal field production')
            fig.tight_layout()

            fig = plt.figure() # Toroidal
            ax = fig.add_subplot(111)
            ax.semilogy(self.time, self.DynTorRms, label='Induction')
            ax.semilogy(self.time, self.DifTorRms, label='Diffusion')
            ax.semilogy(self.time, self.omEffect*self.DynTorRms,
                        label='Omega effect')
            ax.semilogy(self.time, self.dtBtorRms, label='Time derivative', )
            ax.set_xlim(self.time[0], self.time[-1])
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Toroidal field production')
            fig.tight_layout()
        elif self.field in ('SRIC'):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(self.time, self.viscTorq, label='Viscous')
            ax.semilogy(self.time, self.LorTorq, label='Lorentz')
            ax.semilogy(self.time, self.totTorq, label='Total')

            ax.set_xlim(self.time[0], self.time[-1])
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Torque')
            fig.tight_layout()
        elif self.field in ('am_mag_pol', 'am_mag_tor', 'am_kin_pol', 'am_kin_tor'):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for k in range(self.coeffs.shape[1]):
                ax.semilogy(self.time, self.coeffs[:, k], label='m={}'.format(k))
            ax.set_xlabel('Time')
            if self.coeffs.shape[1] < 20:
                ax.legend(loc='best', frameon=False)
            if self.field == 'am_mag_pol':
                ax.set_ylabel('Emag poloidal')
            elif self.field == 'am_mag_tor':
                ax.set_ylabel('Emag toroidal')
            elif self.field == 'am_kin_pol':
                ax.set_ylabel('Ekin poloidal')
            elif self.field == 'am_kin_tor':
                ax.set_ylabel('Ekin toroidal')
            ax.set_xlim(self.time[0], self.time[-1])
            fig.tight_layout()


class TsLookUpTable:
    """
    The purpose of this class is to create a lookup table between the numpy
    array that comes from the reading of the time series and the corresponding
    column.
    """

    def __init__(self, data, field):
        """
        :param data: numpy array that contains the data
        :type data: numpy.ndarray
        :param field: name of the field (i.e. 'eKinR', 'eMagR', 'powerR', ...)
        :type field: str
        """

        self.field = field

        if self.field == 'e_kin':
            self.time = data[:, 0]
            self.ekin_pol = data[:, 1]
            self.ekin_tor = data[:, 2]
            self.ekin_pol_axi = data[:, 3]
            self.ekin_tor_axi = data[:, 4]
            self.ekin_pol_es = data[:, 5]
            self.ekin_tor_es = data[:, 6]
            self.ekin_pol_es_axi = data[:, 7]
            self.ekin_tor_es_axi = data[:, 8]
            self.ekin_tot = self.ekin_pol + self.ekin_tor
            self.ekin_axi = self.ekin_pol_axi + self.ekin_tor_axi
            self.ekin_es = self.ekin_pol_es + self.ekin_tor_es
            self.ekin_es_axi = self.ekin_pol_es_axi + self.ekin_tor_es_axi
            self.ekin_pol_naxi = self.ekin_pol-self.ekin_pol_axi
            self.ekin_tor_naxi = self.ekin_tor-self.ekin_tor_axi
            self.ekin_es_naxi = self.ekin_es-self.ekin_es_axi
            self.ekin_naxi = self.ekin_tot-self.ekin_pol_axi-self.ekin_tor_axi
        elif self.field == 'e_mag_oc':
            self.time = data[:, 0]
            self.emagoc_pol = data[:, 1]
            self.emagoc_tor = data[:, 2]
            self.emagoc_pol_axi = data[:, 3]
            self.emagoc_tor_axi = data[:, 4]
            self.ext_nrj_pol = data[:, 5]
            self.ext_nrj_pol_axi = data[:, 6]
            self.emagoc_pol_eas = data[:, 7]
            self.emagoc_tor_eas = data[:, 8]
            self.emagoc_pol_eas_axi = data[:, 9]
            self.emagoc_tor_eas_axi = data[:, 10]
            self.emag_tot = self.emagoc_pol + self.emagoc_tor
            self.emag_axi = self.emagoc_pol_axi + self.emagoc_tor_axi
            self.emag_eas = self.emagoc_pol_eas + self.emagoc_tor_eas
            self.emag_eas_axi = self.emagoc_pol_eas_axi + self.emagoc_tor_eas_axi
        elif self.field == 'e_mag_ic':
            self.time = data[:, 0]
            self.emagic_pol = data[:, 1]
            self.emagic_tor = data[:, 2]
            self.emagic_pol_axi = data[:, 3]
            self.emagic_tor_axi = data[:, 4]
            self.emagic_tot = self.emagic_pol + self.emagic_tor
        elif self.field == 'timestep':
            self.time = data[:, 0]
            self.dt = data[:, 1]
        elif self.field == 'dipole':
            self.time = data[:, 0]
            self.theta_dip = data[:, 1]
            self.phi_dip = data[:, 2]
            self.dipolarity = data[:, 3]
            self.dip_cmb = data[:, 4]
            self.dip_l11 = data[:, 5] # Cut at l=11
            self.dipTot = data[:, 6] # Also non axisymmetric dipole
            self.dipTot_cmb = data[:, 7] # Non-axi at the CMB
            self.dipTot_l11 = data[:, 8] # Cut at l=11
            self.e_dip_cmb = data[:, 9]
            self.e_dip_ax_cmb = data[:, 10]
            self.e_dip = data[:, 11]
            self.e_dip_ax = data[:, 12]
            self.ecmb = data[:, 13]
            self.egeo = data[:, 14]
            self.ratio_cmb_as = data[:, 16] # (e_cmb-e_es_cmb)/e_cmb
            self.ratio_cmb_naxi = data[:, 17] # (e_cmb-e_axi_cmb)/e_cmb
            self.ratio_l11_cmb_as = data[:, 18] # (e_geo-e_es_geo)/e_geo
            self.ratio_l11_cmb_naxi = data[:, 19] # (e_geo-e_axi_geo)/e_geo
            self.epol_axi_cmb = (-self.ratio_cmb_naxi*self.ecmb+self.ecmb)
            self.epol_rel_cmb = self.epol_axi_cmb/self.ecmb
            self.fdip = np.sqrt(self.dip_l11)
        elif self.field == 'AM':
            self.time = data[:, 0]
            self.am_oc_x = data[:, 1]
            self.am_oc_y = data[:, 2]
            self.am_oc_z = data[:, 3]
            self.am_ic = data[:, 4]
            self.am_ma = data[:, 5]
            self.amz = data[:, 6]
            self.damzdt = data[:, 7]
        elif self.field == 'rot':
            self.time = data[:, 0]
            self.omega_ic = data[:, 1]
            self.lorentz_torque_ic = data[:, 2]
            self.viscous_torque_ic = data[:, 3]
            self.omega_ma = data[:, 4]
            self.lorentz_torque_ma = data[:, 5]
            self.viscous_torque_ma = data[:, 6]
        elif self.field == 'Tay':
            self.time = data[:, 0]
            self.ekin_tora_rel = data[:, 1]
            self.egeos_rel = data[:, 2]
            self.tay = data[:, 3]
            self.tayR = data[:, 4]
            self.tayV = data[:, 5]
            self.ekin_cyl = data[:, 6]
        elif self.field == 'par':
            self.time = data[:, 0]
            self.rm = data[:, 1]
            self.elsasser = data[:, 2]
            self.rossby_l = data[:, 3]
            self.geos = data[:, 4]
            self.dipolarity = data[:, 5]
            self.dip_cmb = data[:, 6]
            self.dlV = data[:, 7]
            self.dmV = data[:, 8]
            self.lvDiss = data[:, 11]
            self.lbDiss = data[:, 12]
            self.dlB = data[:, 13]
            self.dmB = data[:, 14]
            self.els_cmb = data[:, 15]
            if data.shape[-1] > 16:
                self.rolc = data[:, 16]
                self.dlVc = data[:, 17]
                self.reEquat = data[:, 18]
                self.dlPolPeak = np.zeros_like(self.time)
                if data.shape[-1] == 20:
                    self.dlPolPeak = data[:, 18]
                    self.reEquat = data[:, 19]
            else:
                self.rolc = np.zeros_like(self.time)
                self.dlVc = np.zeros_like(self.time)
                self.reEquat = np.zeros_like(self.time)
                self.dlPolPeak = np.zeros_like(self.time)
        elif self.field == 'misc':
            self.time = data[:, 0]
            self.botnuss = data[:, 1]
            self.topnuss = data[:, 2]
            self.bottemp = data[:, 3] / np.sqrt(4.*np.pi)
            self.toptemp = data[:, 4] / np.sqrt(4.*np.pi)
            self.helrms = data[:, 8]
            self.helN = data[:, 5]*self.helrms
            self.helS = data[:, 6]*self.helrms
            try:
                self.botflux = data[:, 16]
                self.topflux = data[:, 17]
            except IndexError:
                self.botflux = np.zeros_like(self.time)
                self.topflux = np.zeros_like(self.time)
                pass
        elif self.field == 'geos':
            self.time = data[:, 0]
            self.geos = data[:, 1]
            self.ekin_ntc_rel = data[:, 2]
            self.ekin_stc_rel = data[:, 3]
            self.ekin = data[:, 4]
            self.corr_vz_otc = data[:, 5]
            self.corr_vortz_otc = data[:, 6]
            self.corr_hel_otc = data[:, 7]
            if data.shape[-1] == 8:
                self.geosA= np.zeros_like(self.time)
                self.geosZ= np.zeros_like(self.time)
                self.geosM= np.zeros_like(self.time)
                self.geosNAP = np.zeros_like(self.time)
            elif data.shape[-1] == 12:
                self.geosA = data[:, 8]
                self.geosZ = data[:, 9]
                self.geosM = data[:, 10]
                self.geosNAP = data[:, 11]
        elif self.field == 'heat':
            self.time = data[:, 0]
            self.botnuss = data[:, 1]
            self.topnuss = data[:, 2]
            self.deltaTnuss = data[:, 3]
            self.bottemp = data[:, 4]
            self.toptemp = data[:, 5]
            self.bots = data[:, 6]
            self.tops = data[:, 7]
            self.topflux = data[:, 8]
            self.botflux = data[:, 9]
            self.toppress = data[:, 10]
            self.mass = data[:, 11]
            try:
                self.botsherwood = data[:, 12]
                self.topsherwood = data[:, 13]
                self.deltasherwood = data[:, 14]
                self.botxi = data[:, 15]
                self.topxi = data[:, 16]
            except IndexError:
                self.topsherwood = np.ones_like(self.time)
                self.botsherwood = np.ones_like(self.time)
                self.deltasherwood = np.ones_like(self.time)
                self.botxi = np.zeros_like(self.time)
                self.topxi = np.zeros_like(self.time)
                pass
        elif self.field == 'helicity':
            self.time = data[:, 0]
            self.helN = data[:, 1]
            self.helS = data[:, 2]
            self.helRMSN = data[:, 3]
            self.helRMSS = data[:, 4]
            self.helnaN = data[:, 5]
            self.helnaS = data[:, 6]
            self.helnaRMSN = data[:, 7]
            self.helnaRMSS = data[:, 8]
        elif self.field == 'earth_like':
            self.time = data[:, 0]
            self.axial_dipole = data[:, 1]
            self.symmetry = data[:, 2]
            self.zonality = data[:, 3]
            self.flux_concentration = data[:, 4]
            self.chi_square = ((np.log(self.axial_dipole)-np.log(1.4))/np.log(2.))**2+\
                              ((np.log(self.symmetry)-np.log(1.))/np.log(2.))**2+\
                              ((np.log(self.zonality)-np.log(0.15))/np.log(2.5))**2+\
                              ((np.log(self.flux_concentration)-np.log(1.5))/np.log(1.75))**2
        elif self.field == 'u_square':
            self.time = data[:, 0]
            self.ekin_pol = data[:, 1]
            self.ekin_tor = data[:, 2]
            self.ekin_pol_axi = data[:, 3]
            self.ekin_tor_axi = data[:, 4]
            self.ekin_tot = self.ekin_pol + self.ekin_tor
            self.ro = data[:, 5]
            self.rm = data[:, 6]
            self.rossby_l = data[:, 7]
            self.dl = data[:, 8]
        elif self.field == 'perpPar':
            self.time = data[:, 0]
            self.eperp = data[:, 1]
            self.epar = data[:, 2]
            self.eperp_axi = data[:, 3]
            self.epar_axi = data[:, 4]
            self.ekin_tot = self.eperp+self.epar
        elif self.field == 'phase':
            self.time = data[:, 0]
            self.rmelt = data[:, 1]
            self.trmelt = data[:, 2]
            self.volS = data[:, 3]
            self.ekinS = data[:, 4]
            self.ekinL = data[:, 5]
            self.flux_cmb = data[:, 6]
            self.flux_icb = data[:, 7]
            self.dEnthdt = data[:, 8]
        elif self.field == 'hemi':
            self.time = data[:, 0]
            self.hemi_vr = data[:, 1]
            self.hemi_ekin = data[:, 2]
            self.hemi_br = data[:, 3]
            self.hemi_emag = data[:, 4]
            self.hemi_cmb = data[:, 5]
            self.ekin = data[:, 6]
            self.emag = data[:, 7]
        elif self.field == 'dtVrms':
            self.time = data[:, 0]
            self.InerRms = data[:, 1]
            self.CorRms = data[:, 2]
            self.LFRms = data[:, 3]
            self.AdvRms = data[:, 4]
            self.DifRms = data[:, 5]
            self.BuoRms = data[:, 6]

            if data.shape[1] == 14:
                self.PreRms = data[:, 7]
                self.geos = data[:, 8] # geostrophic balance
                self.mageos = data[:, 9] # magnetostrophic balance
                self.arcMag = data[:, 10] # Coriolis/Pressure/Buoyancy/Lorentz
                self.corLor = data[:, 11] # Coriolis/Lorentz
                self.preLor = data[:, 12] # Pressure/Lorentz
                self.cia = data[:, 13] # Coriolis/Inertia/Archmedean
                self.arc = np.zeros_like(self.geos)
                self.ChemRms = np.zeros_like(self.geos)
            elif data.shape[1] == 15:
                self.PreRms = data[:, 7]
                self.geos = data[:, 8] # geostrophic balance
                self.mageos = data[:, 9] # magnetostrophic balance
                self.arc    = data[:, 10] # Coriolis/Pressure/Buoyancy
                self.arcMag = data[:, 11] # Coriolis/Pressure/Buoyancy/Lorentz
                self.corLor = data[:, 12] # Coriolis/Lorentz
                self.preLor = data[:, 13] # Pressure/Lorentz
                self.cia = data[:, 14] # Coriolis/Inertia/Archmedean
                self.ChemRms = np.zeros_like(self.geos)
            elif data.shape[1] == 16:
                self.ChemRms = data[:, 7]
                self.PreRms = data[:, 8]
                self.geos = data[:, 9] # geostrophic balance
                self.mageos = data[:, 10] # magnetostrophic balance
                self.arc    = data[:, 11] # Coriolis/Pressure/Buoyancy
                self.arcMag = data[:, 12] # Coriolis/Pressure/Buoyancy/Lorentz
                self.corLor = data[:, 13] # Coriolis/Lorentz
                self.preLor = data[:, 14] # Pressure/Lorentz
                self.cia = data[:, 15] # Coriolis/Inertia/Archmedean
            elif data.shape[1] == 18:
                self.ChemRms = data[:, 7]
                self.PreRms = data[:, 8]
                self.MagTensRms = data[:, 9] # Magnetic tension
                self.MagPreRms = data[:, 10] # Magnetic pressure
                self.geos = data[:, 11] # geostrophic balance
                self.mageos = data[:, 12] # magnetostrophic balance
                self.arc    = data[:, 13] # Coriolis/Pressure/Buoyancy
                self.arcMag = data[:, 14] # Coriolis/Pressure/Buoyancy/Lorentz
                self.corLor = data[:, 15] # Coriolis/Lorentz
                self.preLor = data[:, 16] # Pressure/Lorentz
                self.cia = data[:, 17] # Coriolis/Inertia/Archmedean

        elif self.field == 'dtBrms':
            self.time = data[:, 0]
            self.dtBpolRms = data[:, 1]
            self.dtBtorRms = data[:, 2]
            self.DynPolRms = data[:, 3]
            self.DynTorRms = data[:, 4]
            self.DifPolRms = data[:, 5]
            self.DifTorRms = data[:, 6]
            self.omEffect = data[:, 7]
            self.omega = data[:, 8]
            self.DynDipRms = data[:, 9]
            self.DynDipAxRms = data[:, 10]
        elif self.field == 'dtE':
            self.time = data[:, 0]
            self.dEdt = data[:, 1]
            self.intdEdt = data[:, 2]
            self.reldEdt = data[:, 3]
        elif self.field == 'power':
            self.time = data[:, 0]
            self.buoPower = data[:, 1]
            if data.shape[1] == 11:
                self.buoPower_chem = data[:, 2]
                self.icrotPower = data[:, 3]
                self.mantlerotPower = data[:, 4]
                self.viscDiss = data[:, 5]
                self.ohmDiss = data[:, 6]
                self.icPower = data[:, 7]
                self.mantlePower = data[:, 8]
            elif data.shape[1] == 10:
                self.buoPower_chem = np.zeros_like(self.time)
                self.icrotPower = data[:, 2]
                self.mantlerotPower = data[:, 3]
                self.viscDiss = data[:, 4]
                self.ohmDiss = data[:, 5]
                self.icPower = data[:, 6]
                self.mantlePower = data[:, 7]
            if abs(self.ohmDiss).max() != 0:
                 self.fohm = -self.ohmDiss/(self.buoPower+self.buoPower_chem)
                 self.fvis = -self.viscDiss/(self.buoPower+self.buoPower_chem)
        elif self.field == 'SRIC':
            self.time = data[:,0]
            self.omega_ic = data[:,1]
            self.viscPower = data[:,2]
            self.totPower = data[:,3]
            self.LorPower = data[:,4]
            self.viscTorq = abs(self.viscPower/self.omega_ic)
            self.totTorq = abs(self.totPower/self.omega_ic)
            self.LorTorq = abs(self.LorPower/self.omega_ic)
        elif self.field in ('am_mag_pol', 'am_mag_tor', # Tayler instability
                            'am_kin_pol', 'am_kin_tor'):
            self.time = data[:, 0]
            self.coeffs = data[:, 1:]
        else:
            print('The field "{}" is not know'.format(self.field))

    def __add__(self, new):
        """
        This method allows to sum two look up tables together. This is python
        built-in method.
        """

        out = copy.deepcopy(new)
        timeOld = self.time[-1]
        timeNew = new.time[0]

        for attr in new.__dict__.keys():
            if attr == 'coeffs':
                out.__dict__[attr] = np.vstack((self.__dict__[attr],
                                                out.__dict__[attr][1:, :]))
            elif attr != 'field':
                if attr in self.__dict__.keys():  # If the argument already existed
                    if timeOld != timeNew:
                            out.__dict__[attr] = np.hstack((self.__dict__[attr],
                                                            out.__dict__[attr]))
                    else: # Same time
                        out.__dict__[attr] = np.hstack((self.__dict__[attr],
                                                        out.__dict__[attr][1:]))

        return out
