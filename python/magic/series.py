#-*- coding: utf-8 -*-
import os
import re
import copy
import matplotlib.pyplot as plt
import numpy as np
from .log import MagicSetup
import glob
from .libmagic import (fast_read, scanDir, avgField,
                       timeder,secondtimeder, ReadBinaryTimeseries)


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
            # Copy look-up table arguments into MagicRadial object
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
            ax.plot(self.time, self.ekin_pol, ls='-', c='#30a2da',
                    label='ekin pol')
            ax.plot(self.time, self.ekin_tor, ls='-', c='#fc4f30',
                    label='ekin tor')
            ax.plot(self.time, self.ekin_pol_axi, ls='--', c='#30a2da',
                    label='ekin pol axi')
            ax.plot(self.time, self.ekin_tor_axi, ls='--', c='#fc4f30',
                 label='ekin tor axi')
            ax.plot(self.time, self.ekin_tot, ls='-', c='#31363B')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Ekin')
            fig.tight_layout()
        elif self.field == 'e_mag_oc':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.emagoc_pol, ls='-', c='#30a2da',
                    label='emag pol')
            ax.plot(self.time, self.emagoc_tor, ls='-', c='#fc4f30',
                    label='emag tor')
            ax.plot(self.time, self.emagoc_pol_axi, ls='--', c='#30a2da',
                    label='emag pol axi')
            ax.plot(self.time, self.emagoc_tor_axi, ls='--', c='#fc4f30',
                    label='emag tor axi')
            ax.plot(self.time, self.emag_tot, ls='-', c='#31363B')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Emag')
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
            ax.plot(self.time, self.emagic_pol, ls='-', c='#30a2da',
                    label='emagic pol')
            ax.plot(self.time, self.emagic_tor, ls='-', c='#fc4f30',
                    label='emagic tor')
            ax.plot(self.time, self.emagic_pol_axi, ls='--', c='#30a2da',
                    label='emagic pol axi')
            ax.plot(self.time, self.emagic_tor_axi, ls='--', c='#fc4f30',
                    label='emagic tor axi')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('emag inner core')
            fig.tight_layout()
        elif self.field == 'rot':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.omega_ic, ls='-', c='#30a2da',
                    label='Omega IC')
            ax.set_xlabel('Time')
            ax.set_ylabel('Rotation inner core')
            ax.legend(loc='best', frameon=False)
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
            ax1.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
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
                fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.dipTot, label='Total dipolarity')
            ax.plot(self.time, self.dipolarity, ls='--', c='#30a2da',
                    label='Axisym dipolarity')
            ax.plot(self.time, self.dipTot_cmb, ls='-', c='#6d904f',
                    label='Total dipolarity CMB')
            ax.plot(self.time, self.dip_cmb, ls='--', c='#6d904f',
                    label='Axisym dipolarity')
            ax.plot(self.time, self.dip_l11, ls='-', c='#fc4f30',
                    label='Axisym dip l=11')
            ax.plot(self.time, self.dipTot_l11, ls='--', c='#fc4f30',
                    label='Total dip l=11')
            # ax.plot(self.time, self.dip3, ls='-', c='#e5ae38',
            #         label='Epol axi/Ecmb')
            ax.legend(loc='best', frameon=False)
            ax.set_ylabel('Dipolarity')
            ax.set_xlabel('Time')
            ax.set_ylim(0,1)
            fig.tight_layout()
        elif self.field == 'AM':
            fig = plt.figure()
            ax = fig.add_subplot(211)
            ax.plot(self.time, self.am_oc_z, label='Outer core')
            ax.plot(self.time, self.am_ic, label='Inner core')
            ax.plot(self.time, self.amz, ls='-', c='#31363b', label='Total')
            ax.legend(loc='best', frameon=False)
            ax.set_ylabel('AM')
            ax = fig.add_subplot(212)
            ax.semilogy(self.time[1:], np.abs(self.damzdt[1:]))
            ax.set_xlabel('Time')
            ax.set_ylabel('dAmz / dt')
            fig.tight_layout()
        elif self.field == 'par':
            fig = plt.figure()
            ax = fig.add_subplot(111)
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
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()

            if self.dipolarity.max() > 0.:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(self.time, self.dipolarity, label='Dipolarity')
                ax.plot(self.time, self.dip_cmb, label='Dipolarity CMB')
                ax.legend(loc='upper right', frameon=False)
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
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.corr_vz_otc, label='uz')
            ax.plot(self.time, self.corr_vortz_otc, label='z vorticity')
            ax.plot(self.time, self.corr_hel_otc, label='Helicity')
            ax.set_xlabel('Time')
            ax.set_ylabel('z correlations')
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()

        elif self.field == 'earth_like':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.axial_dipole, label='AD/NAD')
            ax.plot(self.time, self.symmetry, label='O/E')
            ax.plot(self.time, self.zonality, label='Z/NZ')
            ax.plot(self.time, self.flux_concentration, label='FCF')
            ax.set_xlabel('Time')
            ax.set_ylabel('Rating parameters')
            ax.legend(loc='upper right', frameon=False)
            fig.tight_layout()

            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            ax1.plot(self.time, self.chi_square)
            ax1.set_xlabel('Time')
            ax1.set_ylabel('Chi square')
            fig1.tight_layout()
        elif self.field == 'misc':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.topnuss, label='Top Nusselt')
            ax.plot(self.time, self.botnuss, label='Bottom Nusselt')
            ax.legend(loc='lower right', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Nusselt number')
            fig.tight_layout()
            if self.helrms.max() != 0.:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(self.time, self.helrms)
                ax.set_xlabel('Time')
                ax.set_ylabel('Helicity')
                fig.tight_layout()
        elif self.field == 'heat':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            if self.kbots==2 and self.ktops==2:
                ax.plot(self.time, self.deltaTnuss, label=r'$Nu_{\Delta T}$')
            else:
                ax.plot(self.time, self.topnuss, label='Top Nusselt')
                ax.plot(self.time, self.botnuss, label='Bottom Nusselt')
            ax.legend(loc='lower right', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Nusselt number')
            fig.tight_layout()

            if self.topsherwood.max() != 1.0:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(self.time, self.topsherwood, label='Top Sherwood')
                ax.plot(self.time, self.botsherwood, label='Bottom Sherwood')
                ax.legend(loc='lower right', frameon=False)
                ax.set_xlabel('Time')
                ax.set_ylabel('Sherwood number')
                fig.tight_layout()
        elif self.field == 'helicity':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.helRMSN, label='Northern Hemisphere')
            ax.plot(self.time, self.helRMSS, label='Southern Hemisphere')
            ax.legend(loc='lower right')
            ax.set_xlabel('Time')
            ax.set_ylabel('Helicity')
            fig.tight_layout()
        elif self.field == 'u_square':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.ekin_pol, ls='-', c='#30a2da',
                    label='ekin pol')
            ax.plot(self.time, self.ekin_tor, ls='-', c='#fc4f30',
                    label='ekin tor')
            ax.plot(self.time, self.ekin_pol_axi, ls='--', c='#30a2da',
                    label='ekin pol axi')
            ax.plot(self.time, self.ekin_tor_axi, ls='--', c='#fc4f30',
                    label='ekin tor axi')
            ax.plot(self.time, self.ekin_tot, ls='-', c='#31363B')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('u**2')
            fig.tight_layout()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(self.time, self.rm, label='Magnetic Reynolds')
            ax.semilogy(self.time, self.ro, label='Rossby')
            ax.semilogy(self.time, self.rossby_l, label='Rossby l')
            ax.semilogy(self.time, self.dl, label='l')
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
            ax.semilogy(self.time, self.ChemRms, label='Chemical Buoyancy')
            ax.semilogy(self.time, self.InerRms, label='Inertia')
            ax.semilogy(self.time, self.DifRms, label='Diffusion')

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
            fig.tight_layout()
        elif self.field == 'perpPar':
            fig = plt.figure()
            ax= fig.add_subplot(111)
            ax.plot(self.time, self.eperp, ls='-', c='#30a2da',
                    label='ekin perp')
            ax.plot(self.time, self.epar, ls='-', c='#fc4f30',
                    label='ekin par')
            ax.plot(self.time, self.eperp_axi, ls='--', c='#30a2da',
                    label='ekin perp axi')
            ax.plot(self.time, self.epar_axi, ls='--', c='#fc4f30',
                    label='ekin par axi')
            ax.plot(self.time, self.ekin_tot, ls='-', c='#31363B')
            ax.plot(self.time, self.ekin_tot, 'k-')
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
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Power')
            fig.tight_layout()

            if hasattr(self,'fohm'):
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(self.time, self.fohm)
                ax.set_xlabel('Time')
                ax.set_ylabel('fohm')
                ax.set_ylim(0., 1.)
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

            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Torque')
            fig.tight_layout()
        elif self.field in ('am_mag_pol', 'am_mag_tor', 'am_kin_pol', 'am_kin_tor'):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            print(self.coeffs.shape)
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
            self.ekin_pol_symeq = data[:, 5]
            self.ekin_tor_symeq = data[:, 6]
            self.ekin_pol_asymeq = data[:, 7]
            self.ekin_tor_asymeq = data[:, 8]
            self.ekin_tot = self.ekin_pol + self.ekin_tor
            self.ekin_es = self.ekin_pol_symeq + self.ekin_tor_symeq
            self.ekin_eas = self.ekin_pol_asymeq + self.ekin_tor_asymeq
            self.ekin_pol_naxi=self.ekin_pol-self.ekin_pol_axi
            self.ekin_tor_naxi=self.ekin_tor-self.ekin_tor_axi
        elif self.field == 'e_mag_oc':
            self.time = data[:, 0]
            self.emagoc_pol = data[:, 1]
            self.emagoc_tor = data[:, 2]
            self.emagoc_pol_axi = data[:, 3]
            self.emagoc_tor_axi = data[:, 4]
            self.ext_nrj_pol = data[:, 5]
            self.ext_nrj_pol_axi = data[:, 6]
            self.emagoc_pol_es = data[:, 7]
            self.emagoc_tor_es = data[:, 8]
            self.emagoc_pol_eas = data[:, 9]
            self.emagoc_tor_eas = data[:, 10]
            self.emag_tot = self.emagoc_pol + self.emagoc_tor
            self.emag_es = self.emagoc_pol_es + self.emagoc_tor_es
            self.emag_eas = self.emagoc_pol_eas + self.emagoc_tor_eas
        elif self.field == 'e_mag_ic':
            self.time = data[:, 0]
            self.emagic_pol = data[:, 1]
            self.emagic_tor = data[:, 2]
            self.emagic_pol_axi = data[:, 3]
            self.emagic_tor_axi = data[:, 4]
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
            self.ratio = data[:, 17] # (e_cmb-e_as_cmb)/e_cmb
            self.epol_axi_cmb = (-self.ratio*self.ecmb+self.ecmb)
            self.dip3 = self.epol_axi_cmb/self.ecmb
            self.e_tot = self.e_dip/self.dipTot
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
            self.ekin = data[:, 6]
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
            else:
                self.ChemRms = data[:, 7]
                self.PreRms = data[:, 8]
                self.geos = data[:, 9] # geostrophic balance
                self.mageos = data[:, 10] # magnetostrophic balance
                self.arc    = data[:, 11] # Coriolis/Pressure/Buoyancy
                self.arcMag = data[:, 12] # Coriolis/Pressure/Buoyancy/Lorentz
                self.corLor = data[:, 13] # Coriolis/Lorentz
                self.preLor = data[:, 14] # Pressure/Lorentz
                self.cia = data[:, 15] # Coriolis/Inertia/Archmedean

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
                if timeOld != timeNew:
                    out.__dict__[attr] = np.hstack((self.__dict__[attr],
                                                    out.__dict__[attr]))
                else: # Same time
                    out.__dict__[attr] = np.hstack((self.__dict__[attr],
                                                    out.__dict__[attr][1:]))

        return out


class AvgField:
    """
    This class calculates the time-average properties from time series. It will
    store the input starting time in a small file named ``tInitAvg``, such that
    the next time you use it you don't need to give ``tstart`` again.

    >>> # Average from t=2.11 and also store the additional dipole.TAG informations
    >>> a = AvgField(tstart=2.11, dipExtra=True)
    >>> # Average only the files that match the pattern N0m2[a-c]
    >>> a = AvgField(tstart=2.11, tag='N0m2[a-c]')
    >>> # Average only the files that match the pattern N0m2Z*
    >>> a = AvgField(tstart=2.11, tag='N0m2Z*')
    >>> print(a) # print the formatted output
    """

    def __init__(self, tstart=None, tag=None, dipExtra=False, perpPar=False,
                 std=False):
        """
        :param tstart: the starting time for averaging
        :type tstart: float
        :param tag: if you specify an input tag (generic regExp pattern),
                    the averaging process will only happen on the time series
                    that match this input pattern
        :type tag: str
        :param dipExtra: if this parameter is set to ``True``, then additional
                         values extracted from :ref:`dipole.TAG <secDipoleFile>`
                         are also computed
        :type dipExtra: bool

        :param perpPar: additional values extracted from :ref:`perpPar.TAG <secperpParFile>`
                        are also computed
        :type perpPar: bool
        :type std: compute the standard deviation when set to True
        :type std: bool
        """

        if os.path.exists('tInitAvg') and tstart is None:
            file = open('tInitAvg', 'r')
            st = file.readline().strip('\n')
            tstart = float(st)
            file.close()
        elif tstart is not None:
            file = open('tInitAvg', 'w')
            file.write('{}'.format(tstart))
            file.close()
        self.std = std
        self.dipExtra = dipExtra
        self.perpPar = perpPar

        # e_kin file
        ts = MagicTs(field='e_kin', all=True, tag=tag, iplot=False)
        mask = np.where(abs(ts.time-tstart) == min(abs(ts.time-tstart)), 1, 0)
        ind = np.nonzero(mask)[0][0]

        self.integration_time = ts.time[ind:][-1]-ts.time[ind:][0]

        if self.std:
            self.ekin_pol_avg, self.ekin_pol_std = avgField(ts.time[ind:],
                                                   ts.ekin_pol[ind:], std=True)
            self.ekin_tor_avg, self.ekin_tor_std = avgField(ts.time[ind:],
                                                   ts.ekin_tor[ind:], std=True)
            self.ekin_pola_avg, self.ekin_pola_std = avgField(ts.time[ind:],
                                                     ts.ekin_pol_axi[ind:], std=True)
            self.ekin_tora_avg, self.ekin_tora_std = avgField(ts.time[ind:],
                                                     ts.ekin_tor_axi[ind:], std=True)
            self.ekin_tot_avg, self.ekin_tot_std = avgField(ts.time[ind:],
                                                            ts.ekin_pol[ind:]+ts.ekin_tor[ind:],std=True)
        else:
            self.ekin_pol_avg = avgField(ts.time[ind:], ts.ekin_pol[ind:])
            self.ekin_tor_avg = avgField(ts.time[ind:], ts.ekin_tor[ind:])
            self.ekin_pola_avg = avgField(ts.time[ind:], ts.ekin_pol_axi[ind:])
            self.ekin_tora_avg = avgField(ts.time[ind:], ts.ekin_tor_axi[ind:])
            self.ekin_tot_avg =  avgField(ts.time[ind:], ts.ekin_pol[ind:]+ts.ekin_tor[ind:])

        self.tavg = ts.time[-1]-ts.time[ind] # Averaging time

        self.ra = ts.ra
        self.prmag = ts.prmag
        self.pr = ts.pr
        self.ek = ts.ek
        if hasattr(ts, 'strat'):
            self.strat = ts.strat
        if hasattr(ts, 'DissNb'):
            self.strat = ts.DissNb
        self.mode = ts.mode

        # par file
        ts2 = MagicTs(field='par', all=True, iplot=False, tag=tag)
        mask = np.where(abs(ts2.time-tstart) == min(abs(ts2.time-tstart)), 1, 0)
        ind = np.nonzero(mask)[0][0]

        if self.std:
            self.dip, self.dip_std = avgField(ts2.time[ind:],
                                     ts2.dipolarity[ind:], std=True)
            self.dipCMB, self.dipCMB_std = avgField(ts2.time[ind:],
                                           ts2.dip_cmb[ind:], std=True)
            self.els, self.els_std = avgField(ts2.time[ind:],
                                     ts2.elsasser[ind:], std=True)
            self.elsCMB, self.elsCMB_std = avgField(ts2.time[ind:],
                                           ts2.els_cmb[ind:], std=True)
            self.rol, self.rol_std = avgField(ts2.time[ind:],
                                     ts2.rossby_l[ind:], std=True)
            self.reynolds, self.reynolds_std = avgField(ts2.time[ind:],
                                               ts2.rm[ind:], std=True)
            self.dlB, self.dlB_std = avgField(ts2.time[ind:],
                                     ts2.dlB[ind:], std=True)
            self.dmB, self.dmB_std = avgField(ts2.time[ind:],
                                     ts2.dmB[ind:], std=True)
            self.dlV, self.dlV_std = avgField(ts2.time[ind:],
                                     ts2.dlV[ind:], std=True)
            self.dmV, self.dmV_std = avgField(ts2.time[ind:],
                                     ts2.dmV[ind:], std=True)
            self.dlVc, self.dlVc_std = avgField(ts2.time[ind:],
                                       ts2.dlVc[ind:], std=True)
            self.lvDiss, self.lvDiss_std = avgField(ts2.time[ind:],
                                       ts2.lvDiss[ind:], std=True)
            self.lbDiss, self.lbDiss_std = avgField(ts2.time[ind:],
                                       ts2.lbDiss[ind:], std=True)
            self.elsassermod, self.elsassermod_std = avgField(ts2.time[ind:],
                                                              0.5*(ts2.elsasser/(ts2.rm*ts2.lbDiss))[ind:],
                                                              std=True)


        else:
            self.dip = avgField(ts2.time[ind:], ts2.dipolarity[ind:])
            self.dipCMB = avgField(ts2.time[ind:], ts2.dip_cmb[ind:])
            self.els = avgField(ts2.time[ind:], ts2.elsasser[ind:])
            self.elsCMB = avgField(ts2.time[ind:], ts2.els_cmb[ind:])
            self.rol = avgField(ts2.time[ind:], ts2.rossby_l[ind:])
            self.reynolds = avgField(ts2.time[ind:], ts2.rm[ind:])
            self.dlB = avgField(ts2.time[ind:], ts2.dlB[ind:])
            self.dmB = avgField(ts2.time[ind:], ts2.dmB[ind:])
            self.dlV = avgField(ts2.time[ind:], ts2.dlV[ind:])
            self.dmV = avgField(ts2.time[ind:], ts2.dmV[ind:])
            self.dlVc = avgField(ts2.time[ind:], ts2.dlVc[ind:])
            self.lvDiss = avgField(ts2.time[ind:], ts2.lvDiss[ind:])
            self.lbDiss = avgField(ts2.time[ind:], ts2.lbDiss[ind:])
            self.elsassermod = avgField(ts2.time[ind:],
                                        0.5*(ts2.elsasser/(ts2.rm*ts2.lbDiss))[ind:])

        # heat.TAG file
        if len(glob.glob('heat.*')) > 0:
            ts3 = MagicTs(field='heat', all=True, tag=tag, iplot=False)
        elif len(glob.glob('misc.*')) > 0:
            ts3 = MagicTs(field='misc', all=True, tag=tag, iplot=False)
        else:
            ts3 = None

        if ts3 != None:
            if (self.mode != 7 and self.mode != 8):
                mask = np.where(abs(ts3.time-tstart) == min(abs(ts3.time-tstart)), 1, 0)
                ind = np.nonzero(mask)[0][0]
                nuss = 0.5*(ts3.botnuss+ts3.topnuss)

                if self.std:
                    self.nuss, self.nuss_std = avgField(ts3.time[ind:], nuss[ind:], std=True)
                    try:
                        self.deltaTnuss, self.deltaTnuss_std = avgField(ts3.time[ind:], ts3.deltaTnuss[ind:], std=True)
                    except AttributeError:
                        pass
                else:
                    self.nuss = avgField(ts3.time[ind:], nuss[ind:])
                    try:
                        self.deltaTnuss = avgField(ts3.time[ind:], ts3.deltaTnuss[ind:], std=True)
                    except AttributeError:
                        pass
                    self.nubot = avgField(ts3.time[ind:],ts3.botnuss[ind:])
                    self.nutop = avgField(ts3.time[ind:],ts3.topnuss[ind:])

        if self.mode == 0 or self.mode == 8:
            # Emag OC file
            ts4 = MagicTs(field='e_mag_oc', all=True, iplot=False,
                          tag=tag)
            mask = np.where(abs(ts4.time-tstart) == min(abs(ts4.time-tstart)),
                           1, 0)
            ind = np.nonzero(mask)[0][0]
            emag_es = ts4.emagoc_pol_es+ts4.emagoc_tor_es

            if self.std:
                self.emag_pol_avg, self.emag_pol_std = avgField(ts4.time[ind:],
                                              ts4.emagoc_pol[ind:], std=True)
                self.emag_tor_avg, self.emag_tor_std = avgField(ts4.time[ind:],
                                              ts4.emagoc_tor[ind:], std=True)
                self.emag_pola_avg, self.emag_pola_std = avgField(ts4.time[ind:],
                                              ts4.emagoc_pol_axi[ind:], std=True)
                self.emag_tora_avg, self.emag_tora_std = avgField(ts4.time[ind:],
                                              ts4.emagoc_tor_axi[ind:], std=True)
                self.emag_es_avg, self.emag_es_std = avgField(ts4.time[ind:],
                                              emag_es[ind:], std=True)
                self.emag_tot_avg, self.emag_tot_std = avgField(ts4.time[ind:],
                                                                ts4.emagoc_pol[ind:]+ts4.emagoc_tor[ind:],std=True)
            else:
                self.emag_pol_avg = avgField(ts4.time[ind:], ts4.emagoc_pol[ind:])
                self.emag_tor_avg = avgField(ts4.time[ind:], ts4.emagoc_tor[ind:])
                self.emag_pola_avg = avgField(ts4.time[ind:],
                                              ts4.emagoc_pol_axi[ind:])
                self.emag_tora_avg = avgField(ts4.time[ind:],
                                              ts4.emagoc_tor_axi[ind:])
                self.emag_es_avg = avgField(ts4.time[ind:], emag_es[ind:])

                self.emag_tot_avg = avgField(ts4.time[ind:], ts4.emagoc_pol[ind:]+ts4.emagoc_tor[ind:])

            Emag_Ekin = (ts4.emagoc_pol[ind:]+ts4.emagoc_tor[ind:])/(ts.ekin_pol[ind:]+ts.ekin_tor[ind:])
            self.Emag_Ekin, self.Emag_Ekin_std = avgField(ts4.time[ind:], Emag_Ekin,std=True)

            if self.dipExtra:
                # dipole.TAG files
                ts5 = MagicTs(field='dipole', all=True, iplot=False, tag=tag)
                mask = np.where(abs(ts5.time-tstart) == min(abs(ts5.time-tstart)),
                               1, 0)
                ind = np.nonzero(mask)[0][0]

                if self.std:
                    self.dipTot, self.dipTot_std = avgField(ts5.time[ind:],
                                                   ts5.dipTot_cmb[ind:], std=True)
                    self.dipTotl11, self.dipTotl11_std = avgField(ts5.time[ind:],
                                                   ts5.dipTot_l11[ind:], std=True)
                    self.dipl11, self.dipl11_std = avgField(ts5.time[ind:],
                                                   ts5.dip_l11[ind:], std=True)
                    self.dip3, self.dip3_std = avgField(ts5.time[ind:],
                                                   ts5.dip3[ind:], std=True)
                    self.e_dip, self.e_dip_std = avgField(ts5.time[ind:],
                                                          ts5.e_dip[ind:], std=True)
                    self.e_dip_ax, self.e_dip_ax_std = avgField(ts5.time[ind:],
                                                                ts5.e_dip_ax[ind:],std=True)


                else:
                    self.dipTot = avgField(ts5.time[ind:], ts5.dipTot_cmb[ind:])
                    self.dipTotl11 = avgField(ts5.time[ind:],ts5.dipTot_l11[ind:])
                    self.dipl11 = avgField(ts5.time[ind:],ts5.dip_l11[ind:])
                    self.dip3 = avgField(ts5.time[ind:],ts5.dip3[ind:])
                    self.e_dip = avgField(ts5.time[ind:],ts5.e_dip[ind:])
                    self.e_dip_ax = avgField(ts5.time[ind:], ts5.e_dip_ax[ind:])

        # if len(glob.glob('dtVrms.*')) > 0:
        #     # dtVrms.TAG files
        #     tsrms = MagicTs(field='dtVrms', all=True, iplot=False,
        #                     tag=tag)
        #     mask = np.where(abs(tsrms.time-tstart) == min(abs(tsrms.time-tstart)),
        #                     1, 0)
        #     ind = np.nonzero(mask)[0][0]

        #     if self.std:
        #         self.dtVRms, self.dtVRms_std = avgField(tsrms.time[ind:], tsrms.dtVRms[ind:],std=True)
        #         self.CorRms, self.CorRms_std = avgField(tsrms.time[ind:], tsrms.CorRms[ind:],std=True)
        #         self.LFRms,  self.LFRms_std  = avgField(tsrms.time[ind:], tsrms.LFRms [ind:],std=True)
        #         self.AdvRms, self.AdvRms_std = avgField(tsrms.time[ind:], tsrms.AdvRms[ind:],std=True)
        #         self.DifRms, self.DifRms_std = avgField(tsrms.time[ind:], tsrms.DifRms[ind:],std=True)
        #         self.BuoRms, self.BuoRms_std = avgField(tsrms.time[ind:], tsrms.BuoRms[ind:],std=True)
        #         self.PreRms, self.PreRms_std = avgField(tsrms.time[ind:], tsrms.PreRms[ind:],std=True)
        #         self.geos,   self.geos_std   = avgField(tsrms.time[ind:], tsrms.geos[ind:]  ,std=True)
        #         self.mageos, self.mageos_std = avgField(tsrms.time[ind:], tsrms.mageos[ind:],std=True)
        #         self.arc,    self.arc_std    = avgField(tsrms.time[ind:], tsrms.arc[ind:]   ,std=True)
        #         self.arcMag, self.arcMag_std = avgField(tsrms.time[ind:], tsrms.arcMag[ind:],std=True)
        #         self.corLor, self.corLor_std = avgField(tsrms.time[ind:], tsrms.corLor[ind:],std=True)
        #         self.preLor, self.preLor_std = avgField(tsrms.time[ind:], tsrms.preLor[ind:],std=True)
        #         self.cia,    self.cia_std    = avgField(tsrms.time[ind:], tsrms.cia[ind:]   ,std=True)
        #         self.Elsasser_rms, self.Elsasser_rms_std = avgField(tsrms.time[ind:],
        #                                                             tsrms.LFRms[ind:]/tsrms.CorRms[ind:],std=True)
        #     else:
        #         self.dtVRms = avgField(tsrms.time[ind:], tsrms.dtVRms[ind:])
        #         self.CorRms = avgField(tsrms.time[ind:], tsrms.CorRms[ind:])
        #         self.LFRms  = avgField(tsrms.time[ind:], tsrms.LFRms [ind:])
        #         self.AdvRms = avgField(tsrms.time[ind:], tsrms.AdvRms[ind:])
        #         self.DifRms = avgField(tsrms.time[ind:], tsrms.DifRms[ind:])
        #         self.BuoRms = avgField(tsrms.time[ind:], tsrms.BuoRms[ind:])
        #         self.PreRms = avgField(tsrms.time[ind:], tsrms.PreRms[ind:])
        #         self.geos   = avgField(tsrms.time[ind:], tsrms.geos  [ind:])
        #         self.mageos = avgField(tsrms.time[ind:], tsrms.mageos[ind:])
        #         self.arc    = avgField(tsrms.time[ind:], tsrms.arc   [ind:])
        #         self.arcMag = avgField(tsrms.time[ind:], tsrms.arcMag[ind:])
        #         self.corLor = avgField(tsrms.time[ind:], tsrms.corLor[ind:])
        #         self.preLor = avgField(tsrms.time[ind:], tsrms.preLor[ind:])
        #         self.cia    = avgField(tsrms.time[ind:], tsrms.cia   [ind:])
        #         self.Elsasser_rms = avgField(tsrms.time[ind:],
        #                                      tsrms.LFRms[ind:]/tsrms.CorRms[ind:])

        if len(glob.glob('power.*')) > 0:
            # power.TAG files
            tspow = MagicTs(field='power', all=True, iplot=False,
                            tag=tag)
            mask = np.where(abs(tspow.time-tstart) == min(abs(tspow.time-tstart)),
                           1, 0)
            ind = np.nonzero(mask)[0][0]

            if self.std:
                self.viscDiss, self.viscDiss_std = avgField(tspow.time[ind:],
                                            -tspow.viscDiss[ind:], std=True)
                self.buoPower, self.buoPower_std = avgField(tspow.time[ind:],
                                             tspow.buoPower[ind:], std=True)
                if self.mode == 0 or self.mode == 8:
                    self.ohmDiss, self.ohmDiss_std = avgField(tspow.time[ind:],
                                             -tspow.ohmDiss[ind:], std=True)
                    self.fohm, self.fohm_std = avgField(tspow.time[ind:],
                                             tspow.fohm[ind:], std=True)
            else:
                self.viscDiss = avgField(tspow.time[ind:], -tspow.viscDiss[ind:])
                self.buoPower = avgField(tspow.time[ind:], tspow.buoPower[ind:])
                if self.mode == 0 or self.mode == 8:
                    self.ohmDiss = avgField(tspow.time[ind:], -tspow.ohmDiss[ind:])
                    if self.mode == 0:
                        self.fohm = avgField(tspow.time[ind:], tspow.fohm[ind:])

        else:
            self.ohmDiss = -1.
            self.viscDiss = -1.
            self.buoPower = 1.
            self.fohm = 1.

            if self.std:
                self.ohmDiss_std = 0.
                self.viscDiss_std = 0.
                self.buoPower_std = 0.
                self.fohm_std = 0.


        if len(glob.glob('u_square.*')) > 0 and self.strat > 0:
            # u_square.TAG files
            ts = MagicTs(field='u_square', all=True, iplot=False)
            mask = np.where(abs(ts.time-tstart) == min(abs(ts.time-tstart)), 1, 0)
            ind = np.nonzero(mask)[0][0]

            if self.std:
                self.ureynolds, self.ureynolds_std = avgField(ts.time[ind:],
                                                 ts.rm[ind:], std=True)
                self.urossby, self.urossby_std = avgField(ts.time[ind:],
                                                           ts.ro[ind:], std=True)
                self.o_urossby, self.o_urossby_std = avgField(ts.time[ind:],
                                                              1./ts.ro[ind:], std=True)
                self.urol, self.urol_std = avgField(ts.time[ind:],
                                           ts.rossby_l[ind:], std=True)
                self.udlV, self.udlV_std = avgField(ts.time[ind:],
                                           ts.dl[ind:], std=True)
                self.u2_pol, self.u2_pol_std = avgField(ts.time[ind:],
                                           ts.ekin_pol[ind:], std=True)
                self.u2_tor, self.u2_tor_std = avgField(ts.time[ind:],
                                           ts.ekin_tor[ind:], std=True)
                self.u2_pola, self.u2_pola_std = avgField(ts.time[ind:],
                                           ts.ekin_pol_axi[ind:], std=True)
                self.u2_tora, self.u2_tora_std = avgField(ts.time[ind:],
                                           ts.ekin_tor_axi[ind:], std=True)
            else:
                self.ureynolds = avgField(ts.time[ind:], ts.rm[ind:])
                self.urossby = avgField(ts.time[ind:], ts.ro[ind:])
                self.o_urossby =  1./self.urossby
                self.urol = avgField(ts.time[ind:], ts.rossby_l[ind:])
                self.udlV = avgField(ts.time[ind:], ts.dl[ind:])
                self.u2_pol = avgField(ts.time[ind:], ts.ekin_pol[ind:])
                self.u2_tor = avgField(ts.time[ind:], ts.ekin_tor[ind:])
                self.u2_pola = avgField(ts.time[ind:], ts.ekin_pol_axi[ind:])
                self.u2_tora = avgField(ts.time[ind:], ts.ekin_tor_axi[ind:])

        else:
            self.ureynolds = self.reynolds
            self.urol = self.rol
            self.udlV = self.dlV
            self.u2_pol = self.ekin_pol_avg
            self.u2_tor = self.ekin_tor_avg
            self.u2_pola = self.ekin_pola_avg
            self.u2_tora = self.ekin_tora_avg

            if self.std:
                self.ureynolds_std = self.reynolds_std
                self.urol_std = self.rol_std
                self.udlV_std = self.dlV_std
                self.u2_pol_std = self.ekin_pol_std
                self.u2_tor_std = self.ekin_tor_std
                self.u2_pola_std = self.ekin_pola_std
                self.u2_tora_std = self.ekin_tora_std

        if self.perpPar:
            if len(glob.glob('perpPar.*')) > 0:
                # perpPar.TAG files
                ts = MagicTs(field='perpPar', all=True, iplot=False)
                mask = np.where(abs(ts.time-tstart) == min(abs(ts.time-tstart)),
                       1, 0)
                ind = np.nonzero(mask)[0][0]

                if self.std:
                    self.eperp, self.eperp_std = avgField(ts.time[ind:],
                                                 ts.eperp[ind:], std=True)
                    self.epar, self.epar_std = avgField(ts.time[ind:],
                                               ts.epar[ind:], std=True)
                    self.eperp_axi, self.eperp_axi_std = avgField(ts.time[ind:],
                                               ts.eperp_axi[ind:], std=True)
                    self.epar_axi, self.epar_axi_std = avgField(ts.time[ind:],
                                               ts.epar_axi[ind:], std=True)
                else:
                    self.eperp = avgField(ts.time[ind:], ts.eperp[ind:])
                    self.epar = avgField(ts.time[ind:], ts.epar[ind:])
                    self.eperp_axi = avgField(ts.time[ind:], ts.eperp_axi[ind:])
                    self.epar_axi = avgField(ts.time[ind:], ts.epar_axi[ind:])

            else:
                self.eperp = 0.
                self.epar = 0.
                self.eperp_axi = 0.
                self.epar_axi = 0.

                if self.std:
                    self.eperp_std = 0.
                    self.epar_std = 0.
                    self.eperp_axi_std = 0.
                    self.epar_axi_std = 0.

    def __str__(self):
        """
        Formatted output
        """
        st_std=''
        if self.ek == -1:
            ek = 0. # to avoid the -1 for the non-rotating cases
        else:
            ek = self.ek
        if self.mode == 0 or self.mode == 8:
            st = '{:.3e}{:9.2e}{:9.2e}{:9.2e}{:5.2f}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}'.format(\
               self.ra, ek, self.pr, self.prmag, self.strat, self.ekin_pol_avg, \
               self.ekin_tor_avg, self.ekin_pola_avg, self.ekin_tora_avg, \
               self.u2_pol, self.u2_tor, self.u2_pola, self.u2_tora, \
               self.emag_pol_avg, self.emag_tor_avg,  self.emag_pola_avg, \
               self.emag_tora_avg, self.emag_es_avg)

            if self.std:
                st_std = '{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}{:12.5e}'.format(\
                   self.ekin_pol_std, \
                   self.ekin_tor_std, self.ekin_pola_std, self.ekin_tora_std, \
                   self.u2_pol_std, self.u2_tor_std, self.u2_pola_std, \
                   self.u2_tora_std, self.emag_pol_std, self.emag_tor_std, \
                   self.emag_pola_std, self.emag_tora_std, self.emag_es_std)

            st +='{:8.2f}{:8.2f}{:9.2e}{:9.2e}{:9.2e}{:9.2e}{:9.2e}{:9.2e}{:7.3f}{:9.2e}{:9.2e}{:9.2e}{:9.2e}{:9.2e}{:9.2e}'.format(\
                 self.reynolds, self.ureynolds, self.rol, self.urol, \
                 self.dip, self.dipCMB, self.els, self.elsCMB, self.nuss, \
                 self.dlV, self.dmV, self.udlV, self.dlVc, self.dlB, self.dmB)

            if self.std:
                st_std +='{:9.2e}{:9.2e}{:9.2e}{:9.2e}{:9.2e}{:9.2e}{:9.2e}{:9.2e}{:7.3f}{:9.2e}{:9.2e}{:9.2e}{:9.2e}{:9.2e}{:9.2e}'.format(\
                     self.reynolds_std, self.ureynolds_std, self.rol_std, \
                     self.urol_std, self.dip_std, self.dipCMB_std, self.els_std, \
                     self.elsCMB_std, self.nuss_std, self.dlV_std, self.dmV_std, \
                     self.udlV_std, self.dlVc_std, self.dlB_std, self.dmB_std)

            if self.dipExtra:
                st +='{:9.2e}{:9.2e}{:9.2e}{:9.2e}'.format(self.dipTot, self.dipl11, \
                                               self.dipTotl11, self.dip3)

                if self.std:
                    st_std +='{:9.2e}{:9.2e}{:9.2e}{:9.2e}'.format(self.dipTot_std, \
                             self.dipl11_std, self.dipTotl11_std, self.dip3_std)

            st += '{:12.5e}{:12.5e}{:12.5e}{:9.2e}'.format(self.buoPower, self.ohmDiss,\
                                               self.viscDiss, self.fohm)

            if self.std:
                st_std += '{:12.5e}{:12.5e}{:12.5e}{:9.2e}'.format(self.buoPower_std, \
                          self.ohmDiss_std, self.viscDiss_std, self.fohm_std)

        else:
            st = '{:.3e}{:12.5e}{:5.2f}{:6.2f}{:12.5e}{:12.5e}{:12.5e}{:12.5e}'.format(\
               self.ra, ek, self.strat, self.pr, self.ekin_pol_avg, \
               self.ekin_tor_avg, self.ekin_pola_avg, self.ekin_tora_avg)
            if self.std:
                st_std = '{:12.5e}{:12.5e}{:12.5e}{:12.5e}'.format(self.ekin_pol_std, \
                   self.ekin_tor_std, self.ekin_pola_std, self.ekin_tora_std)

            if self.strat == 0:
                self.u2_pol = self.ekin_pol_avg
                self.u2_tor = self.ekin_tor_avg
                self.u2_pola = self.ekin_pola_avg
                self.u2_tora = self.ekin_tora_avg
                self.urol = self.rol
                self.ureynolds = self.reynolds

                if self.std:
                    self.u2_pol_std = self.ekin_pol_std
                    self.u2_tor_std = self.ekin_tor_std
                    self.u2_pola = self.ekin_pola_std
                    self.u2_tora = self.ekin_tora_std
                    self.urol_std = self.rol_std
                    self.ureynolds_std = self.reynolds_std

            st += '{:12.5e}{:12.5e}{:12.5e}{:12.5e}'.format(\
                  self.u2_pol, self.u2_tor, self.u2_pola, self.u2_tora)
            if self.std:
                st_std += '{:12.5e}{:12.5e}{:12.5e}{:12.5e}'.format(\
                      self.u2_pol_std, self.u2_tor_std, self.u2_pola_std, \
                       self.u2_tora_std)
            st +='{:8.2f}{:8.2f}{:9.2e}{:9.2e}{:12.5e}{:9.2e}{:9.2e}{:9.2e}{:9.2e}'.format(\
               self.reynolds, self.ureynolds, self.rol, self.urol, \
               self.nuss, self.dlV, self.dmV, self.udlV, self.dlVc)
            if self.std:
                st_std +='{:9.2e}{:9.2e}{:9.2e}{:9.2e}{:12.5e}{:9.2e}{:9.2e}{:9.2e}{:9.2e}'.format(\
                   self.reynolds_std, self.ureynolds_std, self.rol_std, \
                   self.urol_std, self.nuss_std, self.dlV_std, self.dmV_std, \
                   self.udlV_std, self.dlVc_std)

            #st += '{:12.5e}{:12.5e}'.format(self.buoPower, self.viscDiss)

            #if self.std:
            #    st_std += '{:12.5e}{:12.5e}'.format(self.buoPower_std, self.viscDiss_std)

        if self.perpPar:
            st += '{:12.5e}{:12.5e}{:12.5e}{:12.5e}'.format(self.eperp, self.epar, \
                                                self.eperp_axi, self.epar_axi)

            if self.std:
                st_std += '{:12.5e}{:12.5e}{:12.5e}{:12.5e}'.format(self.eperp_std, \
                          self.epar_std, self.eperp_axi_std, self.epar_axi_std)
        st += st_std
        st += '\n'

        return st
