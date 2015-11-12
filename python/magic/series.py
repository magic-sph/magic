#-*- coding: utf-8 -*-
import os, re
import matplotlib.pyplot as P
import numpy as N
from .log import MagicSetup
import glob
from .libmagic import fast_read, scanDir
from scipy.integrate import trapz


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
       * Miscellaneous: :ref:`misc.TAG <secMiscFile>`
       * Velocity square: :ref:`u_square.TAG <secu_squareFile>`
       * Angular momentum: :ref:`AM.TAG <secAMFile>`
       * Power budget: :ref:`power.TAG <secpowerFile>`
       * Parallel and perpendicular decomposition: :ref:`perpPar.TAG <secperpParFile>`
       * RMS force balance: :ref:`dtVrms.TAG <secdtVrmsFile>`
       * RMS induction terms: :ref:`dtBrms.TAG <secdtBrmsFile>`

    Here are a couple of examples of how to use this function.

    >>> # plot the most recent e_kin.TAG file found in the directoy
    >>> MagicTs(field='e_kin')
    >>>
    >>> # stack **all** the power.TAG file found in the directory
    >>> ts = MagicTs(field='power', all=True)
    >>> print(ts.time, ts.buoPower) # print time and buoyancy power 
    >>>
    >>> # If you only want to read the file ``misc.N0m2z``
    >>> ts = MagicTs(field='misc', tag='N0m2z', iplot=False)
    """

    def __init__(self, datadir='.', field='e_kin', iplot=True, all=False, tag=None):
        """
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
        logFiles = scanDir('log.*')

        if tag is not None:
            files = scanDir('%s.%s' % (self.field, tag))

            # Either the log.tag directly exists and the setup is easy to obtain
            if os.path.exists('log.%s' % tag):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
            # Or the tag is a bit more complicated and we need to find 
            # the corresponding log file
            else:
                mask = re.compile(r'%s\.(.*)' % self.field)
                if mask.match(files[-1]):
                    ending = mask.search(files[-1]).groups(0)[0]
                    if logFiles.__contains__('log.%s' % ending):
                        MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml='log.%s' % ending)

            # Concatenate the files that correspond to the tag
            for k,file in enumerate(files):
                filename = os.path.join(datadir, file)
                if self.field in ('am_mag_pol','am_mag_tor','am_kin_pol','am_kin_tor'):
                    datanew = fast_read(filename, binary=True)
                else:
                    datanew = fast_read(filename)
                if k == 0:
                    data = datanew.copy()
                    ncolRef = data.shape[1]
                else:
                    ncol = datanew.shape[1]
                    if ncol == ncolRef:
                        data = N.vstack((data, datanew))
                    else: # If the number of columns has changed
                        data = N.vstack((data, datanew[:, 0:ncolRef]))

        # If no tag is specified, the most recent is plotted
        elif not all:
            if len(logFiles) != 0:
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml=logFiles[-1])
                name = '%s.%s' % (self.field, self.tag)
                filename = os.path.join(datadir, name)
                if self.field in ('am_mag_pol','am_mag_tor','am_kin_pol','am_kin_tor'):
                    data = fast_read(filename, binary=True)
                else:
                    data = fast_read(filename)
            else:
                mot = '%s.*' % (self.field)
                dat = [(os.stat(i).st_mtime, i) for i in glob.glob(mot)]
                dat.sort()
                filename = dat[-1][1]
                if self.field in ('am_mag_pol','am_mag_tor','am_kin_pol','am_kin_tor'):
                    data = fast_read(filename, binary=True)
                else:
                    data = fast_read(filename)

        # If no tag is specified but all=True, all the directory is plotted
        else:
            if len(logFiles) != 0:
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml=logFiles[-1])
            files = scanDir('%s.*' % (self.field))
            for k,file in enumerate(files):
                filename = os.path.join(datadir, file)
                if self.field in ('am_mag_pol','am_mag_tor','am_kin_pol','am_kin_tor'):
                    datanew = fast_read(filename, binary=True)
                else:
                    datanew = fast_read(filename)
                if k == 0:
                    data = datanew.copy()
                    ncolRef = data.shape[1]
                else:
                    ncol = datanew.shape[1]
                    if ncol == ncolRef:
                        data = N.vstack((data, datanew))
                    else: # If the number of columns has changed
                        data = N.vstack((data, datanew[:, 0:ncolRef]))

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
            self.emag_tot = self.emagoc_pol+self.emagoc_tor
        elif self.field == 'e_mag_ic':
            self.time = data[:, 0]
            self.emagic_pol = data[:, 1]
            self.emagic_tor = data[:, 2]
            self.emagic_pol_axi = data[:, 3]
            self.emagic_tor_axi = data[:, 4]
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
            self.ratio = data[:, 17] # (e_cmb-e_as_cmb)/e_cmb
            self.epol_axi_cmb = (-self.ratio*self.ecmb+self.ecmb)
            self.dip3 = self.epol_axi_cmb/self.ecmb
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
            self.ic_rot = data[:, 1]
            self.lo_ic_rot = data[:, 2]
            self.visc_ic_rot = data[:, 3]
            self.mantle_rot = data[:, 4]
            self.lo_mantle_rot = data[:, 5]
            self.visc_mantle_rot = data[:, 6]
            self.Bpol_re = data[:, 7]
            self.Bpol_im = data[:, 8]
            self.ang_mom = data[:, 9]
            self.ic_ang_mom = data[:, 10]
            self.oc_ang_mom = data[:, 11]
            self.mantle_ang_mom = data[:, 12]
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
            self.dlB = data[:, 13]
            self.dmB = data[:, 14]
            self.els_cmb = data[:, 15]
            try:
                self.rolc = data[:, 16]
                self.dlVc = data[:, 17]
                self.reEquat = data[:, 18]
            except IndexError:
                pass
        elif self.field == 'misc':
            self.time = data[:, 0]
            self.botnuss = data[:, 1]
            self.topnuss = data[:, 2]
            self.helrms = data[:, 8]
            self.helN = data[:, 5]*self.helrms
            self.helS = data[:, 6]*self.helrms
            try:
                self.botflux = data[:, 16]
                self.topflux = data[:, 17]
            except IndexError:
                pass
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
        elif self.field in ('dtVrms'):
            self.time = data[:, 0]
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
        elif self.field in ('dtBrms'):
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
        elif self.field in ('power'):
            self.time = data[:, 0]
            self.buoPower = data[:, 1]
            self.icrotPower = data[:, 2]
            self.mantelrotPower = data[:, 3]
            self.viscDiss = data[:, 4]
            self.ohmDiss = data[:, 5]
            self.icPower = data[:, 6]
            self.mantlePower = data[:, 7]
            self.fohm = -self.ohmDiss/self.buoPower
        elif self.field in ('am_mag_pol', 'am_mag_tor', # Tayler instability
                            'am_kin_pol', 'am_kin_tor'):
            self.time = data[:, 0]
            self.coeffs = data[:, 1:]

        if iplot:
            self.plot()

    def plot(self):
        """
        Plotting subroutines. Only called if 'iplot=True'
        """
        if self.field == 'e_kin':
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.ekin_pol, 'b-', label='ekin pol')
            ax.plot(self.time, self.ekin_tor, 'r-', label='ekin tor')
            ax.plot(self.time, self.ekin_pol_axi, 'b--', label='ekin pol axi')
            ax.plot(self.time, self.ekin_tor_axi, 'r--', label='ekin tor axi')
            ax.plot(self.time, self.ekin_tot, 'k-')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Ekin')
        elif self.field == 'e_mag_oc':
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.emagoc_pol, 'b-', label='emag pol')
            ax.plot(self.time, self.emagoc_tor, 'r-', label='emag tor')
            ax.plot(self.time, self.emagoc_pol_axi, 'b--', label='emag pol axi')
            ax.plot(self.time, self.emagoc_tor_axi, 'r--', label='emag tor axi')
            ax.plot(self.time, self.emag_tot, 'k-')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Emag')
        elif self.field == 'e_mag_ic':
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.emagic_pol, 'b-', label='emag_ic pol')
            ax.plot(self.time, self.emagic_tor, 'r-', label='emag_ic tor')
            ax.plot(self.time, self.emagic_pol_axi, 'b--',
                   label='emag_ic pol axi')
            ax.plot(self.time, self.emagic_tor_axi, 'r--',
                   label='emag_ic tor axi')
            ax.legend(loc='lower right')
            ax.set_xlabel('Time')
            ax.set_ylabel('emag inner core')
        elif self.field == 'dipole':
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.theta_dip, 'b-', label='theta_dip')
            #ax.plot(self.time, self.phi_dip, 'r-', label='phi_dip')
            ax.set_ylabel('Dipole angle')
            ax.set_xlabel('Time')
            ax.set_ylim(-1., 181)

            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.dipTot, 'b-', label='Total dipolarity')
            ax.plot(self.time, self.dipolarity, 'b--', label='Axisym dipolarity')
            ax.plot(self.time, self.dipTot_cmb, 'g-', label='Total dipoloarity CMB')
            ax.plot(self.time, self.dip_cmb, 'g--', label='Axisym dipolarity')
            ax.plot(self.time, self.dip_l11, 'k-', label='Axisym dip l=11')
            ax.plot(self.time, self.dipTot_l11, 'k--', label='Total dip l=11')
            ax.plot(self.time, self.dip3, 'r-', label='Epol axi/Ecmb')
            ax.legend(loc='best', frameon=False)
            ax.set_ylabel('Dipolarity')
            ax.set_xlabel('Time')
        elif self.field == 'AM':
            fig = P.figure()
            ax = fig.add_subplot(211)
            ax.plot(self.time, self.am_oc_z, label='Outer core')
            ax.plot(self.time, self.amz, 'k-', label='Total')
            ax.legend(loc='best')
            ax.set_ylabel('AM')
            ax = fig.add_subplot(212)
            ax.semilogy(self.time[1:], N.abs(self.damzdt[1:]))
            ax.set_xlabel('Time')
            ax.set_ylabel('dAmz / dt')
        elif self.field == 'par':
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(self.time, self.rm, label='Magnetic Reynolds')
            ax.semilogy(self.time, self.elsasser, label='Elsasser')
            ax.semilogy(self.time, self.els_cmb, label='Elsasser CMB')
            ax.semilogy(self.time, self.rossby_l, label='Rossby l')
            if hasattr(self, 'rolc'):
                ax.semilogy(self.time, self.rolc, label='Roc l')
            ax.legend(loc='lower right')
            ax.set_xlabel('Time')
            ax.set_ylabel('Params')
            if self.dipolarity.max() > 0.:
                fig = P.figure()
                ax = fig.add_subplot(111)
                ax.plot(self.time, self.dipolarity, label='Dipolarity')
                ax.plot(self.time, self.dip_cmb, label='Dipolarity CMB')
                ax.legend(loc='upper right')
                ax.set_xlabel('Time')
                ax.set_ylabel('Dipolarity')
        elif self.field == 'misc':
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.topnuss, label='Top Nusselt')
            ax.plot(self.time, self.botnuss, label='Bottom Nusselt')
            ax.legend(loc='lower right')
            ax.set_xlabel('Time')
            ax.set_ylabel('Nusselt number')
            if self.helrms.max() != 0.:
                fig = P.figure()
                ax = fig.add_subplot(111)
                ax.plot(self.time, self.helrms)
                ax.set_xlabel('Time')
                ax.set_ylabel('Helicity')

        elif self.field == 'u_square':
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.ekin_pol, 'b-', label='ekin pol')
            ax.plot(self.time, self.ekin_tor, 'r-', label='ekin tor')
            ax.plot(self.time, self.ekin_pol_axi, 'b--', label='ekin pol axi')
            ax.plot(self.time, self.ekin_tor_axi, 'r--', label='ekin tor axi')
            ax.plot(self.time, self.ekin_tot, 'k-')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('u**2')

            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(self.time, self.rm, label='Magnetic Reynolds')
            ax.semilogy(self.time, self.ro, label='Rossby')
            ax.semilogy(self.time, self.rossby_l, label='Rossby l')
            ax.semilogy(self.time, self.dl, label='l')
            ax.legend(loc='lower right')
            ax.set_xlabel('Time')
            ax.set_ylabel('Params')
        elif self.field in ('dtVrms'):
            fig = P.figure() # Poloidal forces
            ax = fig.add_subplot(111)
            ax.semilogy(self.time, self.dtVRms, label='Time derivative')
            ax.semilogy(self.time, self.CorRms, label='Coriolis')
            ax.semilogy(self.time, self.PreRms, label='Pressure')
            ax.semilogy(self.time, self.LFRms, label='Lorentz')
            ax.semilogy(self.time, self.BuoRms, label='Buoyancy')
            ax.semilogy(self.time, self.AdvRms, label='Inertia')
            ax.semilogy(self.time, self.DifRms, label='Diffusion')

            ax.legend(loc='best', frameon=False, ncol=2)
            ax.set_xlabel('Time')
            ax.set_ylabel('RMS forces')

            fig = P.figure() # Toroidal forces
            ax = fig.add_subplot(111)
            ax.semilogy(self.time, self.geos, label='Geostrophic balance')
            ax.semilogy(self.time, self.mageos, label='Magnetostrophic')
            ax.semilogy(self.time, self.arc, label='Archimedean')
            ax.semilogy(self.time, self.corLor, label='Coriolis/Lorentz')
            ax.semilogy(self.time, self.preLor, label='Pressure/Lorentz')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('RMS balances')

        elif self.field == 'perpPar':
            fig = P.figure()
            ax= fig.add_subplot(111)

            ax.plot(self.time, self.eperp, 'b-', label='Ekin perp.')
            ax.plot(self.time, self.epar, 'r-', label='Ekin par.')
            ax.plot(self.time, self.eperp_axi, 'b--', label='Ekin perp. axi.')
            ax.plot(self.time, self.epar_axi, 'r--', label='Ekin par. axi.')
            ax.plot(self.time, self.ekin_tot, 'k-')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Ekin')

            ax.set_xlim(self.time[0], self.time[-1])

        elif self.field in ('power'):
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(self.time, self.buoPower, label='Buoyancy')
            ax.semilogy(self.time, -self.ohmDiss, label='Ohmic diss.')
            ax.semilogy(self.time, -self.viscDiss, label='Viscous diss.')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Power')

            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.fohm)
            ax.set_xlabel('Time')
            ax.set_ylabel('fohm')
        elif self.field in ('dtBrms'):
            fig = P.figure() # Poloidal
            ax = fig.add_subplot(111)
            ax.semilogy(self.time, self.dtBpolRms, 'k-', label='time derivative')
            ax.semilogy(self.time, self.DynPolRms, 'r-', label='Induction')
            ax.semilogy(self.time, self.DifPolRms, 'b-', label='Diffusion')

            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Poloidal field production')

            fig = P.figure() # Toroidal
            ax = fig.add_subplot(111)
            ax.semilogy(self.time, self.dtBtorRms, 'k-', label='time derivative', )
            ax.semilogy(self.time, self.DynTorRms, 'r-', label='Induction')
            ax.semilogy(self.time, self.DifTorRms, 'b-', label='Diffusion')
            ax.semilogy(self.time, self.omEffect*self.DynTorRms, 'r--',
                        label='Omega effect')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Toroidal field production')

        elif self.field in ('am_mag_pol', 'am_mag_tor', 'am_kin_pol', 'am_kin_tor'):
            fig = P.figure()
            ax = fig.add_subplot(111)
            print(self.coeffs.shape)
            for k in range(self.coeffs.shape[1]):
                ax.semilogy(self.time, self.coeffs[:, k], label='m=%i'%k)
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

    def __init__(self, tstart=None, tag=None, dipExtra=False):
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
        """

        if os.path.exists('tInitAvg') and tstart is None:
            file = open('tInitAvg', 'r')
            st = file.readline().strip('\n')
            tstart = float(st)
            file.close()
        elif tstart is not None:
            file = open('tInitAvg', 'w')
            file.write('%f' % tstart)
            file.close()
        self.dipExtra = dipExtra
        ts = MagicTs(field='e_kin', all=True, tag=tag, iplot=False)
        mask = N.where(abs(ts.time-tstart) == min(abs(ts.time-tstart)), 1, 0)
        ind = N.nonzero(mask)[0][0]
        fac = 1./(ts.time.max()-ts.time[ind])
        self.ekin_pol_avg = fac * trapz(ts.ekin_pol[ind:], ts.time[ind:])
        self.ekin_tor_avg = fac * trapz(ts.ekin_tor[ind:], ts.time[ind:])
        self.ekin_pola_avg = fac * trapz(ts.ekin_pol_axi[ind:], ts.time[ind:])
        self.ekin_tora_avg = fac * trapz(ts.ekin_tor_axi[ind:], ts.time[ind:])

        self.ra = ts.ra
        self.prmag = ts.prmag
        self.pr = ts.pr
        self.ek = ts.ek
        self.strat = ts.strat
        self.mode = ts.mode

        ts2 = MagicTs(field='par', all=True, iplot=False, tag=tag)
        mask = N.where(abs(ts2.time-tstart) == min(abs(ts2.time-tstart)), 1, 0)
        ind = N.nonzero(mask)[0][0]
        fac = 1./(ts2.time.max()-ts2.time[ind])
        self.dip = fac * trapz(ts2.dipolarity[ind:], ts2.time[ind:])
        self.dipCMB = fac * trapz(ts2.dip_cmb[ind:], ts2.time[ind:])
        self.els = fac * trapz(ts2.elsasser[ind:], ts2.time[ind:])
        self.elsCMB = fac * trapz(ts2.els_cmb[ind:], ts2.time[ind:])
        self.rol = fac * trapz(ts2.rossby_l[ind:], ts2.time[ind:])
        self.reynolds = fac * trapz(ts2.rm[ind:], ts2.time[ind:])
        self.dlB = fac * trapz(ts2.dlB[ind:], ts2.time[ind:])
        self.dmB = fac * trapz(ts2.dmB[ind:], ts2.time[ind:])
        self.dlV = fac * trapz(ts2.dlV[ind:], ts2.time[ind:])
        self.dmV = fac * trapz(ts2.dmV[ind:], ts2.time[ind:])

        ts3 = MagicTs(field='misc', all=True, tag=tag, iplot=False)
        mask = N.where(abs(ts3.time-tstart) == min(abs(ts3.time-tstart)), 1, 0)
        ind = N.nonzero(mask)[0][0]
        fac = 1./(ts3.time.max()-ts3.time[ind])
        nussb = fac * trapz(ts3.botnuss[ind:], ts3.time[ind:])
        nusst = fac * trapz(ts3.topnuss[ind:], ts3.time[ind:])
        self.nuss = 0.5*(nussb+nusst)

        if self.mode == 0 or self.mode == 8:
            ts4 = MagicTs(field='e_mag_oc', all=True, iplot=False, 
                          tag=tag)
            mask = N.where(abs(ts4.time-tstart) == min(abs(ts4.time-tstart)), 
                           1, 0)
            ind = N.nonzero(mask)[0][0]
            fac = 1./(ts4.time.max()-ts4.time[ind])
            self.emag_pol_avg = fac * trapz(ts4.emagoc_pol[ind:], ts4.time[ind:])
            self.emag_tor_avg = fac * trapz(ts4.emagoc_tor[ind:], ts4.time[ind:])
            self.emag_pola_avg = fac* trapz(ts4.emagoc_pol_axi[ind:], 
                                            ts4.time[ind:])
            self.emag_tora_avg = fac* trapz(ts4.emagoc_tor_axi[ind:], 
                                            ts4.time[ind:])
            if self.dipExtra:
                ts5 = MagicTs(field='dipole', all=True, iplot=False, 
                              tag=tag)
                fac = 1./(ts5.time.max()-ts5.time[ind])
                self.dipTot = fac*trapz(ts5.dipTot_cmb[ind:], ts5.time[ind:])
                self.dipTotl11 = fac*trapz(ts5.dipTot_l11[ind:], ts5.time[ind:])
                self.dipl11 = fac*trapz(ts5.dip_l11[ind:], ts5.time[ind:])
                self.dip3 = fac*trapz(ts5.dip3[ind:], ts5.time[ind:])

            if len(glob.glob('power.*')) > 0:
                tspow = MagicTs(field='power', all=True, iplot=False,
                                tag=tag)
                mask = N.where(abs(tspow.time-tstart) == min(abs(tspow.time-tstart)), 
                               1, 0)
                ind = N.nonzero(mask)[0][0]
                fac = 1./(tspow.time.max()-tspow.time[ind])
                self.ohmDiss = fac*trapz(tspow.ohmDiss[ind:], tspow.time[ind:])
                self.buoPower = fac*trapz(tspow.buoPower[ind:], tspow.time[ind:])
                self.fohm = fac*trapz(tspow.fohm[ind:], tspow.time[ind:])
            else:
                self.ohmDiss = -1.
                self.buoPower = 1.
                self.fohm = 1.

        if len(glob.glob('u_square.*')) > 0 and self.strat > 0:
            ts = MagicTs(field='u_square', all=True, iplot=False)
            mask = N.where(abs(ts.time-tstart) == min(abs(ts.time-tstart)), 1, 0)
            ind = N.nonzero(mask)[0][0]
            fac = 1./(ts.time.max()-ts.time[ind])
            self.ureynolds = fac * trapz(ts.rm[ind:], ts.time[ind:])
            self.urol = fac * trapz(ts.rossby_l[ind:], ts.time[ind:])
            self.udlV = fac * trapz(ts.dl[ind:], ts.time[ind:])
            self.u2_pol = fac * trapz(ts.ekin_pol[ind:], ts.time[ind:])
            self.u2_tor = fac * trapz(ts.ekin_tor[ind:], ts.time[ind:])
            self.u2_pola = fac * trapz(ts.ekin_pol_axi[ind:], ts.time[ind:])
            self.u2_tora = fac * trapz(ts.ekin_tor_axi[ind:], ts.time[ind:])

        else:
            self.ureynolds = self.reynolds
            self.urol = self.rol
            self.udlV = self.dlV
            self.u2_pol = self.ekin_pol_avg
            self.u2_tor = self.ekin_tor_avg
            self.u2_pola = self.ekin_pola_avg
            self.u2_tora = self.ekin_tora_avg

    def __str__(self):
        """
        Formatted output
        """
        if self.ek == -1:
            ek = 0. # to avoid the -1 for the non-rotating cases
        else:
            ek = self.ek
        if self.mode == 0 or self.mode == 8:
            st = '%.3e%9.2e%9.2e%9.2e%5.2f%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e' % \
              (self.ra, ek, self.pr, self.prmag, self.strat, self.ekin_pol_avg, \
               self.ekin_tor_avg, self.ekin_pola_avg, self.ekin_tora_avg, \
               self.u2_pol, self.u2_tor, self.u2_pola, self.u2_tora, \
               self.emag_pol_avg, self.emag_tor_avg,  self.emag_pola_avg, \
               self.emag_tora_avg)
             
            st +='%8.2f%8.2f%9.2e%9.2e%9.2e%9.2e%9.2e%9.2e%7.3f%9.2e%9.2e%9.2e%9.2e' % \
                (self.reynolds, self.ureynolds, self.rol, self.urol, \
                 self.dip, self.dipCMB, self.els, self.elsCMB, self.nuss, \
                 self.dlV, self.udlV, self.dlB, self.dmB)
            if self.dipExtra:
                st +='%9.2e%9.2e%9.2e%9.2e' % (self.dipTot, self.dipl11, \
                                               self.dipTotl11, self.dip3)

            st += '%12.5e%12.5e%9.2e\n' % (self.buoPower, -self.ohmDiss, self.fohm)
        else:
            st = '%.3e%12.5e%5.2f%6.2f%12.5e%12.5e%12.5e%12.5e' % \
              (self.ra, ek, self.strat, self.pr, self.ekin_pol_avg, \
               self.ekin_tor_avg, self.ekin_pola_avg, self.ekin_tora_avg)
            if self.strat == 0:
                self.u2_pol = self.ekin_pol_avg
                self.u2_tor = self.ekin_tor_avg
                self.u2_pola = self.ekin_pola_avg
                self.u2_tora = self.ekin_tora_avg
                self.urol = self.rol
                self.ureynolds = self.reynolds
            st += '%12.5e%12.5e%12.5e%12.5e' % \
                  (self.u2_pol, self.u2_tor, self.u2_pola, self.u2_tora)
            st +='%8.2f%8.2f%9.2e%9.2e%12.5e%9.2e%9.2e%9.2e\n' % \
              (self.reynolds, self.ureynolds, self.rol, self.urol, \
               self.nuss, self.dlV, self.dmV, self.udlV)
        return st
