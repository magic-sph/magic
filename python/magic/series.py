# -*- coding: utf-8 -*-
import os, re
import pylab as P
import numpy as N
from setup import MagicSetup
import string
import glob
from libmagic import fast_read, scanDir
from scipy.integrate import trapz

__author__  = "$Author$"
__date__   = "$Date$"
__version__ = "$Revision$"


class MagicTs(MagicSetup):

    def __init__(self, datadir='.', field='e_kin', iplot=True, all=False,
                 tag=None):
        """
        A class to plot time series of the MagIC code

        :param field: the file you want to plot
        :param iplot: display/hide the plot
        :param all: a boolean if you want to get the complete time series
        :param tag: if you specify a tag, it tries to build the
                        corresponding time series
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
                data = fast_read(filename)
            else:
                mot = '%s.*' % (self.field)
                dat = [(os.stat(i).st_mtime, i) for i in glob.glob(mot)]
                dat.sort()
                filename = dat[-1][1]
                data = fast_read(filename)

        # If no tag is specified but all=True, all the directory is plotted
        else:
            if len(logFiles) != 0:
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml=logFiles[-1])
            files = scanDir('%s.*' % (self.field))
            for k,file in enumerate(files):
                filename = os.path.join(datadir, file)
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
        elif self.field in ('dtVrms', 'dtVAsRms'):
            self.time = data[:, 0]
            self.dtVPolRms = data[:, 1]
            self.dtVTorRms = data[:, 2]
            self.CorPolRms = data[:, 3]
            self.CorTorRms = data[:, 4]
            self.LFPolRms = data[:, 5]
            self.LFTorRms = data[:, 6]
            self.AdvPolRms = data[:, 7]
            self.AdvTorRms = data[:, 8]
            self.DifPolRms = data[:, 9]
            self.DifTorRms = data[:, 10]
            self.BuoRms = data[:, 11]
            self.PreRms = data[:, 12]
            self.geos = data[:, 13] # geostrophic balance
            self.mgeos = data[:, 14] # magnetostrophic balance
            self.archim = data[:, 15] # archimedean balance
        elif self.field in ('dtBrms'):
            self.time = data[:, 0]
            self.dtBpolRms = data[:, 1]
            self.dtBtorRms = data[:, 2]
            self.StrPolRms = data[:, 3]
            self.StrTorRms = data[:, 4]
            self.AdvPolRms = data[:, 5]
            self.AdvTorRms = data[:, 6]
            self.DifPolRms = data[:, 7]
            self.DifTorRms = data[:, 8]
            self.omEffect = data[:, 9]
            self.omega = data[:, 10]
            self.DynPolRms = data[:, 11]
            self.DynTorRms = data[:, 12]
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
            iplot = False

        if iplot:
            self.plot()

    def plot(self):
        """
        Plotting subroutines. Only called if 'iplot=True'
        """
        if self.field == 'e_kin':
            P.figure()
            P.plot(self.time, self.ekin_pol, 'b-', label='ekin pol')
            P.plot(self.time, self.ekin_tor, 'r-', label='ekin tor')
            P.plot(self.time, self.ekin_pol_axi, 'b--', label='ekin pol axi')
            P.plot(self.time, self.ekin_tor_axi, 'r--', label='ekin tor axi')
            P.plot(self.time, self.ekin_tot, 'k-')
            P.legend(loc='best', frameon=False)
            P.xlabel('Time', fontsize=18)
            P.ylabel('Ekin', fontsize=18)
        elif self.field == 'e_mag_oc':
            P.figure()
            P.plot(self.time, self.emagoc_pol, 'b-', label='emag pol')
            P.plot(self.time, self.emagoc_tor, 'r-', label='emag tor')
            P.plot(self.time, self.emagoc_pol_axi, 'b--', label='emag pol axi')
            P.plot(self.time, self.emagoc_tor_axi, 'r--', label='emag tor axi')
            P.plot(self.time, self.emag_tot, 'k-')
            P.legend(loc='best', frameon=False)
            P.xlabel('Time', fontsize=18)
            P.ylabel('Emag', fontsize=18)
        elif self.field == 'e_mag_ic':
            P.figure()
            P.plot(self.time, self.emagic_pol, 'b-', label='emag_ic pol')
            P.plot(self.time, self.emagic_tor, 'r-', label='emag_ic tor')
            P.plot(self.time, self.emagic_pol_axi, 'b--',
                   label='emag_ic pol axi')
            P.plot(self.time, self.emagic_tor_axi, 'r--',
                   label='emag_ic tor axi')
            P.legend(loc='lower right')
            P.xlabel('Time', fontsize=18)
            P.ylabel('emag inner core', fontsize=18)
        elif self.field == 'dipole':
            P.figure()
            P.xlabel('Time')
            P.plot(self.time, self.theta_dip, 'b-', label='theta_dip')
            #P.plot(self.time, self.phi_dip, 'r-', label='phi_dip')
            P.ylabel('Dipole angle')
            P.ylim(-1., 181)
            P.figure()
            P.plot(self.time, self.dipTot, 'b-', label='Total dipolarity')
            P.plot(self.time, self.dipolarity, 'b--', label='Axisym dipolarity')
            P.plot(self.time, self.dipTot_cmb, 'g-', label='Total dipoloarity CMB')
            P.plot(self.time, self.dip_cmb, 'g--', label='Axisym dipolarity')
            P.plot(self.time, self.dip_l11, 'k-', label='Axisym dip l=11')
            P.plot(self.time, self.dipTot_l11, 'k--', label='Total dip l=11')
            P.plot(self.time, self.dip3, 'r-', label='Epol axi/Ecmb')
            P.legend(loc='best', frameon=False)
            P.ylabel('Dipolarity', fontsize=18)
            P.xlabel('Time', fontsize=18)
            P.figure()
            P.plot(self.time, self.epol_axi_cmb, 'b-')
            P.xlabel('Time', fontsize=18)
        elif self.field == 'AM':
            P.figure()
            #P.plot(self.time, self.am_ic, label='Inner core')
            P.subplot(211)
            P.plot(self.time, self.am_oc_z, label='Outer core')
            P.plot(self.time, self.amz, 'k-', label='Total')
            P.legend(loc='best')
            P.ylabel('Angular momentum')
            P.subplot(212)
            P.semilogy(self.time[1:], N.abs(self.damzdt[1:]))
            P.xlabel('Time')
            P.ylabel('dAmz / dt')
        elif self.field == 'par':
            P.figure()
            P.semilogy(self.time, self.rm, label='Magnetic Reynolds')
            P.semilogy(self.time, self.elsasser, label='Elsasser')
            P.semilogy(self.time, self.els_cmb, label='Elsasser CMB')
            P.semilogy(self.time, self.rossby_l, label='Rossby l')
            P.legend(loc='lower right')
            P.xlabel('Time', fontsize=18)
            P.ylabel('Params', fontsize=18)
            if self.dipolarity.max() > 0.:
                P.figure()
                P.plot(self.time, self.dipolarity, label='Dipolarity')
                P.plot(self.time, self.dip_cmb, label='Dipolarity CMB')
                P.legend(loc='upper right')
                P.xlabel('Time', fontsize=18)
                P.ylabel('Dipolarity', fontsize=18)
        elif self.field == 'misc':
            P.figure()
            P.plot(self.time, self.topnuss, label='Top Nusselt')
            P.plot(self.time, self.botnuss, label='Bottom Nusselt')
            P.legend(loc='lower right')
            P.xlabel('Time', fontsize=18)
            P.ylabel('Nusselt number', fontsize=18)
            #if hasattr(self, 'botflux'):
                #P.figure()
                #P.plot(self.time, self.topflux, label='Top Flux')
                #P.plot(self.time, self.botflux, label='Bottom Flux')
                #P.legend(loc='lower right')
                #P.xlabel('Time', fontsize=18)
                #P.ylabel('Flux', fontsize=18)

            if self.helrms.max() != 0.:
                P.figure()
                P.plot(self.time, self.helrms)
                P.xlabel('Time', fontsize=18)
                P.ylabel('Helicity', fontsize=18)

        elif self.field == 'u_square':
            P.figure()
            P.plot(self.time, self.ekin_pol, 'b-', label='ekin pol')
            P.plot(self.time, self.ekin_tor, 'r-', label='ekin tor')
            P.plot(self.time, self.ekin_pol_axi, 'b--', label='ekin pol axi')
            P.plot(self.time, self.ekin_tor_axi, 'r--', label='ekin tor axi')
            P.plot(self.time, self.ekin_tot, 'k-')
            P.legend(loc='best', frameon=False)
            P.xlabel('Time', fontsize=18)
            P.ylabel('u**2', fontsize=18)

            P.figure()
            P.semilogy(self.time, self.rm, label='Magnetic Reynolds')
            P.semilogy(self.time, self.ro, label='Rossby')
            P.semilogy(self.time, self.rossby_l, label='Rossby l')
            P.semilogy(self.time, self.dl, label='l')
            P.legend(loc='lower right')
            P.xlabel('Time')
            P.ylabel('Params')
        elif self.field in ('dtVrms', 'dtVAsRms'):
            P.figure() # Poloidal forces
            P.subplots_adjust(top=0.95, right=0.95)
            P.semilogy(self.time, self.dtVPolRms, label='Time derivative')
            P.semilogy(self.time, self.CorPolRms, label='Coriolis')
            P.semilogy(self.time, self.PreRms, label='Pressure')
            P.semilogy(self.time, self.LFPolRms, label='Lorentz')
            P.semilogy(self.time, self.BuoRms, label='Buoyancy')
            P.semilogy(self.time, self.AdvPolRms, label='Inertia')
            P.semilogy(self.time, self.DifPolRms, label='Diffusion')

            P.legend(loc='best', frameon=False)
            P.xlabel('Time', fontsize=18)
            P.ylabel('Poloidal RMS forces', fontsize=18)

            P.figure() # Toroidal forces
            P.semilogy(self.time, self.dtVTorRms, label='Time derivative')
            P.semilogy(self.time, self.CorTorRms, label='Coriolis')
            P.semilogy(self.time, self.LFTorRms, label='Lorentz')
            P.semilogy(self.time, self.AdvTorRms, label='Inertia')
            P.semilogy(self.time, self.DifTorRms, label='Diffusion')
            P.legend(loc='best', frameon=False)
            P.xlabel('Time', fontsize=18)
            P.ylabel('Toroidal RMS forces', fontsize=18)

        elif self.field in ('power'):
            P.figure()
            P.semilogy(self.time, self.buoPower, label='Buoyancy')
            P.semilogy(self.time, -self.ohmDiss, label='Ohmic diss.')
            P.semilogy(self.time, -self.viscDiss, label='Viscous diss.')
            P.legend(loc='best', frameon=False)
            P.xlabel('Time', fontsize=18)
            P.ylabel('Power', fontsize=18)

            P.figure()
            P.plot(self.time, self.fohm)
            P.xlabel('Time', fontsize=18)
            P.ylabel('fohm', fontsize=18)
        elif self.field in ('dtBrms'):
            P.figure() # Poloidal forces
            P.subplots_adjust(top=0.95, right=0.95)
            P.semilogy(self.time, self.dtBpolRms, label='time derivative')
            P.semilogy(self.time, self.StrPolRms, label='Stretching')
            P.semilogy(self.time, self.AdvPolRms, label='Advection')
            P.semilogy(self.time, self.DifPolRms, label='Diffusion')
            P.semilogy(self.time, self.DynPolRms, label='Dynamo')

            P.legend(loc='best', frameon=False)
            P.xlabel('Time', fontsize=18)
            P.ylabel('Poloidal field production', fontsize=18)

            P.figure() # Toroidal forces
            P.semilogy(self.time, self.dtBtorRms, label='time derivative')
            P.semilogy(self.time, self.StrTorRms, label='Stretching')
            P.semilogy(self.time, self.AdvTorRms, label='Advection')
            P.semilogy(self.time, self.DifTorRms, label='Diffusion')
            P.semilogy(self.time, self.DynTorRms, label='Dynamo')
            P.legend(loc='best', frameon=False)
            P.xlabel('Time', fontsize=18)
            P.ylabel('Toroidal field production', fontsize=18)

        P.subplots_adjust(top=0.95, right=0.95)


class AvgField:

    def __init__(self, tstart, tag=None, dipExtra=False):
        """
        A class to get average properties from time series

        :param tstart: the starting time for averaging
        :param tag: if you specify a tag, it tries to build the
                        corresponding time series, starting from your pattern
        :param dipExtra: if this parameter is set to true, then additional
                         values extracted from dipole.tag are computed
        """

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
        #mask = N.where(abs(ts2.time-tstart) == min(abs(ts2.time-tstart)), 1, 0)
        #ind = N.nonzero(mask)[0][0]
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

        ts3 = MagicTs(field='misc', all=True, tag=tag, iplot=False)
        #mask = N.where(abs(ts3.time-tstart) == min(abs(ts3.time-tstart)), 1, 0)
        #ind = N.nonzero(mask)[0][0]
        fac = 1./(ts3.time.max()-ts3.time[ind])
        self.nuss = fac * trapz(ts3.botnuss[ind:], ts3.time[ind:])

        if self.mode == 0:
            ts4 = MagicTs(field='e_mag_oc', all=True, iplot=False, 
                          tag=tag)
            #mask = N.where(abs(ts4.time-tstart) == min(abs(ts4.time-tstart)), 
                           #1, 0)
            #ind = N.nonzero(mask)[0][0]
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
                mask = N.where(abs(tspow.time-tstart) == min(abs(tspow.time-tstart)), 1, 0)
                ind = N.nonzero(mask)[0][0]
                fac = 1./(tspow.time.max()-tspow.time[ind])
                self.ohmDiss = fac*trapz(tspow.ohmDiss[ind:], tspow.time[ind:])
                self.buoPower = fac*trapz(tspow.buoPower[ind:], tspow.time[ind:])
                self.fohm = fac*trapz(tspow.fohm[ind:], tspow.time[ind:])
            else:
                self.ohmDiss = -1.
                self.buoPower = 1.
                self.fohm = 1.

        if len(glob.glob('u_square.*')) > 0:
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
        if self.ek == -1:
            ek = 0. # to avoid the -1 for the non-rotating cases
        else:
            ek = self.ek
        if self.mode == 0:
            st = '%.3e%9.2e%9.2e%9.2e%5.2f%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e' % \
              (self.ra, ek, self.pr, self.prmag, self.strat, self.ekin_pol_avg, \
               self.ekin_tor_avg, self.ekin_pola_avg, self.ekin_tora_avg, \
               self.u2_pol, self.u2_tor, self.u2_pola, self.u2_tora, \
               self.emag_pol_avg, self.emag_tor_avg,  self.emag_pola_avg, \
               self.emag_tora_avg)
             
            if self.dipExtra:
                st +='%8.2f%8.2f%9.2e%9.2e%9.2e%9.2e%9.2e%9.2e%7.3f%9.2e%9.2e%9.2e%9.2e' % \
                    (self.reynolds, self.ureynolds, self.rol, self.urol, \
                     self.dip, self.dipCMB, self.els, self.elsCMB, self.nuss, \
                     self.dlV, self.udlV, self.dlB, self.dmB)
                st +='%9.2e%9.2e%9.2e%9.2e' % (self.dipTot, self.dipl11, \
                                               self.dipTotl11, self.dip3)
            else:
                st +='%8.2f%8.2f%9.2e%9.2e%9.2e%9.2e%9.2e%9.2e%7.3f%9.2e%9.2e%9.2e%9.2e' % \
                    (self.reynolds, self.ureynolds, self.rol, self.urol, \
                     self.dip, self.dipCMB, self.els, self.elsCMB, self.nuss, \
                     self.dlV, self.udlV, self.dlB, self.dmB)

            st += '%12.5e%12.5e%9.2e\n' % (self.buoPower, -self.ohmDiss, self.fohm)
        else:
            st = '%.3e%12.5e%5.2f%6.2f%12.5e%12.5e%12.5e%12.5e' % \
              (self.ra, ek, self.strat, self.pr, self.ekin_pol_avg, \
               self.ekin_tor_avg, self.ekin_pola_avg, self.ekin_tora_avg)
            st += '%12.5e%12.5e%12.5e%12.5e' % \
                  (self.u2_pol, self.u2_tor, self.u2_pola, self.u2_tora)
            st +='%8.2f%8.2f%9.2e%9.2e%9.2e%9.2e%9.2e\n' % \
              (self.reynolds, self.ureynolds, self.rol, self.urol, \
               self.nuss, self.dlV, self.udlV)
        return st
