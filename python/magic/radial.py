# -*- coding: utf-8 -*-
import os, re
import matplotlib.pyplot as plt
import numpy as np
from .log import MagicSetup
from .libmagic import fast_read,scanDir


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

    def __init__(self, datadir='.', field='eKin', iplot=True, tag=None, tags=None):
        """
        :param field: the field you want to plot
        :type field: str
        :param iplot: to plot the output, default is True
        :type iplot: bool
        :param tag: a specific tag, default is None
        :type tag: str
        :param tags: a list that contains multiple tags: useful to sum
                     several radial files
        :type tags: list
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

        if tags is None:
            if tag is not None:
                file = '%s.%s' % (self.name, tag)
                filename = os.path.join(datadir, file)
                if os.path.exists('log.%s' % tag):
                    MagicSetup.__init__(self, datadir=datadir, quiet=True, 
                                        nml='log.%s' % tag)
            else:
                pattern = os.path.join(datadir, '%s.*'% self.name)
                files = scanDir(pattern)
                filename = files[-1]
                # Determine the setup
                mask = re.compile(r'%s\.(.*)' % self.name)
                ending = mask.search(files[-1]).groups(0)[0]
                if os.path.exists('log.%s' % ending):
                    try:
                        MagicSetup.__init__(self, datadir=datadir, quiet=True, 
                                        nml='log.%s' % ending)
                    except AttributeError:
                        pass

            if not os.path.exists(filename):
                print('No such file')
                return

            if self.name == 'varCond' or self.name == 'varVisc' or self.name == 'varDiff' \
               or self.name == 'anel':
                data = fast_read(filename, skiplines=1)
            else:
                data = fast_read(filename, skiplines=0)
        else:
            self.nsteps = 0
            for k, tag in enumerate(tags):
                nml = MagicSetup(quiet=True, datadir=datadir, nml='log.%s' % tag)
                self.nsteps += int(nml.steps_gone)
                file = '%s.%s' % (self.name, tag)
                filename = os.path.join(datadir, file)
                if k == 0:
                    self.tstart = nml.start_time
                    if self.name == 'bLayersR':
                        data = fast_read(filename, skiplines=0)
                        for i in [0, 1, 3, 4]:
                            data[:, i] *= (nml.stop_time-nml.start_time)
                        data[:, 2] *= (int(nml.steps_gone)-1)
                        if data.shape[1] == 6: # Patch if the last column is not there
                            self.gradT2 = data[:, 5]*(nml.stop_time-nml.start_time)
                            gradT2init = nml.start_time
                            gradT2finish = nml.stop_time
                    else:
                        data = fast_read(filename, skiplines=0)*(nml.stop_time-nml.start_time)
                else:
                    if self.name == 'bLayersR':
                        if os.path.exists(filename):
                            dat = fast_read(filename, skiplines=0)
                            for i in [0, 1, 3, 4]:
                                data[:, i] += dat[:, i]*(nml.stop_time-nml.start_time)
                            data[:, 2] = data[:, 2] + dat[:, 2]*(int(nml.steps_gone)-1)
                            if dat.shape[1] == 6: # Patch if the last column is not there
                                if hasattr(self, 'gradT2'):
                                    self.gradT2 = self.gradT2+dat[:, 5]*(nml.stop_time-nml.start_time)
                                    gradT2finish = nml.stop_time
                                else:
                                    self.gradT2 = dat[:, 5]*(nml.stop_time-nml.start_time)
                                    gradT2init = nml.start_time
                                    gradT2finish = nml.stop_time
                    elif self.name == 'eKinR':
                        if os.path.exists(filename):
                            dat = fast_read(filename, skiplines=0)
                            for i in [0, 1, 3, 4, 5, 6, 7, 8]:
                                data[:, i] += dat[:, i]*(nml.stop_time-nml.start_time)
                    else:
                        if os.path.exists(filename):
                            data += fast_read(filename, skiplines=0)*(nml.stop_time-nml.start_time)
            self.tstop = nml.stop_time
            if self.name == 'bLayersR':
                for i in [0, 1, 3, 4]:
                    data[:, i] /= (self.tstop-self.tstart)
                data[:, 2] /= (self.nsteps-1)

                if hasattr(self, 'gradT2'):
                    self.gradT2 = self.gradT2/(gradT2finish-gradT2init)
                    dat = np.zeros((data.shape[0], 6), dtype=self.gradT2.dtype)
                    dat[:, 0:5] = data[:, 0:5]
                    dat[:, 5] = self.gradT2
                    data = dat
            else:
                data /= (self.tstop-self.tstart)
            MagicSetup.__init__(self, datadir=datadir, quiet=True, nml='log.%s' % tags[-1])

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
            self.emag_tot = self.emag_pol+self.emag_tor
            self.ecmb = self.emag_tot[0]
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
            if data.shape[1] == 5:
                self.buoPower_chem = data[:, 2]
                self.viscDiss = data[:, 3]
                self.ohmDiss = data[:, 4]
            elif data.shape[1] == 4:
                self.viscDiss = data[:, 2]
                self.ohmDiss = data[:, 3]
        elif self.name == 'parrad':
            self.radius = data[:, 0]
            self.rm = data[:, 1]
            self.rol = data[:, 2]
            self.urol = data[:, 3]
            self.dlV = data[:, 4]
            self.udlV = data[:, 5]
            self.udlVc = data[:, 8]
            self.entropy = data[:, 9]
            self.varS = data[:, 10]
            self.uh = data[:, 11]
            self.duhdr = data[:, 12]
        elif self.name == 'parR':
            self.radius = data[:, 0]
            self.rm = data[:, 1]
            self.rol = data[:, 2]
            self.urol = data[:, 3]
            self.dlV = data[:, 4]
            self.dlVc = data[:, 5]
            self.udlV = data[:, 6]
            self.udlVc = data[:, 7]
        elif self.name == 'bLayersR':
            self.radius = data[:, 0]
            self.entropy = data[:, 1]
            self.varS = data[:, 2]
            self.uh = data[:, 3]
            self.duhdr = data[:, 4]
            try:
                self.dissS = data[:, 5]
            except IndexError:
                pass
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
        elif self.name == 'heatR':
            self.radius = data[:, 0]
            self.entropy = data[:, 1]
            self.temperature = data[:, 2]
            self.pressure = data[:, 3]
            self.density = data[:, 4]
        elif self.name == 'perpParR':
            self.radius = data[:, 0]
            self.Eperp = data[:, 1]
            self.Epar = data[:, 2]
            self.Eperp_axi = data[:, 3]
            self.Epar_axi = data[:, 4]

        if iplot:
            self.plot()

    def plot(self):
        """
        Display the result when ``iplot=True``
        """
        if self.name == 'eKinR':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.ekin_pol, ls='-', c='#30a2da',
                    label='ekin pol')
            ax.plot(self.radius, self.ekin_tor, ls='-', c='#fc4f30',
                    label='ekin tor')
            ax.plot(self.radius, self.ekin_pol_axi, ls='--', c='#30a2da',
                    label='ekin pol axi')
            ax.plot(self.radius, self.ekin_tor_axi, ls='--', c='#fc4f30',
                    label='ekin tor axi')
            ax.plot(self.radius, self.ekin_pol+self.ekin_tor, ls='-', c='#31363B',
                    label='ekin tot')
            ax.set_xlabel('Radius')
            ax.set_ylabel('Kinetic energy')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
        elif self.name == 'eMagR':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.emag_pol, ls='-', c='#30a2da',
                    label='emag pol')
            ax.plot(self.radius, self.emag_tor, ls='-', c='#fc4f30',
                    label='emag tor')
            ax.plot(self.radius, self.emag_pol_axi, ls='--', c='#30a2da',
                    label='emag pol axi')
            ax.plot(self.radius, self.emag_tor_axi, ls='--', c='#fc4f30',
                    label='emag tor axi')
            ax.plot(self.radius, self.emag_pol+self.emag_tor, ls='-', c='#31363B',
                    label='emag tot')

            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Radius')
            ax.set_ylabel('Magnetic energy')
            ax.set_xlim(self.radius.min(), self.radius.max())
            if hasattr(self, 'con_radratio'):
                if self.nVarCond == 2:
                    ax.axvline(self.con_radratio*self.radius[0], color='k',
                              linestyle='--')

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.dip_ratio)
            ax.set_xlabel('Radius')
            ax.set_ylabel('Dipolarity')
            ax.set_xlim(self.radius.min(), self.radius.max())
        elif self.name == 'anel':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(self.radius, self.temp0, label='Temperature')
            ax.semilogy(self.radius, self.rho0, label='Density')
            ax.set_ylabel('Reference state')
            ax.set_xlabel('Radius')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.beta, label='beta')
            ax.plot(self.radius, self.dbeta, label='dbeta/dr')
            ax.set_xlabel('Radius')
            ax.set_ylabel('Derivatives of rho')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.grav)
            ax.set_xlabel('Radius')
            ax.set_xlabel('Gravity')
            ax.set_xlim(self.radius.min(), self.radius.max())
        elif self.name == 'varDiff':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(self.radius, self.conduc, label='conductivity')
            ax.semilogy(self.radius, self.kappa, label='diffusivity')
            ax.semilogy(self.radius, self.prandtl, label='Prandtl')
            ax.set_ylabel('Thermal properties')
            ax.set_xlabel('Radius')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
            
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.dLkappa, label='dLkappa')
            ax.set_xlabel('Radius')
            ax.set_ylabel('$d\ln\kappa / dr$')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
        elif self.name == 'varVisc':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(self.radius, self.dynVisc, label='dyn. visc')
            ax.semilogy(self.radius, self.kinVisc, label='kin. visc')
            ax.semilogy(self.radius, self.prandtl, label='Prandtl')
            ax.semilogy(self.radius, self.ekman, label='Ekman')
            ax.set_ylabel('Thermal properties')
            ax.set_xlabel('Radius')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
            
            fig = plt.figure()
            ax = fig.add_subplot(111)
            axplot(self.radius, self.dLvisc, label='dLvisc')
            ax.set_xlabel('Radius')
            ax.set_ylabel(r'$d\ln\nu / dr$')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
        elif self.name == 'varCond':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(self.radius, self.conduc, label='conductivity')
            ax.set_ylabel('Electrical conductivity')
            ax.set_xlabel('Radius')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
        elif self.name == 'powerR':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.buoPower, label='Power')
            ax.plot(self.radius, self.viscDiss, label='visc diss')
            if self.ohmDiss.max() != 0.:
                ax.plot(self.radius, self.ohmDiss, label='ohm diss')
            ax.set_xlabel('Radius')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
        elif self.name == 'parrad' or self.name == 'bLayersR':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.entropy, label='entropy')
            ax.twinx()
            ax.plot(self.radius, self.varS/self.varS.max(), 
                    label='entropy variance')
            ax.set_xlabel('Radius')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.uh, label='uh')
            ax.plot(self.radius, self.duhdr, label='duhdr')
            ax.set_xlabel('Radius')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
        elif self.name == 'parR':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.rm)
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.set_xlabel('Radius')
            ax.set_ylabel('Rm')

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.rol, label='Rol')
            ax.plot(self.radius, self.urol, label='u Rol')
            ax.set_xlabel('Radius')
            ax.set_ylabel('Rol')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
        elif self.name == 'fluxesR':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.fcond/self.fcond[0], label='Fcond')
            ax.plot(self.radius, self.fconv/self.fcond[0], label='Fconv')
            ax.plot(self.radius, self.fkin/self.fcond[0], label='Fkin')
            ax.plot(self.radius, self.fvisc/self.fcond[0], label='Fvisc')
            if self.prmag != 0:
                ax.plot(self.radius, self.fpoyn/self.fcond[0], label='Fpoyn')
                ax.plot(self.radius, self.fres/self.fcond[0], label='Fres')

            ax.set_xlabel('Radius')
            ax.set_ylabel('Fluxes')
            ax.plot(self.radius, self.ftot/self.fcond[0])
            ax.legend(loc='best', frameon=False)
            ax.set_xlim(self.radius[-1], self.radius[0])
        elif self.name == 'heatR':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.entropy, label='s')
            ax.plot(self.radius, self.temperature, label='T')
            ax.set_xlabel('Radius')
            ax.set_ylabel('Temperature, Entropy')
            ax.legend(loc='best', frameon=False)
            ax.set_xlim(self.radius[-1], self.radius[0])

            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            ax1.plot(self.radius, self.pressure, label='p')
            ax1.set_xlabel('Radius')
            ax1.set_ylabel('Pressure')
            ax1.set_xlim(self.radius[-1], self.radius[0])

        elif self.name == 'perpParR':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.Eperp, 'b-', label='E perp.')
            ax.plot(self.radius, self.Epar, 'r-', label='E par.')
            ax.plot(self.radius, self.Eperp_axi, 'b--', label='E eperp. ax.')
            ax.plot(self.radius, self.Epar_axi, 'r--', label='E par. ax.')
            ax.set_xlabel('Radius')
            ax.set_ylabel('Kinetic energy')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)

        if hasattr(self, 'con_radratio'):
            if self.nVarCond == 2:
                ax.axvline(self.con_radratio*self.radius[0], color='k',
                          linestyle='--')
