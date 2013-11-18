# -*- coding: utf-8 -*-
import os, re
import pylab as P
import numpy as N
from log import MagicSetup
from libmagic import fast_read,scanDir


__author__  = "$Author$"
__date__   = "$Date$"
__version__ = "$Revision$"


class MagicRadial(MagicSetup):
    """
    Read the radial files of the Magic Code: eKinR.TAG, eMagR.TAG,
    anel.TAG, varDiff.TAG
    """

    def __init__(self, datadir='.', field='eKin', iplot=True, tag=None):
        """
        :param field: the field you want to plot
        :param iplot: to plot the output, default is True
        :param tag: a specific tag, default is None
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
        else:
            print 'No corresponding radial profiles... Try again'

        if tag is not None:
            file = '%s.%s' % (self.name, tag)
            filename = os.path.join(datadir, file)
            if os.path.exists('log.%s' % tag):
                MagicSetup.__init__(self, datadir=datadir, quiet=True, 
                                    nml='log.%s' % tag)
        else:
            files = scanDir('%s.*'% self.name)
            filename = os.path.join(datadir, files[-1])
            # Determine the setup
            mask = re.compile(r'%s\.(.*)' % self.name)
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists('log.%s' % ending):
                MagicSetup.__init__(self, datadir=datadir, quiet=True, 
                                    nml='log.%s' % ending)

        if not os.path.exists(filename):
            print 'No such file'
            return

        if self.name == 'varCond' or self.name == 'varVisc' or self.name == 'varDiff' \
           or self.name == 'anel':
            data = fast_read(filename, skiplines=1)
        else:
            data = fast_read(filename, skiplines=0)

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

        if iplot:
            self.plot()

    def plot(self):
        if self.name == 'eKinR':
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.ekin_pol+self.ekin_tor, 'k-', 
                   label='ekin tot')
            ax.plot(self.radius, self.ekin_pol, 'b-', label='ekin pol')
            ax.plot(self.radius, self.ekin_tor, 'r-', label='ekin tor')
            ax.plot(self.radius, self.ekin_pol_axi, 'b--', label='ekin pol axi')
            ax.plot(self.radius, self.ekin_tor_axi, 'r--', label='ekin tor axi')
            ax.set_xlabel('Radius')
            ax.set_ylabel('Kinetic energy')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
        elif self.name == 'eMagR':
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.emag_pol+self.emag_tor, 'k-', 
                   label='emag tot')
            ax.plot(self.radius, self.emag_pol, 'b-', label='emag pol')
            ax.plot(self.radius, self.emag_tor, 'r-', label='emag tor')
            ax.plot(self.radius, self.emag_pol_axi, 'b--', label='emag pol axi')
            ax.plot(self.radius, self.emag_tor_axi, 'r--', label='emag tor axi')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Radius')
            ax.set_ylabel('Magnetic energy')
            ax.set_xlim(self.radius.min(), self.radius.max())
            if hasattr(self, 'con_radratio'):
                if self.nVarCond == 2:
                    ax.axvline(self.con_radratio*self.radius[0], color='k',
                              linestyle='--')

            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.dip_ratio)
            ax.set_xlabel('Radius')
            ax.set_ylabel('Dipolarity')
            ax.set_xlim(self.radius.min(), self.radius.max())
        elif self.name == 'anel':
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(self.radius, self.temp0, 'b-', label='temp')
            ax.semilogy(self.radius, self.rho0, 'r-', label='rho')
            ax.set_ylabel('Reference state')
            ax.set_xlabel('Radius')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)

            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.beta, 'b-', label='beta')
            ax.plot(self.radius, self.dbeta, 'r-', label='dbeta/dr')
            ax.set_xlabel('Radius')
            ax.set_ylabel('Derivatives of rho')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)

            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.grav)
            ax.set_xlabel('Radius')
            ax.set_xlabel('Gravity')
            ax.set_xlim(self.radius.min(), self.radius.max())
        elif self.name == 'varDiff':
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(self.radius, self.conduc, 'b-', label='conductivity')
            ax.semilogy(self.radius, self.kappa, 'r-', label='diffusivity')
            ax.semilogy(self.radius, self.prandtl, 'g-', label='Prandtl')
            ax.set_ylabel('Thermal properties')
            ax.set_xlabel('Radius')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
            
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.dLkappa, 'b-', label='dLkappa')
            ax.set_xlabel('Radius')
            ax.set_ylabel('$d\ln\kappa / dr$')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
        elif self.name == 'varVisc':
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(self.radius, self.dynVisc, 'b-', label='dyn. visc')
            ax.semilogy(self.radius, self.kinVisc, 'r-', label='kin. visc')
            ax.semilogy(self.radius, self.prandtl, 'g-', label='Prandtl')
            ax.semilogy(self.radius, self.ekman, 'k-', label='Ekman')
            ax.set_ylabel('Thermal properties')
            ax.set_xlabel('Radius')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
            
            fig = P.figure()
            ax = fig.add_subplot(111)
            axplot(self.radius, self.dLvisc, 'b-', label='dLvisc')
            ax.set_xlabel('Radius')
            ax.set_ylabel(r'$d\ln\nu / dr$')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
        elif self.name == 'varCond':
            P.figure()
            P.semilogy(self.radius, self.conduc, 'b-', label='conductivity')
            P.ylabel('Electrical conductivity')
            P.xlabel('Radius')
            P.xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
        elif self.name == 'powerR':
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.buoPower, 'b-', label='Power')
            ax.plot(self.radius, self.viscDiss, 'g-', label='visc diss')
            if self.ohmDiss.max() != 0.:
                ax.plot(self.radius, self.ohmDiss, 'r-', label='ohm diss')
            ax.set_xlabel('Radius')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
        elif self.name == 'parrad' or self.name == 'bLayersR':
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.entropy, 'b-', label='entropy')
            ax.twinx()
            ax.plot(self.radius, self.varS, 'g-', label='entropy variance')
            ax.set_xlabel('Radius')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)

            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.uh, 'b-', label='uh')
            ax.plot(self.radius, self.duhdr, 'g-', label='duhdr')
            ax.set_xlabel('Radius')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)
        elif self.name == 'parR':
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.radius, self.rol, label='Rol')
            ax.plot(self.radius, self.urol, label='u Rol')
            ax.set_xlabel('Radius')
            ax.set_ylabel('Rol')
            ax.set_xlim(self.radius.min(), self.radius.max())
            ax.legend(loc='best', frameon=False)

        if hasattr(self, 'con_radratio'):
            if self.nVarCond == 2:
                ax.axvline(self.con_radratio*self.radius[0], color='k',
                          linestyle='--')
