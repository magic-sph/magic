# -*- coding: utf-8 -*-
import os, re
import pylab as P
import numpy as N
from setup import MagicSetup
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
            self.emag_pols_surf = data[:, 5]
            self.emag_pol_axis_surf = data[:, 6]
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
            self.entropy = data[:, 9]
            self.varS = data[:, 10]
            self.uh = data[:, 11]
            self.duhdr = data[:, 12]
        elif self.name == 'bLayersR':
            self.radius = data[:, 0]
            self.entropy = data[:, 1]
            self.varS = data[:, 2]
            self.uh = data[:, 3]
            self.duhdr = data[:, 4]

        if iplot:
            self.plot()

    def plot(self):
        P.rc('figure.subplot', right=0.95, top=0.95)
        P.rc('axes', labelsize=20)
        P.rc('xtick', labelsize=14)
        P.rc('ytick', labelsize=14)
        if self.name == 'eKinR':
            P.figure()
            P.plot(self.radius, self.ekin_pol+self.ekin_tor, 'k-', 
                   label='ekin tot')
            P.plot(self.radius, self.ekin_pol, 'b-', label='ekin pol')
            P.plot(self.radius, self.ekin_tor, 'r-', label='ekin tor')
            P.plot(self.radius, self.ekin_pol_axi, 'b--', label='ekin pol axi')
            P.plot(self.radius, self.ekin_tor_axi, 'r--', label='ekin tor axi')
            P.xlabel('Radius')
            P.ylabel('Kinetic energy')
            P.xlim(self.radius.min(), self.radius.max())
        elif self.name == 'eMagR':
            P.figure()
            P.plot(self.radius, self.emag_pol+self.emag_tor, 'k-', 
                   label='emag tot')
            P.plot(self.radius, self.emag_pol, 'b-', label='emag pol')
            P.plot(self.radius, self.emag_tor, 'r-', label='emag tor')
            P.plot(self.radius, self.emag_pol_axi, 'b--', label='emag pol axi')
            P.plot(self.radius, self.emag_tor_axi, 'r--', label='emag tor axi')
            P.xlabel('Radius')
            P.ylabel('Magnetic energy')
            P.xlim(self.radius.min(), self.radius.max())
            if hasattr(self, 'nVarCond'):
                if self.nVarCond == 2:
                    P.axvline(self.con_radratio*self.radius[0], color='k',
                              linestyle='--')
            P.figure()
            P.plot(self.radius, self.dip_ratio)
            P.xlabel('Radius')
            P.ylabel('Dipolarity')
            P.xlim(self.radius.min(), self.radius.max())
        elif self.name == 'anel':
            P.figure()
            P.semilogy(self.radius, self.temp0, 'b-', label='temp')
            P.semilogy(self.radius, self.rho0, 'r-', label='rho')
            P.ylabel('Reference state')
            P.xlabel('Radius')
            P.xlim(self.radius.min(), self.radius.max())
            P.legend(loc='best', frameon=False)

            P.figure()
            P.plot(self.radius, self.beta, 'b-', label='beta')
            P.plot(self.radius, self.dbeta, 'r-', label='dbeta/dr')
            P.xlabel('Radius')
            P.ylabel('Derivatives of rho')
            P.xlim(self.radius.min(), self.radius.max())

            P.figure()
            P.plot(self.radius, self.grav)
            P.xlabel('Radius')
            P.ylabel('Gravity')
            P.xlim(self.radius.min(), self.radius.max())
        elif self.name == 'varDiff':
            P.figure()
            P.semilogy(self.radius, self.conduc, 'b-', label='conductivity')
            P.semilogy(self.radius, self.kappa, 'r-', label='diffusivity')
            P.semilogy(self.radius, self.prandtl, 'g-', label='Prandtl')
            P.ylabel('Thermal properties')
            P.xlabel('Radius')
            P.xlim(self.radius.min(), self.radius.max())
            P.legend(loc='best', frameon=False)
            
            P.figure()
            P.plot(self.radius, self.dLkappa, 'b-', label='dLkappa')
            P.xlabel('Radius')
            P.ylabel('$d\ln\kappa / dr$')
            P.xlim(self.radius.min(), self.radius.max())
        elif self.name == 'varVisc':
            P.figure()
            P.semilogy(self.radius, self.dynVisc, 'b-', label='dyn. visc')
            P.semilogy(self.radius, self.kinVisc, 'r-', label='kin. visc')
            P.semilogy(self.radius, self.prandtl, 'g-', label='Prandtl')
            P.semilogy(self.radius, self.ekman, 'k-', label='Ekman')
            P.ylabel('Thermal properties')
            P.xlabel('Radius')
            P.xlim(self.radius.min(), self.radius.max())
            P.legend(loc='best', frameon=False)
            
            P.figure()
            P.plot(self.radius, self.dLvisc, 'b-', label='dLvisc')
            P.xlabel('Radius')
            P.ylabel(r'$d\ln\nu / dr$')
            P.xlim(self.radius.min(), self.radius.max())
        elif self.name == 'varCond':
            P.figure()
            P.semilogy(self.radius, self.conduc, 'b-', label='conductivity')
            P.ylabel('Electrical conductivity')
            P.xlabel('Radius')
            P.xlim(self.radius.min(), self.radius.max())
        elif self.name == 'powerR':
            P.figure()
            P.plot(self.radius, self.buoPower, 'b-', label='Power')
            P.plot(self.radius, self.viscDiss, 'g-', label='visc diss')
            if self.ohmDiss.max() != 0.:
                P.plot(self.radius, self.ohmDiss, 'r-', label='ohm diss')
            P.xlabel('Radius')
            P.xlim(self.radius.min(), self.radius.max())
        elif self.name == 'parrad' or self.name == 'bLayersR':
            P.figure()
            P.plot(self.radius, self.entropy, 'b-', label='entropy')
            P.twinx()
            P.plot(self.radius, self.varS, 'g-', label='entropy variance')
            P.xlabel('Radius')
            P.xlim(self.radius.min(), self.radius.max())
            P.legend(loc='best', frameon=False)

            P.figure()
            P.plot(self.radius, self.uh, 'b-', label='uh')
            P.plot(self.radius, self.duhdr, 'g-', label='duhdr')
            P.xlabel('Radius')
            P.xlim(self.radius.min(), self.radius.max())

        if hasattr(self, 'nVarCond'):
            if self.nVarCond == 2:
                P.axvline(self.con_radratio*self.radius[0], color='k',
                          linestyle='--')
        P.legend(loc='best', frameon=False)
