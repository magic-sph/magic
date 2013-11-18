# -*- coding: utf-8 -*-
import pylab as P
import numpy as N
from magic import MagicRadial, matder, intcheb, MagicSetup, scanDir
from magic.setup import labTex
from scipy.signal import argrelextrema
from scipy.integrate import simps

__author__  = "$Author$"
__date__   = "$Date$"
__version__ = "$Revision$"

def integBulkBc(rad, field, ri, ro, lambdai, lambdao):
    # Dissipation in the boundary layers
    mask = (rad >= ro-lambdao)
    y = field[mask]
    x = rad[mask]
    outerB = simps(y, x)
    mask = (rad<=ri+lambdai)
    y = field[mask]
    x = rad[mask]
    innerB = simps(y, x)
    integBc = -(outerB+innerB)

    # Dissipation in the bulk
    mask = (rad <= ro-lambdao)*(rad>=ri+lambdai)
    y = field[mask]
    x = rad[mask]
    integBulk = -simps(y, x)

    return integBc, integBulk

def getMaxima(rr, field):
    maxS = []
    for k in range(len(field)):
        if k > 3 and k < len(field)-3:
            if field[k] > field[k-1] and field[k] > field[k+1]:
                maxS.append(rr[k])
    return maxS

class BLayers(MagicSetup):

    def __init__(self, iplot=False, quiet=False):
        """
        :param iplot: a boolean to toggle plotting
        :param quiet: a boolean to (not) display the output
        """
        logFiles = scanDir('log.*')
        MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])
        par = MagicRadial(field='bLayersR', iplot=False)
        self.varS = par.varS
        self.ss = par.entropy
        self.uh = par.uh
        self.duh = par.duhdr
        self.rad = par.radius
        self.ro = self.rad[0]
        self.ri = self.rad[-1]
        if hasattr(par, 'dissS'):
            self.dissS = par.dissS
            self.epsT = -4.*N.pi*intcheb(self.rad**2*self.dissS, len(self.rad)-1, self.ro,
                                         self.ri)

            self.epsTR = 4.*N.pi*self.rad**2*self.dissS
            rrMax = getMaxima(self.rad, -abs(self.epsTR-self.epsT))

            self.dissTopS = self.ro-rrMax[0]
            self.dissBotS = rrMax[-1]-self.ri

            self.dissEpsTbl, self.dissEpsTbulk = integBulkBc(self.rad, self.epsTR, 
                             self.ri, self.ro, self.dissBotS, self.dissTopS)

            print 'thDiss bl, bulk',  self.dissEpsTbl/self.epsT, self.dissEpsTbulk/self.epsT
        # First way of defining the thermal boundary layers: with var(S)
        #rThLayer = getMaxima(self.rad, self.varS)
        ind = argrelextrema(self.varS, N.greater)[0]
        if len(ind) != 0:
            self.bcTopVarS = self.ro-self.rad[ind[0]]
            self.bcBotVarS = self.rad[ind[-1]]-self.ri
        else:
            self.bcTopVarS = 1.
            self.bcBotVarS = 1.
        if hasattr(self, 'epsT'):
            self.varSEpsTbl, self.varSEpsTbulk = integBulkBc(self.rad, self.epsTR, 
                         self.ri, self.ro, self.bcBotVarS, self.bcTopVarS)
            print 'var(S) bl, bulk', self.varSEpsTbl/self.epsT, self.varSEpsTbulk/self.epsT

        # Second way of defining the thermal boundary layers: intersection of the slopes
        d1 = matder(len(self.rad)-1, self.ro, self.ri)
        self.ttm = 3.*intcheb(self.ss*self.rad**2, len(self.rad)-1, self.ri, self.ro) \
                   /(self.ro**3-self.ri**3)
        dsdr = N.dot(d1, self.ss)
        slopeTop = dsdr[0]*(self.rad-self.ro)+self.ss[0]
        slopeBot = dsdr[-1]*(self.rad-self.ri)+self.ss[-1]

        self.dtdrm = dsdr[len(self.ss)/2]
        slopeMid = self.dtdrm*(self.rad-(self.ri+self.ro)/2.)+self.ss[len(self.ss)/2]

        mask = N.where(abs(slopeTop-self.ttm) == abs(slopeTop-self.ttm).min(), 1, 0)
        mask = N.where(abs(slopeTop-slopeMid) == abs(slopeTop-slopeMid).min(), 1, 0)
        ind = N.nonzero(mask)[0][0]
        self.bcTopSlope = self.ro-self.rad[ind]
        mask = N.where(abs(slopeBot-self.ttm) == abs(slopeBot-self.ttm).min(), 1, 0)
        mask = N.where(abs(slopeBot-slopeMid) == abs(slopeBot-slopeMid).min(), 1, 0)
        ind = N.nonzero(mask)[0][0]
        self.bcBotSlope = self.rad[ind]-self.ri

        if hasattr(self, 'epsT'):
            self.slopeEpsTbl, self.slopeEpsTbulk = integBulkBc(self.rad, self.epsTR, 
                         self.ri, self.ro, self.bcBotSlope, self.bcTopSlope)

            print 'slopes bl, bulk', self.slopeEpsTbl/self.epsT, self.slopeEpsTbulk/self.epsT
            
        pow = MagicRadial(field='powerR', iplot=False)
        self.vi = pow.viscDiss
        self.epsV = intcheb(self.vi, len(self.rad)-1, self.ri, self.ro)
        ind = argrelextrema(-abs(self.vi-self.epsV), N.greater)[0]
        rrMax = getMaxima(self.rad, -abs(self.vi-self.epsV))
        self.dissTopV = self.ro-self.rad[ind[0]]
        self.dissBotV = self.rad[ind[-1]]-self.ri
        self.dissEpsVbl, self.dissEpsVbulk = integBulkBc(self.rad, self.vi, 
                         self.ri, self.ro, self.dissBotV, self.dissTopV)
        print 'visc Diss bl, bulk', self.dissEpsVbl/self.epsV, self.dissEpsVbulk/self.epsV


        # First way of defining the viscous boundary layers: with duhdr
        #rViscousLayer = getMaxima(self.rad, self.duh)
        if self.kbotv == 1 and self.ktopv == 1:
            ind = argrelextrema(self.duh, N.greater)[0]
            if len(ind) == 0:
                self.bcTopduh = 1.
                self.bcBotduh = 1.
            else:
                if ind[0] < 4:
                    self.bcTopduh = self.ro-self.rad[ind[1]]
                else:
                    self.bcTopduh = self.ro-self.rad[ind[0]]
                if len(self.rad)-ind[-1] < 4:
                    self.bcBotduh = self.rad[ind[-2]]-self.ri
                else:
                    self.bcBotduh = self.rad[ind[-1]]-self.ri
        else:
            ind = argrelextrema(self.uh, N.greater)[0]
            if len(ind) == 0:
                self.bcTopduh = 1.
                self.bcBotduh = 1.
            else:
                if ind[0] < 4:
                    self.bcTopduh = self.ro-self.rad[ind[1]]
                else:
                    self.bcTopduh = self.ro-self.rad[ind[0]]
                if len(self.rad)-ind[-1] < 4:
                    self.bcBotduh = self.rad[ind[-2]]-self.ri
                else:
                    self.bcBotduh = self.rad[ind[-1]]-self.ri
        self.uhEpsVbl, self.uhEpsVbulk = integBulkBc(self.rad, self.vi, 
                         self.ri, self.ro, self.bcBotduh, self.bcTopduh)
        print 'uh bl, bulk', self.uhEpsVbl/self.epsV, self.uhEpsVbulk/self.epsV

        # Second way of defining the viscous boundary layers: with 
        # the viscous heating profile
        #rViscousLayer = getMaxima(self.rad, self.vi)
        if self.kbotv == 1 and self.ktopv == 1:
            ind = argrelextrema(self.vi, N.greater)[0]
            if ind[0] < 4:
                self.bcTopDiss = self.ro-self.rad[ind[1]]
            else:
                self.bcTopDiss = self.ro-self.rad[ind[0]]
            if len(self.rad)-ind[-1] < 4:
                self.bcBotDiss = self.rad[ind[-2]]-self.ri
            else:
                self.bcBotDiss = self.rad[ind[-1]]-self.ri
        else:
            ind = argrelextrema(-self.vi, N.greater)[0]
            if ind[0] < 4:
                self.bcTopDiss = self.ro-self.rad[ind[1]]
            else:
                self.bcTopDiss = self.ro-self.rad[ind[0]]
            if len(self.rad)-ind[-1] < 4:
                self.bcBotDiss = self.rad[ind[-2]]-self.ri
            else:
                self.bcBotDiss = self.rad[ind[-1]]-self.ri


        # Convective Rol in the thermal boundary Layer
        par = MagicRadial(field='parR', iplot=False)
        kin = MagicRadial(field='eKinR', iplot=False)
        ekinNas = kin.ekin_pol+kin.ekin_tor-kin.ekin_pol_axi-kin.ekin_tor_axi
        ReR = N.sqrt(2.*abs(ekinNas)/par.radius**2/(4.*N.pi))
        RolC = ReR*par.ek/par.dlVc

        y = RolC[par.radius >= self.ro-self.bcTopSlope]
        x = par.radius[par.radius >= self.ro-self.bcTopSlope]
        self.rolTop = simps(3.*y*x**2, x)/(self.ro**3-(self.ro-self.bcTopSlope)**3)

        y = RolC[par.radius <= self.ri+self.bcBotSlope]
        x = par.radius[par.radius <= self.ri+self.bcBotSlope]
        self.rolBot = simps(3.*y*x**2, x)/((self.ri+self.bcBotSlope)**3-self.ri**3)

        if iplot:
            self.plot(slopeTop, slopeBot, slopeMid)

        if not quiet:
            print self

    def plot(self, slopeTop, slopeBot, slopeMid):
        #P.rcdefaults()
        fig = P.figure()
        ax = fig.add_subplot(211)
        ax.plot(self.rad, self.ss)
        ax.axhline(self.ttm, color='gray', linestyle='-')
        ax.plot(self.rad, slopeTop, 'k--')
        ax.plot(self.rad, slopeBot, 'k--')
        ax.plot(self.rad, slopeMid, 'k--')
        ax.set_ylim(self.ss[0], self.ss[-1])
        ax.set_ylabel('Entropy')
        ax1 = ax.twinx()
        ax1.plot(self.rad, self.varS/self.varS.max(), 'g-')
        ax1.axvline(self.ro-self.bcTopVarS, color='k', linestyle=':')
        ax1.axvline(self.ri+self.bcBotVarS, color='k', linestyle=':')
        ax1.set_ylim(0, 1)
        ax1.set_ylabel('var(s)')
        ax.set_xlabel('Radius')
        ax.set_xlim(self.rad[-1], self.rad[0])

        ax = fig.add_subplot(212)
        if self.kbotv == 1 and self.ktopv == 1:
            ax.plot(self.rad, self.duh/self.duh.max())
            if labTex:
                ax.set_ylabel(r'$\partial u_h/\partial r$')
            else:
                ax.set_ylabel('duh/dr')
        else:
            ax.plot(self.rad, self.uh/self.uh.max())
            if labTex:
                ax.set_ylabel(r'$u_h$')
            else:
                ax.set_ylabel('uh')
        ax.axvline(self.ro-self.bcTopduh, color='k', linestyle='--')
        ax.axvline(self.ri+self.bcBotduh, color='k', linestyle='--')
        ax.plot(self.rad, self.vi/self.vi.max())
        ax.axvline(self.ro-self.bcTopDiss, color='k', linestyle=':')
        ax.axvline(self.ri+self.bcBotDiss, color='k', linestyle=':')
        ax.set_xlim(self.rad[-1], self.rad[0])
        ax.set_xlabel('Radius')

        fig = P.figure()
        ax = fig.add_subplot(211)
        ax.semilogy(self.rad, self.vi)
        ax.axhline(self.epsV, color='k', linestyle='--')
        ax.axvline(self.ro-self.dissTopV, color='k', linestyle='--')
        ax.axvline(self.ri+self.dissBotV, color='k', linestyle='--')
        ax.set_xlim(self.rad[-1], self.rad[0])
        ax.set_xlabel('Radius')
        ax.set_ylabel('Viscous dissipation')
        if hasattr(self, 'dissS'):
            ax = fig.add_subplot(212)
            ax.semilogy(self.rad, self.epsTR)
            ax.axhline(self.epsT, color='k', linestyle=':')
            ax.axvline(self.ro-self.dissTopS, color='k', linestyle='--')
            ax.axvline(self.ri+self.dissBotS, color='k', linestyle='--')
            ax.set_xlim(self.rad[-1], self.rad[0])
            ax.set_xlabel('Radius')
            ax.set_ylabel('Thermal Dissipation')

    def __str__(self):
        if self.ek == -1:
            ek = 0. # to avoid the -1 for the non-rotating cases
        else:
            ek = self.ek
        if self.mode == 0:
            st ='.3e%9.2e%9.2e%9.2e%5.2f' % (self.ra, ek, self.pr, self.prmag, 
                                             self.strat)
        else:
            st = '%.3e%12.5e%5.2f%6.2f' % (self.ra, ek, self.strat, self.pr)

        st += '%12.5e%12.5e%12.5e%12.5e' % (self.bcTopVarS, self.bcTopSlope,
                                            self.bcBotVarS, self.bcBotSlope)
        st += '%12.5e%12.5e%12.5e%12.5e' % (self.bcTopduh, self.bcTopDiss,
                                            self.bcBotduh, self.bcBotDiss)
        st += '%12.5e%12.5e' % (abs(self.rolTop), abs(self.rolBot))
        return st


if __name__ == '__main__':

    b = BLayers(iplot=True)
    print b
    P.show()
