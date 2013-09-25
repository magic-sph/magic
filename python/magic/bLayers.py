# -*- coding: utf-8 -*-
import pylab as P
import numpy as N
from magic import MagicRadial, matder, intcheb, MagicSetup, scanDir
from scipy.signal import argrelextrema
from scipy.integrate import simps

def getMaxima(rr, field):
    maxS = []
    for k in range(len(field[1:])):
        if field[k] > field[k-1] and field[k] > field[k+1]:
            maxS.append(rr[k])
    return maxS

class BLayers(MagicSetup):

    def __init__(self, iplot=False, quiet=False):
        logFiles = scanDir('log.*')
        MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])
        par = MagicRadial(field='bLayersR', iplot=False)
        self.varS = par.varS
        self.ss = par.entropy
        self.duh = par.duhdr
        self.rad = par.radius
        self.ro = self.rad[0]
        self.ri = self.rad[-1]
        pow = MagicRadial(field='powerR', iplot=False)
        self.vi = pow.viscDiss

        # First way of defining the thermal boundary layers: with var(S)
        #rThLayer = getMaxima(self.rad, self.varS)
        ind = argrelextrema(self.varS, N.greater)[0]
        self.bcTopVarS = self.ro-self.rad[ind[0]]
        self.bcBotVarS = self.rad[ind[-1]]-self.ri

        # Second way of defining the thermal boundary layers: intersection of the slopes
        d1 = matder(len(self.rad)-1, self.ro, self.ri)
        self.ttm = 3.*intcheb(self.ss*self.rad**2, len(self.rad)-1, self.ri, self.ro) \
                   /(self.ro**3-self.ri**3)
        dsdr = N.dot(d1, self.ss)
        slopeTop = dsdr[0]*(self.rad-self.ro)+self.ss[0]
        slopeBot = dsdr[-1]*(self.rad-self.ri)+self.ss[-1]

        mask = N.where(abs(slopeTop-self.ttm) == abs(slopeTop-self.ttm).min(), 1, 0)
        ind = N.nonzero(mask)[0][0]
        self.bcTopSlope = self.ro-self.rad[ind]
        mask = N.where(abs(slopeBot-self.ttm) == abs(slopeBot-self.ttm).min(), 1, 0)
        ind = N.nonzero(mask)[0][0]
        self.bcBotSlope = self.rad[ind]-self.ri

        # First way of defining the viscous boundary layers: with duhdr
        #rViscousLayer = getMaxima(self.rad, self.duh)
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

        # Second way of defining the viscous boundary layers: with the viscous heating profile
        #rViscousLayer = getMaxima(self.rad, self.vi)
        ind = argrelextrema(self.vi, N.greater)[0]
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
        ReR = N.sqrt(2.*ekinNas/par.radius**2/(4.*N.pi))
        RolC = ReR*par.ek/par.dlVc

        y = RolC[par.radius >= self.ro-self.bcTopSlope]
        x = par.radius[par.radius >= self.ro-self.bcTopSlope]
        self.rolTop = simps(3.*y*x**2, x)/(self.ro**3-(self.ro-self.bcTopSlope)**3)

        y = RolC[par.radius <= self.ri+self.bcBotSlope]
        x = par.radius[par.radius <= self.ri+self.bcBotSlope]
        self.rolBot = simps(3.*y*x**2, x)/((self.ri+self.bcBotSlope)**3-self.ri**3)

        if iplot:
            self.plot(slopeTop, slopeBot)

        if not quiet:
            print self

    def plot(self, slopeTop, slopeBot):
        #P.rcdefaults()
        fig = P.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.rad, self.ss)
        ax.axhline(self.ttm, color='k', linestyle='--')
        ax.plot(self.rad, slopeTop, 'k--')
        ax.plot(self.rad, slopeBot, 'k--')
        ax.set_ylim(self.ss[0], self.ss[-1])
        ax.set_ylabel('Entropy', fontsize=18)
        ax1 = ax.twinx()
        ax1.plot(self.rad, self.varS/self.varS.max(), 'g-')
        ax1.axvline(self.ro-self.bcTopVarS, color='k', linestyle=':')
        ax1.axvline(self.ri+self.bcBotVarS, color='k', linestyle=':')
        ax1.set_ylim(0, 1)
        ax1.set_ylabel('var(s)', fontsize=18)
        ax.set_xlabel('Radius', fontsize=18)
        ax.set_xlim(self.rad[-1], self.rad[0])

        fig1 = P.figure()
        ax = fig1.add_subplot(111)
        ax.plot(self.rad, self.duh/self.duh.max())
        ax.axvline(self.ro-self.bcTopduh, color='k', linestyle='--')
        ax.axvline(self.ri+self.bcBotduh, color='k', linestyle='--')
        ax.plot(self.rad, self.vi/self.vi.max())
        ax.axvline(self.ro-self.bcTopDiss, color='k', linestyle=':')
        ax.axvline(self.ri+self.bcBotDiss, color='k', linestyle=':')
        ax.set_xlim(self.rad[-1], self.rad[0])
        ax.set_xlabel('Radius', fontsize=18)
        ax.set_ylabel('duh/dr', fontsize=18)

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
