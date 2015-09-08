# -*- coding: utf-8 -*-
import matplotlib.pyplot as P
import numpy as N
from magic import scanDir, MagicSetup, Movie, matder, chebgrid, rderavg, AvgField
import os, pickle
from scipy.integrate import simps, trapz

__author__  = "$Author$"
__date__   = "$Date$"
__version__ = "$Revision$"

class ThetaHeat(MagicSetup):
    """
    This routines allows to conduct some analysis of the latitudinal
    variation of the heat transfer
    """

    def __init__(self, iplot=False, angle=10, pickleName='thHeat.pickle'):
        """
        :param iplot: a boolean to toggle the plots on/off
        :param angle: the integration angle
        """

        angle = angle * N.pi / 180

        if os.path.exists('tInitAvg'):
            file = open('tInitAvg', 'r')
            tstart = float(file.readline())
            file.close()
            logFiles = scanDir('log.*')
            tags = []
            for lg in logFiles:
                nml = MagicSetup(quiet=True, nml=lg)
                if nml.start_time >  tstart:
                    if os.path.exists('bLayersR.%s' % nml.tag):
                        tags.append(nml.tag)
            if len(tags) == 0:
                tags = [nml.tag]
                print('Only 1 tag: %s' % tags)
            MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])

            a = AvgField()
            self.nuss = a.nuss
        else:
            logFiles = scanDir('log.*')
            MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])
        

        if not os.path.exists(pickleName):
            # reading ATmov
            k = 0
            for tag in tags:
                file = 'ATmov.%s' % tag
                if os.path.exists(file):
                    if k == 0:
                        m = Movie(file=file, iplot=False)
                        print(file)
                    else:
                        m += Movie(file=file, iplot=False)
                        print(file)
                    k += 1

            # reading AHF_mov
            kk = 0
            for tag in tags:
                file = 'AHF_mov.%s' % tag
                if os.path.exists(file):
                    if kk == 0:
                        m1 = Movie(file=file, iplot=False)
                        print(file)
                    else:
                        m1 += Movie(file=file, iplot=False)
                        print(file)
                    kk += 1

            self.tempmean = m.data.mean(axis=0)
            self.colat = m.theta


            if kk > 0: # i.e. at least one AHF_mov file has been found
                self.flux = m1.data.mean(axis=0)
            else:
                self.flux = rderavg(self.tempmean, eta=self.radratio, 
                                exclude=False, spectral=False)

            # Pickle saving
            file = open(pickleName, 'wb')
            pickle.dump([self.colat, self.tempmean, self.flux], file)
            file.close()
        else:
            file = open(pickleName, 'r')
            self.colat, self.tempmean, self.flux = pickle.load(file)
            file.close()

        self.ri = self.radratio/(1.-self.radratio)
        self.ro = 1./(1.-self.radratio)

        self.ntheta, self.nr = self.tempmean.shape
        self.radius = chebgrid(self.nr-1, self.ro, self.ri)
        th2D = N.zeros((self.ntheta, self.nr), dtype=self.radius.dtype)
        #self.colat = N.linspace(0., N.pi, self.ntheta)

        for i in range(self.ntheta):
            th2D[i, :] = self.colat[i]

        self.temprm = 0.5*simps(self.tempmean*N.sin(th2D), th2D, axis=0)
        sinTh = N.sin(self.colat)
        d1 = matder(self.nr-1, self.ro, self.ri)

        # Conducting temperature profile (Boussinesq only!)
        self.tcond = self.ri*self.ro/self.radius-self.ri+self.temprm[0]
        self.fcond = -self.ri*self.ro/self.radius**2
        self.nusstop = self.flux[:, 0] / self.fcond[0]
        self.nussbot = self.flux[:, -1] / self.fcond[-1]

        # Close to the equator
        mask2D = (th2D>=N.pi/2.-angle/2.)*(th2D<=N.pi/2+angle/2.)
        mask = (self.colat>=N.pi/2.-angle/2.)*(self.colat<=N.pi/2+angle/2.)
        fac = 1./simps(sinTh[mask], self.colat[mask])
        self.nussBotEq = fac*simps(self.nussbot[mask]*sinTh[mask], self.colat[mask])
        self.nussTopEq = fac*simps(self.nusstop[mask]*sinTh[mask], self.colat[mask])
        sinC = sinTh.copy()
        sinC[~mask] = 0.
        fac = 1./simps(sinC, self.colat)
        tempC = self.tempmean.copy()
        tempC[~mask2D] = 0.
        self.tempEq = fac*simps(tempC*N.sin(th2D), th2D, axis=0)

        dtempEq = N.dot(d1, self.tempEq)
        self.betaEq = dtempEq[self.nr/2]

        # 45\deg inclination
        mask2D = (th2D>=N.pi/4.-angle/2.)*(th2D<=N.pi/4+angle/2.)
        mask = (self.colat>=N.pi/4.-angle/2.)*(self.colat<=N.pi/4+angle/2.)
        fac = 1./simps(N.sin(self.colat[mask]), self.colat[mask])
        nussBot45NH = fac*simps(self.nussbot[mask]*sinTh[mask], self.colat[mask])
        nussTop45NH = fac*simps(self.nusstop[mask]*sinTh[mask], self.colat[mask])
        sinC = sinTh.copy()
        sinC[~mask] = 0.
        fac = 1./simps(sinC, self.colat)
        tempC = self.tempmean.copy()
        tempC[~mask2D] = 0.
        temp45NH = fac*simps(tempC*N.sin(th2D), th2D, axis=0)

        mask2D = (th2D>=3.*N.pi/4.-angle/2.)*(th2D<=3.*N.pi/4+angle/2.)
        mask = (self.colat>=3.*N.pi/4.-angle/2.)*(self.colat<=3.*N.pi/4+angle/2.)
        fac = 1./simps(N.sin(self.colat[mask]), self.colat[mask])
        nussBot45SH = fac*simps(self.nussbot[mask]*sinTh[mask], self.colat[mask])
        nussTop45SH = fac*simps(self.nusstop[mask]*sinTh[mask], self.colat[mask])
        sinC = sinTh.copy()
        sinC[~mask] = 0.
        fac = 1./simps(sinC, self.colat)
        tempC = self.tempmean.copy()
        tempC[~mask2D] = 0.
        temp45SH = fac*simps(tempC*N.sin(th2D), th2D, axis=0)

        self.nussTop45 = 0.5*(nussTop45NH+nussTop45SH)
        self.nussBot45 = 0.5*(nussBot45NH+nussBot45SH)
        self.temp45 = 0.5*(temp45NH+temp45SH)

        dtemp45 = N.dot(d1, self.temp45)
        self.beta45 = dtemp45[self.nr/2]

        # Polar regions
        mask2D = (th2D<=angle/2.)
        mask = (self.colat<=angle/2.)
        fac = 1./simps(N.sin(self.colat[mask]), self.colat[mask])
        nussBotPoNH = fac*simps(self.nussbot[mask]*sinTh[mask], self.colat[mask])
        nussTopPoNH = fac*simps(self.nusstop[mask]*sinTh[mask], self.colat[mask])
        sinC = sinTh.copy()
        sinC[~mask] = 0.
        fac = 1./simps(sinC, self.colat)
        tempC = self.tempmean.copy()
        tempC[~mask2D] = 0.
        tempPolNH = fac*simps(tempC*N.sin(th2D), th2D, axis=0)

        mask2D = (th2D>=N.pi-angle/2.)
        mask = (self.colat>=N.pi-angle/2.)
        fac = 1./simps(N.sin(self.colat[mask]), self.colat[mask])
        nussBotPoSH = fac*simps(self.nussbot[mask]*sinTh[mask], self.colat[mask])
        nussTopPoSH = fac*simps(self.nusstop[mask]*sinTh[mask], self.colat[mask])
        sinC = sinTh.copy()
        sinC[~mask] = 0.
        fac = 1./simps(sinC, self.colat)
        tempC = self.tempmean.copy()
        tempC[~mask2D] = 0.
        tempPolSH = fac*simps(tempC*N.sin(th2D), th2D, axis=0)

        self.nussBotPo = 0.5*(nussBotPoNH+nussBotPoSH)
        self.nussTopPo = 0.5*(nussTopPoNH+nussTopPoSH)
        self.tempPol = 0.5*(tempPolNH+tempPolSH)

        dtempPol = N.dot(d1, self.tempPol)
        self.betaPol = dtempPol[self.nr/2]

        if iplot:
            self.plot()

        print(self)


    def __str__(self):
        if self.ek == -1:
            ek = 0. # to avoid the -1 for the non-rotating cases
        else:
            ek = self.ek
        if self.mode == 0:
            st ='%9.3e%9.2e%9.2e%9.2e%5.2f' % (self.ra, ek, self.pr, self.prmag,
                                               self.strat)
        else:
            st = '%.3e%12.5e%5.2f%6.2f%6.2f' % (self.ra, ek, self.strat, self.pr, 
                                                self.radratio)

        st += '%12.5e' % (self.nuss)
        st += '%12.5e%12.5e%12.5e' % (self.nussBotEq, self.nussBot45, self.nussBotPo)
        st += '%12.5e%12.5e%12.5e' % (self.nussTopEq, self.nussTop45, self.nussTopPo)
        st += ' %12.5e %12.5e %12.5e' % (self.betaEq, self.beta45, self.betaPol)

        return st



    def plot(self):
        fig = P.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.radius, self.tcond-self.tcond[0], 'k--', label='Cond. temp.')
        ax.plot(self.radius, self.temprm-self.temprm[0], 'k-', lw=2, label='Mean temp.')
        ax.plot(self.radius, self.tempEq-self.temprm[0],
                'r-.', lw=2, label='Temp. equat')
        ax.plot(self.radius, self.temp45-self.temprm[0],
                'b:', lw=2, label='Temp. 45')
        ax.plot(self.radius, self.tempPol-self.temprm[0],
                'g--', lw=2, label='Temp. Pol')
        ax.set_xlim(self.ri, self.ro)
        ax.set_ylim(0., 1.)
        ax.set_ylabel('T')
        ax.set_xlabel('r')
        ax.legend(loc='upper right', frameon=False)

        fig1 = P.figure()
        ax1 = fig1.add_subplot(111)
        ax1.plot(self.colat*180./N.pi, self.nusstop, 'k-', lw=2, label='Top Nu')
        ax1.plot(self.colat*180./N.pi, self.nussbot, 'g--', lw=2, label='Top Nu')
        ax1.set_xlim(0., 180.)
        ax1.set_ylabel('Nu')
        ax1.set_xlabel('Theta')
        ax1.legend(loc='upper right', frameon=False)
        ax1.axvline(180./N.pi*N.arcsin(self.ri/self.ro), color='k', linestyle='--')
        ax1.axvline(180-180./N.pi*N.arcsin(self.ri/self.ro), color='k', linestyle='--')



if __name__ == '__main__':
    t = ThetaHeat(iplot=True)
    P.show()
