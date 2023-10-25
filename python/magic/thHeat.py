# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from magic import scanDir, MagicSetup, Movie, chebgrid, rderavg, AvgField
import os, pickle
from scipy.integrate import simps, trapz

json_model = {'phys_params': [],
              'time_series': { 'heat': ['topnuss', 'botnuss']},
              'spectra': {},
              'radial_profiles': {}}


class ThetaHeat(MagicSetup):
    """
    This class allows to conduct some analysis of the latitudinal
    variation of the heat transfer. It relies on the movie files
    :ref:`ATmov.TAG <secMovieFile>` and :ref:`AHF_mov.TAG <secMovieFile>`.
    As it's a bit time-consuming, the calculations are stored in a
    python.pickle file to quicken future usage of the data.

    Since this function is supposed to use time-averaged quantities, the usual
    procedure is first to define the initial averaging time using
    :py:class:`AvgField <magic.AvgField>`: (this needs to be done only once)

    >>> a = AvgField(tstart=2.58)

    Once the ``tInitAvg`` file exists, the latitudinal heat transfer analysis
    can be done using:

    >>> # For chunk-averages over 10^\degree in the polar and equatorial regions.
    >>> th = ThetaHeat(angle=10)
    >>> # Formatted output
    >>> print(th)
    """

    def __init__(self, iplot=False, angle=10, pickleName='thHeat.pickle',
                 quiet=False):
        """
        :param iplot: a boolean to toggle the plots on/off
        :type iplot: bool
        :param angle: the integration angle in degrees
        :type angle: float
        :pickleName: calculations a
        :param quiet: a boolean to switch on/off verbose outputs
        :type quiet: bool
        """

        angle = angle * np.pi / 180

        if os.path.exists('tInitAvg'):
            f = open('tInitAvg', 'r')
            tstart = float(f.readline())
            f.close()
            logFiles = scanDir('log.*')
            tags = []
            for lg in logFiles:
                nml = MagicSetup(quiet=True, nml=lg)
                if nml.start_time >  tstart:
                    if os.path.exists('ATmov.{}'.format(nml.tag)):
                        tags.append(nml.tag)
            if len(tags) == 0:
                tags = [nml.tag]
                if not quiet:
                    print('Only 1 tag: {}'.format(tags))
            MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])

            a = AvgField(model=json_model, write=False)
            self.nuss = 0.5 * (a.topnuss_av+a.botnuss_av)
        else:
            logFiles = scanDir('log.*')
            if len(logFiles) > 0:
                MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])

        if not os.path.exists(pickleName):
            # reading ATmov
            k = 0
            for tag in tags:
                f = 'ATmov.{}'.format(tag)
                if os.path.exists(f):
                    if k == 0:
                        m = Movie(file=f, iplot=False)
                        print(f)
                    else:
                        m += Movie(file=f, iplot=False)
                        print(f)
                    k += 1

            # reading AHF_mov
            kk = 0
            for tag in tags:
                f = 'AHF_mov.{}'.format(tag)
                if os.path.exists(f):
                    if kk == 0:
                        m1 = Movie(file=f, iplot=False)
                        print(f)
                    else:
                        m1 += Movie(file=f, iplot=False)
                        print(f)
                    kk += 1

            self.tempmean = m.data[0, ...].mean(axis=0)
            self.tempstd = m.data[0, ...].std(axis=0)
            self.colat = m.theta


            if kk > 0: # i.e. at least one AHF_mov file has been found
                self.fluxmean = m1.data[0, ...].mean(axis=0)
                self.fluxstd = m1.data[0, ...].std(axis=0)
            else:
                self.fluxmean = rderavg(self.tempmean, m.radius, exclude=False)
                self.fluxstd = rderavg(self.tempstd, m.radius, exclude=False)

            # Pickle saving
            try:
                with open(pickleName, 'wb') as f:
                    pickle.dump([self.colat, self.tempmean, self.tempstd,
                                 self.fluxmean, self.fluxstd], f)
            except PermissionError:
                print('No write access in the current directory')
        else:
            with open(pickleName, 'rb') as f:
                dat = pickle.load(f)
                if len(dat) == 5:
                    self.colat, self.tempmean, self.tempstd, \
                                self.fluxmean, self.fluxstd = dat
                else:
                    self.colat, self.tempmean, self.fluxmean = dat
                    self.fluxstd = np.zeros_like(self.fluxmean)
                    self.tempstd = np.zeros_like(self.fluxmean)

        self.ri = self.radratio/(1.-self.radratio)
        self.ro = 1./(1.-self.radratio)

        self.ntheta, self.nr = self.tempmean.shape
        if not hasattr(self, 'radial_scheme') or \
           (self.radial_scheme=='CHEB' and self.l_newmap==False):
            # Redefine to get double precision
            self.radius = chebgrid(self.nr-1, self.ro, self.ri)
        else:
            self.radius = m.radius
        th2D = np.zeros((self.ntheta, self.nr), dtype=self.radius.dtype)
        #self.colat = np.linspace(0., np.pi, self.ntheta)

        for i in range(self.ntheta):
            th2D[i, :] = self.colat[i]

        self.temprmmean = 0.5*simps(self.tempmean*np.sin(th2D), th2D, axis=0)
        self.temprmstd = 0.5*simps(self.tempstd*np.sin(th2D), th2D, axis=0)
        sinTh = np.sin(self.colat)

        # Conducting temperature profile (Boussinesq only!)
        if self.ktops == 1 and self.kbots == 1:
            self.tcond = self.ri*self.ro/self.radius-self.ri+self.temprmmean[0]
            self.fcond = -self.ri*self.ro/self.radius**2
        elif self.ktops == 1 and self.kbots != 1:
            qbot = -1.
            ttop = self.temprmmean[0]
            self.fcond = -self.ri**2 / self.radius**2
            self.tcond = self.ri**2/self.radius - self.ri**2/self.ro + ttop
        else:
            if os.path.exists('pscond.dat'):
                dat = np.loadtxt('pscond.dat')
                self.tcond = dat[:, 1]
                self.fcond = rderavg(self.tcond, self.radius)

        self.nusstopmean = self.fluxmean[:, 0] / self.fcond[0]
        self.nussbotmean = self.fluxmean[:, -1] / self.fcond[-1]
        self.nusstopstd = self.fluxstd[:, 0] / self.fcond[0]
        self.nussbotstd = self.fluxstd[:, -1] / self.fcond[-1]


        # Close to the equator
        mask2D = (th2D>=np.pi/2.-angle/2.)*(th2D<=np.pi/2+angle/2.)
        mask = (self.colat>=np.pi/2.-angle/2.)*(self.colat<=np.pi/2+angle/2.)
        fac = 1./simps(sinTh[mask], self.colat[mask])
        self.nussBotEq = fac*simps(self.nussbotmean[mask]*sinTh[mask], self.colat[mask])
        self.nussTopEq = fac*simps(self.nusstopmean[mask]*sinTh[mask], self.colat[mask])
        sinC = sinTh.copy()
        sinC[~mask] = 0.
        fac = 1./simps(sinC, self.colat)
        tempC = self.tempmean.copy()
        tempC[~mask2D] = 0.
        self.tempEqmean = fac*simps(tempC*np.sin(th2D), th2D, axis=0)
        tempC = self.tempstd.copy()
        tempC[~mask2D] = 0.
        self.tempEqstd = fac*simps(tempC*np.sin(th2D), th2D, axis=0)

        dtempEq = rderavg(self.tempEqmean, self.radius)
        self.betaEq = dtempEq[len(dtempEq)//2]

        # 45\deg inclination
        mask2D = (th2D>=np.pi/4.-angle/2.)*(th2D<=np.pi/4+angle/2.)
        mask = (self.colat>=np.pi/4.-angle/2.)*(self.colat<=np.pi/4+angle/2.)
        fac = 1./simps(np.sin(self.colat[mask]), self.colat[mask])
        nussBot45NH = fac*simps(self.nussbotmean[mask]*sinTh[mask], self.colat[mask])
        nussTop45NH = fac*simps(self.nusstopmean[mask]*sinTh[mask], self.colat[mask])
        sinC = sinTh.copy()
        sinC[~mask] = 0.
        fac = 1./simps(sinC, self.colat)
        tempC = self.tempmean.copy()
        tempC[~mask2D] = 0.
        temp45NH = fac*simps(tempC*np.sin(th2D), th2D, axis=0)

        mask2D = (th2D>=3.*np.pi/4.-angle/2.)*(th2D<=3.*np.pi/4+angle/2.)
        mask = (self.colat>=3.*np.pi/4.-angle/2.)*(self.colat<=3.*np.pi/4+angle/2.)
        fac = 1./simps(np.sin(self.colat[mask]), self.colat[mask])
        nussBot45SH = fac*simps(self.nussbotmean[mask]*sinTh[mask], self.colat[mask])
        nussTop45SH = fac*simps(self.nusstopmean[mask]*sinTh[mask], self.colat[mask])
        sinC = sinTh.copy()
        sinC[~mask] = 0.
        fac = 1./simps(sinC, self.colat)
        tempC = self.tempmean.copy()
        tempC[~mask2D] = 0.
        temp45SH = fac*simps(tempC*np.sin(th2D), th2D, axis=0)

        self.nussTop45 = 0.5*(nussTop45NH+nussTop45SH)
        self.nussBot45 = 0.5*(nussBot45NH+nussBot45SH)
        self.temp45 = 0.5*(temp45NH+temp45SH)

        dtemp45 = rderavg(self.temp45, self.radius)
        self.beta45 = dtemp45[len(dtemp45)//2]

        # Polar regions
        mask2D = (th2D<=angle/2.)
        mask = (self.colat<=angle/2.)
        fac = 1./simps(np.sin(self.colat[mask]), self.colat[mask])
        nussBotPoNH = fac*simps(self.nussbotmean[mask]*sinTh[mask], self.colat[mask])
        nussTopPoNH = fac*simps(self.nusstopmean[mask]*sinTh[mask], self.colat[mask])
        sinC = sinTh.copy()
        sinC[~mask] = 0.
        fac = 1./simps(sinC, self.colat)
        tempC = self.tempmean.copy()
        tempC[~mask2D] = 0.
        tempPolNHmean = fac*simps(tempC*np.sin(th2D), th2D, axis=0)
        tempC = self.tempstd.copy()
        tempC[~mask2D] = 0.
        tempPolNHstd = fac*simps(tempC*np.sin(th2D), th2D, axis=0)

        mask2D = (th2D>=np.pi-angle/2.)
        mask = (self.colat>=np.pi-angle/2.)
        fac = 1./simps(np.sin(self.colat[mask]), self.colat[mask])
        nussBotPoSH = fac*simps(self.nussbotmean[mask]*sinTh[mask], self.colat[mask])
        nussTopPoSH = fac*simps(self.nusstopmean[mask]*sinTh[mask], self.colat[mask])
        sinC = sinTh.copy()
        sinC[~mask] = 0.
        fac = 1./simps(sinC, self.colat)
        tempC = self.tempmean.copy()
        tempC[~mask2D] = 0.
        tempPolSHmean = fac*simps(tempC*np.sin(th2D), th2D, axis=0)
        tempC = self.tempstd.copy()
        tempC[~mask2D] = 0.
        tempPolSHstd = fac*simps(tempC*np.sin(th2D), th2D, axis=0)

        self.nussBotPo = 0.5*(nussBotPoNH+nussBotPoSH)
        self.nussTopPo = 0.5*(nussTopPoNH+nussTopPoSH)
        self.tempPolmean = 0.5*(tempPolNHmean+tempPolSHmean)
        self.tempPolstd= 0.5*(tempPolNHstd+tempPolSHstd)

        dtempPol = rderavg(self.tempPolmean, self.radius)
        self.betaPol = dtempPol[len(dtempPol)//2]

        # Inside and outside TC
        angleTC = np.arcsin(self.ri/self.ro)
        mask2D = (th2D<=angleTC)
        mask = (self.colat<=angleTC)
        fac = 1./simps(np.sin(self.colat[mask]), self.colat[mask])
        nussITC_NH = fac*simps(self.nusstopmean[mask]*sinTh[mask], self.colat[mask])

        mask2D = (th2D>=np.pi-angleTC)
        mask = (self.colat>=np.pi-angleTC)
        fac = 1./simps(np.sin(self.colat[mask]), self.colat[mask])
        nussITC_SH = fac*simps(self.nusstopmean[mask]*sinTh[mask], self.colat[mask])

        self.nussITC = 0.5*(nussITC_NH+nussITC_SH)

        mask2D = (th2D>=angleTC)*(th2D<=np.pi-angleTC)
        mask = (self.colat>=angleTC)*(self.colat<=np.pi-angleTC)
        fac = 1./simps(sinTh[mask], self.colat[mask])
        self.nussOTC = fac*simps(self.nusstopmean[mask]*sinTh[mask], self.colat[mask])

        if iplot:
            self.plot()

        if not quiet:
            print(self)


    def __str__(self):
        """
        Formatted outputs

        >>> th = ThetaHeat()
        >>> print(th)
        """
        if self.ek == -1:
            ek = 0. # to avoid the -1 for the non-rotating cases
        else:
            ek = self.ek
        if self.mode == 0:
            st ='{:9.3e}{:9.2e}{:9.2e}{:9.2e}{:5.2f}'.format(self.ra, ek,
                self.pr, self.prmag, self.strat)
        else:
            st = '{:.3e}{:12.5e}{:5.2f}{:6.2f}{:6.2f}'.format(self.ra, ek,
                self.strat, self.pr, self.radratio)

        st += '{:12.5e}'.format(self.nuss)
        st += '{:12.5e}{:12.5e}{:12.5e}'.format(self.nussBotEq, self.nussBot45,
                                                self.nussBotPo)
        st += '{:12.5e}{:12.5e}{:12.5e}'.format(self.nussTopEq, self.nussTop45,
                                                self.nussTopPo)
        st += ' {:12.5e} {:12.5e} {:12.5e}'.format(self.betaEq, self.beta45,
                                                   self.betaPol)

        return st


    def plot(self):
        """
        Plotting function
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.radius, self.tcond-self.tcond[0], 'k--', label='Cond. temp.')
        ax.fill_between(self.radius,
                        self.temprmmean-self.temprmstd-self.temprmmean[0],
                        self.temprmmean+self.temprmstd-self.temprmmean[0],
                        alpha=0.2)
        ax.plot(self.radius, self.temprmmean-self.temprmmean[0], ls='-', c='C0',
                lw=2, label='Mean temp.')
        ax.fill_between(self.radius,
                        self.tempEqmean-self.tempEqstd-self.temprmmean[0],
                        self.tempEqmean+self.tempEqstd-self.temprmmean[0],
                        alpha=0.2)
        ax.plot(self.radius, self.tempEqmean-self.temprmmean[0],
                ls='-.', c='C1', lw=2, label='Temp. equat')
        ax.fill_between(self.radius,
                        self.tempPolmean-self.tempPolstd-self.temprmmean[0],
                        self.tempPolmean+self.tempPolstd-self.temprmmean[0],
                        alpha=0.2)
        ax.plot(self.radius, self.tempPolmean-self.temprmmean[0],
                ls='--', c='C2', lw=2, label='Temp. Pol')
        ax.set_xlim(self.ri, self.ro)
        ax.set_ylim(0., 1.)
        ax.set_ylabel('T')
        ax.set_xlabel('r')
        ax.legend(loc='upper right', frameon=False)

        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax1.fill_between(self.colat*180./np.pi, self.nusstopmean-self.nusstopstd,
                         self.nusstopmean+self.nusstopstd, alpha=0.2)
        ax1.plot(self.colat*180./np.pi, self.nusstopmean, ls='-', color='C0',
                 lw=2, label='Top Nu')

        ax1.fill_between(self.colat*180./np.pi, self.nussbotmean-self.nussbotstd,
                         self.nussbotmean+self.nussbotstd, alpha=0.2)
        ax1.plot(self.colat*180./np.pi, self.nussbotmean, ls='--', c='C1',
                 lw=2, label='Bot Nu')
        ax1.set_xlim(0., 180.)
        ax1.set_ylabel('Nu')
        ax1.set_xlabel('Theta')
        ax1.legend(loc='upper right', frameon=False)
        ax1.axvline(180./np.pi*np.arcsin(self.ri/self.ro), color='k', linestyle='--')
        ax1.axvline(180-180./np.pi*np.arcsin(self.ri/self.ro), color='k', linestyle='--')



if __name__ == '__main__':
    t = ThetaHeat(iplot=True)
    plt.show()
