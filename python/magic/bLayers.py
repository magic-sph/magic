# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import os
from magic import MagicRadial, matder, intcheb, MagicSetup, scanDir, AvgField
from magic.setup import labTex
from scipy.signal import argrelextrema
from scipy.integrate import simps
from scipy.interpolate import splrep, splev

def getAccuratePeaks(rad, uh, uhTop, uhBot, ri, ro):
    """
    This functions performs a spline extrapolation around the maxima
    of the input array uh to define a more accurate location of the
    boundary layer.

    :param rad: radius
    :type rad: numpy.ndarray
    :param uh: the horizontal velocity profile
    :type uh: numpy.ndarray
    :param uhTop: first peak value of uh close to the outer boundary
    :type uhTop: float
    :param uhBot: first peak value of uh close to the inner boundary
    :type uhBot: float
    :param ri: the inner core radius
    :type ri: float
    :param ro: the outer core radius
    :type ro: float
    :returns: four floats: thickness of the bottom boundary layer,
              thickness of the top boundary layer, extrapolated value of uh
              at the bottom boundary layer, extrapolated value of uh at the
              top boundary layer
    :rtype: list
    """
    idxT = np.nonzero(np.where(uh==uhTop, 1, 0))[0][0]
    x = rad[idxT-3:idxT+4]
    y =  uh[idxT-3:idxT+4]
    newx = np.linspace(x[0], x[-1], 128)
    tckp = splrep(x[::-1], y[::-1])
    newy = splev(newx, tckp)
    ind = argrelextrema(newy, np.greater)[0]
    bcTopUh = ro-newx[ind]
    topUh = newy[ind]

    idxT = np.nonzero(np.where(uh==uhBot, 1, 0))[0][0]
    x = rad[idxT-3:idxT+4]
    y =  uh[idxT-3:idxT+4]
    newx = np.linspace(x[0], x[-1], 128)
    tckp = splrep(x[::-1], y[::-1])
    newy = splev(newx, tckp)
    ind = argrelextrema(newy, np.greater)[0]
    bcBotUh = newx[ind]-ri
    botUh = newy[ind]

    return bcBotUh[0], bcTopUh[0], botUh[0], topUh[0]


def integBulkBc(rad, field, ri, ro, lambdai, lambdao, normed=False):
    """
    This function evaluates the radial integral of the input array field
    in the boundary layer and in the bulk separately.

    :param rad: radius
    :type rad: numpy.ndarray
    :param field: the input radial profile
    :type field: numpy.ndarray
    :param ri: the inner core radius
    :type ri: float
    :param ro: the outer core radius
    :type ro: float
    :param lambdai: thickness of the inner boundary layer
    :type lambdai: float
    :param lambdao: thickness of the outer boundary layer
    :type lambdao: float
    :param normed: when set to True, the outputs are normalised by the volumes
                   of the boundary layers and the fluid bulk, respectively. In
                   that case, the outputs are volume-averaged quantities.
    :type normed: bool
    :returns: two floats that contains the boundary layer and the bulk
              integrations (integBc, integBulk)
    :rtype: list
    """
    # Dissipation in the boundary layers
    field2 = field.copy()
    mask = (rad<=ro-lambdao) * (rad>=ri+lambdai)
    field2[mask] = 0.
    integBc = intcheb(field2, len(field2)-1, ri, ro)

    if normed:
        volBc1 = 4./3.*np.pi*(ro**3-(ro-lambdao)**3)
        volBc2 = 4./3.*np.pi*((ri+lambdai)**3-(ri)**3)
        volBc = volBc1+volBc2
        integBc = integBc/volBc

    # Dissipation in the bulk
    field2 = field.copy()
    mask = (rad>ro-lambdao)
    field2[mask] = 0.
    mask = (rad<ri+lambdai)
    field2[mask] = 0.
    integBulk = intcheb(field2, len(field2)-1, ri, ro)

    if normed:
        volBulk = 4./3.*np.pi*((ro-lambdao)**3-(ri+lambdai)**3)
        integBulk = integBulk/volBulk

    return integBc, integBulk

def integBotTop(rad, field, ri, ro, lambdai, lambdao, normed=False):
    """
    This function evaluates the radial integral of the input array field
    in the bottom and top boundary layers separately.

    :param rad: radius
    :type rad: numpy.ndarray
    :param field: the input radial profile
    :type field: numpy.ndarray
    :param ri: the inner core radius
    :type ri: float
    :param ro: the outer core radius
    :type ro: float
    :param lambdai: thickness of the inner boundary layer
    :type lambdai: float
    :param lambdao: thickness of the outer boundary layer
    :type lambdao: float
    :param normed: when set to True, the outputs are normalised by the volumes
                   of the boundary layers. In that case, the outputs are
                   volume-averaged quantities.
    :type normed: bool
    :returns: two floats that contains the bottom and top boundary layers
              integrations (integBot, integTop)
    :rtype: list
    """
    field2 = field.copy()
    mask = (rad<=ro-lambdao)
    field2[mask] = 0.
    integTop = intcheb(field2, len(field2)-1, ri, ro)
    field2 = field.copy()
    mask = (rad>=ri+lambdai)
    field2[mask] = 0.
    integBot = intcheb(field2, len(field2)-1, ri, ro)

    if normed:
        volBc1 = 4./3.*np.pi*(ro**3-(ro-lambdao)**3)
        volBc2 = 4./3.*np.pi*((ri+lambdai)**3-(ri)**3)
        integTop /= volBc1
        integBot /= volBc2

    return integBot, integTop

def getMaxima(field):
    """
    This function determines the local maxima of the input array field

    :param field: the input array
    :type field: numpy.ndarray
    :returns:  a list containing the indices of the local maxima
    :rtype: list
    """
    maxS = []
    for k in range(len(field)):
        if k > 3 and k < len(field)-3:
            if field[k] > field[k-1] and field[k] > field[k+1]:
                maxS.append(k)
    return maxS

class BLayers(MagicSetup):
    """
    This class allows to determine the viscous and thermal boundary layers
    using several classical methods (slope method, peak values, dissipation
    rates, etc.). It uses the following files:

       * Kinetic energy: :ref:`eKinR.TAG <secEkinRFile>`
       * Power budget: :ref:`powerR.TAG <secPowerRfile>`
       * Radial profiles used for boundary layers: :ref:`bLayersR.TAG <secBLayersRfile>`

    This function can thus **only** be used when both
    :ref:`powerR.TAG <secPowerRfile>` and :ref:`bLayersR.TAG <secBLayersRfile>`
    exist in the working directory.

    .. warning:: This function works well as long as rigid boundaries and
                 fixed temperature boundary conditions are employed. Other
                 combination of boundary conditions (fixed fluxes and/or
                 stress-free) might give wrong results, since boundary layers
                 become awkward to define in that case.

    Since this function is supposed to use time-averaged quantities, the usual
    procedure is first to define the initial averaging time using
    :py:class:`AvgField <magic.AvgField>`: (this needs to be done only once)

    >>> a = AvgField(tstart=2.58)

    Once the ``tInitAvg`` file exists, the boundary layer calculation can be
    done:

    >>> bl = BLayers(iplot=True)
    >>> # print the formatted output
    >>> print(bl)
    """

    def __init__(self, iplot=False, quiet=False):
        """
        :param iplot: display the result when set to True (default False)
        :type iplot: bool
        :param quiet: less verbose when set to True (default is False)
        :type quiet: bool
        """
        if os.path.exists('tInitAvg'):
            file = open('tInitAvg', 'r')
            tstart = float(file.readline())
            file.close()
            logFiles = scanDir('log.*')
            tags = []
            for lg in logFiles:
                nml = MagicSetup(quiet=True, nml=lg)
                if nml.start_time >  tstart:
                    if os.path.exists('bLayersR.{}'.format(nml.tag)):
                        tags.append(nml.tag)
            if len(tags) > 0:
                print(tags)
            else:
                tags = None
            MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])

            a = AvgField()
            self.nuss = a.nuss
            self.reynolds = a.reynolds
            e2fluct = a.ekin_pol_avg+a.ekin_tor_avg-a.ekin_pola_avg-a.ekin_tora_avg
        else:
            logFiles = scanDir('log.*')
            MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])
            tags = None
            self.nuss = 1.
            self.reynolds = 1.
            e2fluct = 1.
        par = MagicRadial(field='bLayersR', iplot=False, tags=tags)
        self.varS = abs(par.entropy_SD)
        self.ss = par.entropy

        if os.path.exists('tInitAvg'):
            logFiles = scanDir('log.*', tfix=1409827718.0)
            # Workaround for code mistake before this time
            tfix = 1409827718.0
            tagsFix = []
            for lg in logFiles:
                nml = MagicSetup(quiet=True, nml=lg)
                if nml.start_time >  tstart:
                    if os.path.exists('bLayersR.{}'.format(nml.tag)):
                        tagsFix.append(nml.tag)
            if len(tagsFix) > 0:
                print('Fix temp. tags', tagsFix)
                parFix = MagicRadial(field='bLayersR', iplot=False, tags=tagsFix)
                self.varS = abs(parFix.entropy_SD)
                self.ss = parFix.entropy

            self.tags = tagsFix
        self.uh = par.uh
        self.duh = par.duhdr
        self.rad = par.radius
        self.ro = self.rad[0]
        self.ri = self.rad[-1]

        vol_oc = 4./3.* np.pi * (self.ro**3-self.ri**3)
        self.rey_fluct = np.sqrt(2.*e2fluct/vol_oc)

        self.reh = 4.*np.pi*intcheb(self.rad**2*self.uh, len(self.rad)-1,
                        self.ri, self.ro)/(4./3.*np.pi*(self.ro**3-self.ri**3))

        # Thermal dissipation boundary layer
        if hasattr(par, 'dissS'):
            self.dissS = par.dissS
            self.epsT = -4.*np.pi*intcheb(self.rad**2*self.dissS, len(self.rad)-1,
                                         self.ro, self.ri)
            self.epsTR = 4.*np.pi*self.rad**2*self.dissS
            ind = getMaxima(-abs(self.epsTR-self.epsT))

            try:
                self.dissTopS = self.ro-self.rad[ind[0]]
                self.dissBotS = self.rad[ind[-1]]-self.ri
                self.dissEpsTbl, self.dissEpsTbulk = integBulkBc(self.rad, self.epsTR,
                             self.ri, self.ro, self.dissBotS, self.dissTopS)
            except IndexError:
                self.dissTopS = self.ro
                self.dissBotS = self.ri
                self.dissEpsTbl, self.dissEpsTbulk = 0., 0.


            print('thDiss bl, bulk',  self.dissEpsTbl/self.epsT,
                  self.dissEpsTbulk/self.epsT)
        # First way of defining the thermal boundary layers: with var(S)
        #rThLayer = getMaxima(self.rad, self.varS)
        ind = argrelextrema(self.varS, np.greater)[0]
        if len(ind) != 0:
            self.bcTopVarS = self.ro-self.rad[ind[0]]
            self.bcBotVarS = self.rad[ind[-1]]-self.ri
        else:
            self.bcTopVarS = 1.
            self.bcBotVarS = 1.
        if hasattr(self, 'epsT'):
            self.varSEpsTbl, self.varSEpsTbulk = integBulkBc(self.rad, self.epsTR,
                         self.ri, self.ro, self.bcBotVarS, self.bcTopVarS)
            print('var(S) bl, bulk', self.varSEpsTbl/self.epsT, self.varSEpsTbulk/self.epsT)

        # Second way of defining the thermal boundary layers: intersection of the slopes
        d1 = matder(len(self.rad)-1, self.ro, self.ri)
        self.ttm = 3.*intcheb(self.ss*self.rad**2, len(self.rad)-1, self.ri, self.ro) \
                   /(self.ro**3-self.ri**3)
        dsdr = np.dot(d1, self.ss)
        self.beta = dsdr[len(dsdr)//2]
        print('beta={:.2f}'.format(self.beta))
        self.slopeTop = dsdr[2]*(self.rad-self.ro)+self.ss[0]
        self.slopeBot = dsdr[-1]*(self.rad-self.ri)+self.ss[-1]

        self.dtdrm = dsdr[len(self.ss)//2]
        tmid = self.ss[len(self.ss)//2]
        self.slopeMid = self.dtdrm*(self.rad-self.rad[len(self.rad)//2])+tmid

        self.bcTopSlope = (tmid-self.ss[0])/(self.dtdrm-dsdr[2])
        self.bcBotSlope = -(tmid-self.ss[-1])/(self.dtdrm-dsdr[-1])

        # 2nd round with a more accurate slope
        bSlope = dsdr[self.rad <= self.ri+self.bcBotSlope/4.].mean()
        tSlope = dsdr[self.rad >= self.ro-self.bcTopSlope/4.].mean()
        self.slopeBot = bSlope*(self.rad-self.ri)+self.ss[-1]
        self.slopeTop = tSlope*(self.rad-self.ro)+self.ss[0]
        #self.bcTopSlope = -(self.ttm-self.ss[0])/tSlope
        self.bcTopSlope = -(tmid-self.dtdrm*self.rad[len(self.rad)//2] - self.ss[0] \
                          + tSlope*self.ro)/(self.dtdrm-tSlope)
        self.bcBotSlope = -(tmid-self.dtdrm*self.rad[len(self.rad)//2] - self.ss[-1] \
                          + bSlope*self.ri)/(self.dtdrm-bSlope)
        self.dto = tSlope*(self.bcTopSlope-self.ro)+self.ss[0]
        self.dti = bSlope*(self.bcBotSlope-self.ri)+self.ss[-1]
        self.dto = self.dto-self.ss[0]
        self.dti = self.ss[-1]-self.dti

        self.bcTopSlope = self.ro - self.bcTopSlope
        self.bcBotSlope = self.bcBotSlope - self.ri

        if hasattr(self, 'epsT'):
            self.slopeEpsTbl, self.slopeEpsTbulk = integBulkBc(self.rad, self.epsTR,
                         self.ri, self.ro, self.bcBotSlope, self.bcTopSlope)

            print('slopes bl, bulk', self.slopeEpsTbl/self.epsT,
                  self.slopeEpsTbulk/self.epsT)

        pow = MagicRadial(field='powerR', iplot=False, tags=tags)
        self.vi = pow.viscDiss
        self.buo = pow.buoPower

        self.epsV = -intcheb(self.vi, len(self.rad)-1, self.ro, self.ri)
        ind = getMaxima(-abs(self.vi-self.epsV))
        if len(ind) > 2:
            for i in ind:
                if self.vi[i-1]-self.epsV > 0 and self.vi[i+1]-self.epsV < 0:
                    self.dissTopV = self.ro-self.rad[i]
                elif self.vi[i-1]-self.epsV < 0 and self.vi[i+1]-self.epsV > 0:
                    self.dissBotV = self.rad[i]-self.ri
        else:
            self.dissTopV = self.ro-self.rad[ind[0]]
            self.dissBotV = self.rad[ind[-1]]-self.ri
        try:
            self.dissEpsVbl, self.dissEpsVbulk = integBulkBc(self.rad, self.vi,
                             self.ri, self.ro, self.dissBotV, self.dissTopV)
        except AttributeError:
            self.dissTopV = 0.
            self.dissBotV = 0.
            self.dissEpsVbl = 0.
            self.dissEpsVbulk = 0.

        print('visc Diss bl, bulk', self.dissEpsVbl/self.epsV,
              self.dissEpsVbulk/self.epsV)

        # First way of defining the viscous boundary layers: with duhdr
        #rViscousLayer = getMaxima(self.rad, self.duh)
        if self.kbotv == 1 and self.ktopv == 1:
            ind = argrelextrema(self.duh, np.greater)[0]
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
            self.slopeTopU = 0.
            self.slopeBotU = 0.
            self.uhTopSlope = 0.
            self.uhBotSlope = 0.
            self.slopeEpsUbl = 0.
            self.slopeEpsUbulk = 0.
            self.uhBot = 0.
            self.uhTop = 0.
        else:
            ind = argrelextrema(self.uh, np.greater)[0]
            if len(ind) == 1:
                ind = argrelextrema(self.uh, np.greater_equal)[0]
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

            self.uhTop = self.uh[self.rad==self.ro-self.bcTopduh][0]
            self.uhBot = self.uh[self.rad==self.ri+self.bcBotduh][0]

            self.bcBotduh, self.bcTopduh, self.uhBot, self.uhTop =      \
                        getAccuratePeaks(self.rad, self.uh, self.uhTop, \
                                         self.uhBot, self.ri, self.ro)

            duhdr = np.dot(d1, self.uh)

            #1st round
            mask = (self.rad>=self.ro-self.bcTopduh/4)*(self.rad<self.ro)
            slopeT = duhdr[mask].mean()
            mask = (self.rad<=self.ri+self.bcBotduh/4)*(self.rad>self.ri)
            slopeB = duhdr[mask].mean()
            self.slopeTopU = slopeT*(self.rad-self.ro)+self.uh[0]
            self.slopeBotU = slopeB*(self.rad-self.ri)+self.uh[-1]
            self.uhTopSlope = -self.uhTop/slopeT
            self.uhBotSlope = self.uhBot/slopeB

            #2nd round
            mask = (self.rad>=self.ro-self.uhTopSlope/4.)*(self.rad<self.ro)
            slopeT = duhdr[mask].mean()
            mask = (self.rad<=self.ri+self.uhBotSlope/4)*(self.rad>self.ri)
            slopeB = duhdr[mask].mean()
            self.uhTopSlope = -self.uhTop/slopeT
            self.uhBotSlope = self.uhBot/slopeB

            self.slopeEpsUbl, self.slopeEpsUbulk = integBulkBc(self.rad, self.vi,
                         self.ri, self.ro, self.uhBotSlope, self.uhTopSlope)

        self.uhEpsVbl, self.uhEpsVbulk = integBulkBc(self.rad, self.vi,
                         self.ri, self.ro, self.bcBotduh, self.bcTopduh)
        print('uh bl, bulk', self.uhEpsVbl/self.epsV, self.uhEpsVbulk/self.epsV)

        # Convective Rol in the thermal boundary Layer
        par = MagicRadial(field='parR', iplot=False, tags=tags)
        kin = MagicRadial(field='eKinR', iplot=False, tags=tags)
        ekinNas = kin.ekin_pol+kin.ekin_tor-kin.ekin_pol_axi-kin.ekin_tor_axi
        ReR = np.sqrt(2.*abs(ekinNas)/par.radius**2/(4.*np.pi))
        RolC = ReR*par.ek/par.dlVc

        self.dl = par.dlVc
        y = RolC[par.radius >= self.ro-self.bcTopSlope]
        x = par.radius[par.radius >= self.ro-self.bcTopSlope]
        try:
            self.rolTop = simps(3.*y*x**2, x)/(self.ro**3-(self.ro-self.bcTopSlope)**3)
        except IndexError:
            self.rolTop = 0.

        self.rolbl, self.rolbulk = integBulkBc(self.rad, 4.*np.pi*RolC*self.rad**2,
                                     self.ri, self.ro, self.bcBotSlope, self.bcTopSlope,
                                     normed=True)

        self.rebl, self.rebulk = integBulkBc(self.rad, 4.*np.pi*ReR*self.rad**2,
                                     self.ri, self.ro, self.bcBotSlope, self.bcTopSlope,
                                     normed=True)

        self.lengthbl, self.lengthbulk = integBulkBc(self.rad, self.dl*4.*np.pi*self.rad**2,
                                     self.ri, self.ro, self.bcBotSlope, self.bcTopSlope,
                                     normed=True)

        self.rehbl, self.rehbulk = integBulkBc(self.rad, self.uh*4.*np.pi*self.rad**2,
                                     self.ri, self.ro, self.bcBotduh, self.bcTopduh,
                                     normed=True)

        y = RolC[par.radius <= self.ri+self.bcBotSlope]
        x = par.radius[par.radius <= self.ri+self.bcBotSlope]
        self.rolBot = simps(3.*y*x**2, x)/((self.ri+self.bcBotSlope)**3-self.ri**3)
        print('reynols bc, reynolds bulk', self.rebl, self.rebulk)
        print('reh bc, reh bulk', self.rehbl, self.rehbulk)
        print('rolbc, rolbulk, roltop, rolbot', self.rolbl, self.rolbulk,
              self.rolBot, self.rolTop)

        par.dlVc[0] = 0.
        par.dlVc[-1] = 0.
        self.lBot, self.lTop = integBotTop(self.rad, 4.*np.pi*self.rad**2*par.dlVc,
                         self.ri, self.ro, self.bcBotSlope, self.bcTopSlope, normed=True)

        uhbm, utbm = integBotTop(self.rad, 4.*np.pi*self.uh,
                         self.ri, self.ro, self.bcBotSlope, self.bcTopSlope, normed=True)

        # Convective Rol in the thermal boundary Layer
        if len(scanDir('perpParR.*')) != 0:
            tags = []
            for lg in logFiles:
                nml = MagicSetup(quiet=True, nml=lg)
                if nml.start_time >  tstart:
                    if os.path.exists('perpParR.{}'.format(nml.tag)):
                        tags.append(nml.tag)
            perpPar = MagicRadial(field='perpParR', iplot=False, tags=tags)
            eperpNas = perpPar.Eperp-perpPar.Eperp_axi
            eparNas = perpPar.Epar-perpPar.Epar_axi
            RePerpNas = np.sqrt(2.*abs(eperpNas))
            ReParNas = np.sqrt(2.*abs(eparNas))
            RePerp = np.sqrt(2.*abs(perpPar.Eperp))
            RePar = np.sqrt(2.*abs(perpPar.Epar))

            self.reperpbl, self.reperpbulk = integBulkBc(self.rad,
                                             4.*np.pi*RePerp*self.rad**2,
                                             self.ri, self.ro, self.bcBotSlope,
                                             self.bcTopSlope, normed=True)
            self.reparbl, self.reparbulk = integBulkBc(self.rad,
                                           4.*np.pi*RePar*self.rad**2,
                                           self.ri, self.ro, self.bcBotSlope,
                                           self.bcTopSlope, normed=True)
            self.reperpnasbl, self.reperpnasbulk = integBulkBc(self.rad,
                                                   4.*np.pi*RePerpNas*self.rad**2,
                                                   self.ri, self.ro,
                                                   self.bcBotSlope,
                                                   self.bcTopSlope, normed=True)
            self.reparnasbl, self.reparnasbulk = integBulkBc(self.rad,
                                                 4.*np.pi*ReParNas*self.rad**2,
                                                 self.ri, self.ro,
                                                 self.bcBotSlope, self.bcTopSlope,
                                                 normed=True)
        else:
            self.reperpbl = 0.
            self.reperpbulk = 0.
            self.reparbl = 0.
            self.reparbulk = 0.
            self.reperpnasbl = 0.
            self.reperpnasbulk = 0.
            self.reparnasbl = 0.
            self.reparnasbulk = 0.

        if iplot:
            self.plot()

        if not quiet:
            print(self)

    def plot(self):
        """
        Plotting function
        """
        #plt.rcdefaults()
        fig = plt.figure()
        ax = fig.add_subplot(211)
        ax.plot(self.rad, self.ss)
        ax.axhline(self.ttm, color='gray', linestyle='-')
        ax.plot(self.rad, self.slopeTop, 'k--')
        ax.plot(self.rad, self.slopeBot, 'k--')
        ax.plot(self.rad, self.slopeMid, 'k--')
        ax.axvline(self.ri+self.bcBotSlope, color='r')
        ax.axvline(self.ro-self.bcTopSlope, color='r')
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
            ax.set_ylim(0., 1.1*self.uh.max())
            ax.plot(self.rad, self.uh)
            ax.plot(self.rad, self.slopeTopU, 'k--')
            ax.plot(self.rad, self.slopeBotU, 'k--')
            mask = (np.abs(self.rad-self.ri-self.bcBotduh)==np.abs(self.rad-self.ri-self.bcBotduh).min())
            ax.axhline(self.uh[mask], color='k', linestyle='--', xmin=0., xmax=self.bcBotduh)
            mask = (np.abs(self.rad-self.ro+self.bcTopduh)==np.abs(self.rad-self.ro+self.bcTopduh).min())
            ax.axhline(self.uh[mask], color='k', linestyle='--', xmin=1.-self.bcTopduh, xmax=1.)
            if labTex:
                ax.set_ylabel(r'$u_h$')
            else:
                ax.set_ylabel('uh')
        ax.axvline(self.ro-self.bcTopduh, color='k', linestyle='--')
        ax.axvline(self.ri+self.bcBotduh, color='k', linestyle='--')
        ax.plot(self.rad, self.vi/self.vi.max())
        ax.set_xlim(self.rad[-1], self.rad[0])
        ax.set_xlabel('Radius')

        fig = plt.figure()
        ax = fig.add_subplot(211)
        ax.semilogy(self.rad, self.vi)
        ax.axhline(self.epsV, color='k', linestyle='--')
        ax.axvline(self.ro-self.dissTopV, color='k', linestyle='--')
        ax.axvline(self.ri+self.dissBotV, color='k', linestyle='--')
        ax.axvline(self.ro-self.uhTopSlope, color='g', linestyle='-', lw=1.5)
        ax.axvline(self.ri+self.uhBotSlope, color='g', linestyle='-', lw=1.5)
        ax.set_xlim(self.rad[-1], self.rad[0])
        ax.set_xlabel('Radius')
        ax.set_ylabel('Viscous dissipation')
        if hasattr(self, 'dissS'):
            ax = fig.add_subplot(212)
            ax.semilogy(self.rad, self.epsTR)
            ax.axhline(self.epsT, color='k', linestyle=':')
            ax.axvline(self.ro-self.dissTopS, color='k', linestyle='--')
            ax.axvline(self.ri+self.dissBotS, color='k', linestyle='--')
            ax.axvline(self.ro-self.bcTopSlope, color='g', linestyle='-', lw=1.5)
            ax.axvline(self.ri+self.bcBotSlope, color='g', linestyle='-', lw=1.5)
            ax.set_xlim(self.rad[-1], self.rad[0])
            ax.set_xlabel('Radius')
            ax.set_ylabel('Thermal Dissipation')

    def __str__(self):
        """
        Formatted output
        """
        if self.ek == -1:
            ek = 0. # to avoid the -1 for the non-rotating cases
        else:
            ek = self.ek
        if self.mode == 0:
            st ='{:9.3e}{:9.2e:}{:9.2e}{:9.2e}{:5.2f}'.format(self.ra, ek, self.pr,
                                                              self.prmag, self.strat)
        else:
            st = '{:.3e}{:12.5e}{:5.2f}{:6.2f}{:6.2f}'.format(self.ra, ek, self.strat,
                                                              self.pr, self.radratio)

        st += '{:12.5e}{:12.5e}{:12.5e}'.format(self.nuss, self.reynolds, self.rey_fluct)
        st += '{:12.5e}'.format(self.epsT)
        st += '{:12.5e}{:12.5e}{:5.2f}{:5.2f}'.format(self.bcBotSlope, self.bcTopSlope,
                          self.slopeEpsTbl/self.epsT, self.slopeEpsTbulk/self.epsT)
        st += '{:12.5e}{:12.5e}{:5.2f}{:5.2f}'.format(self.bcBotVarS, self.bcTopVarS,
                            self.varSEpsTbl/self.epsT, self.varSEpsTbulk/self.epsT)
        st += '{:12.5e}{:12.5e}{:5.2f}{:5.2f}'.format(self.dissBotS, self.dissTopS,
                            self.dissEpsTbl/self.epsT, self.dissEpsTbulk/self.epsT)

        st += '{:12.5e}'.format(self.epsV)
        st += '{:12.5e}{:12.5e}{:5.2f}{:5.2f}'.format(self.bcBotduh, self.bcTopduh,
                            self.uhEpsVbl/self.epsV, self.uhEpsVbulk/self.epsV)
        st += '{:12.5e}{:12.5e}{:5.2f}{:5.2f}'.format(self.uhBotSlope, self.uhTopSlope,
                          self.slopeEpsUbl/self.epsV, self.slopeEpsUbulk/self.epsV)
        st += '{:12.5e}{:12.5e}{:5.2f}{:5.2f}'.format(self.dissBotV, self.dissTopV,
                            self.dissEpsVbl/self.epsV, self.dissEpsVbulk/self.epsV)
        st += ' {:12.5e}'.format(self.beta)
        st += '{:12.5e}{:12.5e}'.format(abs(self.rolbl), abs(self.rolbulk))
        st += '{:12.5e}{:12.5e}'.format(self.rebl, self.rebulk)
        st += '{:12.5e}{:12.5e}'.format(self.rehbl, self.rehbulk)
        st += '{:12.5e}{:12.5e}'.format(self.lengthbl, self.lengthbulk)
        st += '{:12.5e}{:12.5e}'.format(self.ss[len(self.ss)//2]-self.ss[0],
                                        self.ttm-self.ss[0])
        st += '{:12.5e}{:12.5e}'.format(self.dti, self.dto)
        st += '{:12.5e}{:12.5e}{:12.5e}'.format(self.reh, self.uhBot, self.uhTop)
        st += '{:12.5e}{:12.5e}'.format(self.lBot, self.lTop)

        st  += '{:12.5e}{:12.5e}{:12.5e}{:12.5e}'.format(self.reperpbl, self.reperpbulk,
            self.reparbl, self.reparbulk)
        st  += '{:12.5e}{:12.5e}{:12.5e}{:12.5e}'.format(self.reperpnasbl,
            self.reperpnasbulk, self.reparnasbl, self.reparnasbulk)

        return st


if __name__ == '__main__':

    b = BLayers(iplot=True)
    print(b)
    plt.show()
