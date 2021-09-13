# -*- coding: utf-8 -*-
from magic import *
from scipy.integrate import trapz
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import numpy as np
import os


class CompSims:
    """
    This class allows to compare an analyse several DNS simultaneously. It is possible
    to compare time-series or :ref:`graphic files <secGraphFile>`. To set it up, you
    first need to create a file that contains the list of directories you want to analyse:

    .. code-block:: bash

       $ cat inputList
       E3e4Eps5e3Q05
       E3e4Eps2e3Q07
       E3e4Eps2e3Q08
       E3e4Eps2e3Q09

    This list thus contains four directories (one run per directory) that can be further
    analysed:

    >>> # Display the time-series of kinetic energy on 2 columns
    >>> CompSims(file='inputList', field='ts', ncol=2)
    >>> # Display the equatorial cuts of v_r
    >>> CompSims(file='inputList', field='vr', type='equat', levels=65, cm='seismic')
    >>> # Display the radial cuts of B_r at r=0.8 r_o
    >>> CompSims(file='inputList', field='br', type='surf', r=0.8)
    >>> # Display the average zonal flow
    >>> CompSims(file='inputList', field='vp', type='avg')
    """

    def __init__(self, file='liste', field='ts', ncol=4, cm='RdYlBu_r', dpi=96,
                 normed=True, levels=16, type=None,
                 r=0.9, bw=False, ave=False, cut=1):
        """
        :param file: the input file that contains the list of directories that one
                     wants to analyse
        :type file: str
        :param field: name of the input field. Possible options are:
                      'ts': displaye the time-series of kinetic energy;
                      'e_mag': display the time-series of magnetic energy;
                      'flux': display the time-series of the Nusselt numbers;
                      'zonal': display the surface zonal flow;
                      'Anything else': try to interpret the field
        :type field: str
        :param type: nature of the plot. Possible values are:
                     'avg' or 'slice': phi-average or phi-slice;
                     'equat': equatorial cut;
                     'surf': radial cut;
                     'ts*: time series
        :type type: str
        :param ncol: number of columns of the figure
        :type ncol: int
        :param ave: when set to True, it tries to read a time-averaged graphic file
        :type ave: bool
        :param r: the radius at which you want to display the input
                  data (in normalised units with the radius of the outer boundary)
        :type r: float
        :param levels: the number of levels in the contour
        :type levels: int
        :param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
        :type cm: str
        :param normed: when set to True, the colormap is centered around zero.
                       Default is True, except for entropy/temperature plots.
        :type normed: bool
        :param dpi: dot per inch when saving PNGs
        :type dpi: int
        :param bw: when set to True, display grey-scaled contour levels
        :type bw: bool
        :param cut: adjust the contour extrema to max(abs(data))*cut
        :type cut: float
        """
        self.dataliste = []
        self.workdir = os.getcwd()
        self.field = field
        self.cm = cm
        self.normed = normed
        self.cut = cut
        self.levels = levels
        self.r = r
        self.bw = bw # for black and white outputs
        self.ave = ave # for G_ave.TAG files
        f = open(file, 'r')
        for line in f.readlines():
            self.dataliste.append(line.strip())
        f.close()
        self.ncol = ncol
        self.nplot = len(self.dataliste)
        if (self.nplot % self.ncol != 0):
                self.nrow = self.nplot/self.ncol + 1
        else:
            self.nrow = self.nplot/self.ncol

        plt.ioff()
        if type == 'avg' or type == 'slice':
            plt.figure(figsize=(self.ncol*1.5, self.nrow*3), dpi=dpi)
        elif type == 'equat':
            plt.figure(figsize=(self.ncol*2.5, self.nrow*2.5), dpi=dpi)
        elif type == 'surf':
            plt.figure(figsize=(self.ncol*3, self.nrow*1.7), dpi=dpi)
        else:
            plt.figure(figsize=(self.ncol*3, self.nrow*3), dpi=dpi)
        if self.nrow == 1:
            if type == 'surf':
                plt.subplots_adjust(left=0.05, right=0.98, top=0.85, bottom=0.05)
            else:
                plt.subplots_adjust(left=0.05, right=0.98, top=0.92, bottom=0.05)
        else:
            if type == 'surf':
                plt.subplots_adjust(left=0.05, right=0.98, top=0.92, bottom=0.05)
            else:
                plt.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.05)

        if self.field == 'ts':
            self.plotTs()
        elif self.field == 'e_mag':
            self.plotEmag()
        elif self.field == 'flux':
            self.plotFlux()
        elif self.field == 'zonal':
            self.plotZonal()
        else:
            if type == 'avg':
                self.plotAvg()
            elif type == 'slice':
                self.plotSlice()
            elif type == 'equat':
                self.plotEquat()
            elif type == 'surf':
                self.plotSurf()
        plt.show()
        plt.ion()

    def plotTs(self):
        """
        Plot time-series of the kinetic energy
        """
        iplot = 1
        #myyfmt = ScalarFormatter(useOffset=True)
        #myyfmt.set_powerlimits((1,1))
        for datadir in self.dataliste:
            os.chdir(self.workdir + '/' + datadir)
            print(datadir)
            ts = MagicTs(field='e_kin', iplot=False, all=True)

            ax = plt.subplot(self.nrow, self.ncol, iplot)
            ax.semilogy(ts.time, ts.ekin_pol, 'b-')
            ax.semilogy(ts.time, ts.ekin_tor, 'r-')
            ax.semilogy(ts.time, ts.ekin_pol_axi, 'b--')
            ax.semilogy(ts.time, ts.ekin_tor_axi, 'r--')
            ax.semilogy(ts.time, ts.ekin_tot, 'k-')

            #ax.yaxis.set_major_formatter(myyfmt)
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(10)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(10)
            ax.set_title('Ra = {:.1e}'.format(ts.ra), fontsize=10)
            ax.set_xlim((ts.time.min(), ts.time.max()))
            iplot += 1
        os.chdir(self.workdir)

    def plotEmag(self):
        """
        Plot time-series of the magnetic energy
        """
        iplot = 1
        #myyfmt = ScalarFormatter(useOffset=True)
        #myyfmt.set_powerlimits((1,1))
        for datadir in self.dataliste:
            os.chdir(self.workdir + '/' + datadir)
            print(datadir)
            ts = MagicTs(field='e_mag_oc', iplot=False, all=True)

            ax = plt.subplot(self.nrow, self.ncol, iplot)
            ax.semilogy(ts.time, ts.emagoc_pol, 'b-')
            ax.semilogy(ts.time, ts.emagoc_tor, 'r-')
            ax.semilogy(ts.time, ts.emagoc_pol_axi, 'b--')
            ax.semilogy(ts.time, ts.emagoc_tor_axi, 'r--')
            ax.semilogy(ts.time, ts.emag_tot, 'k-')

            #ax.yaxis.set_major_formatter(myyfmt)
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(10)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(10)
            ax.set_title('Ra = {:.1e}, Pm= {:.1f}'.format(ts.ra, ts.prmag), fontsize=10)
            ax.set_xlim((ts.time.min(), ts.time.max()))
            iplot += 1
        os.chdir(self.workdir)

    def plotFlux(self):
        """
        Plot time-series of the top and bottom Nusselt numbers
        """
        iplot = 1
        for datadir in self.dataliste:
            os.chdir(self.workdir + '/' + datadir)
            print(datadir)
            ts = MagicTs(field='misc', iplot=False, all=True)

            ax = plt.subplot(self.nrow, self.ncol, iplot)
            ax.plot(ts.time, ts.botnuss, 'b-')
            ax.plot(ts.time, ts.topnuss, 'g-')

            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(10)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(10)
            ax.set_title('Ra = {:.1e}'.format(ts.ra), fontsize=10)
            ax.set_xlim((ts.time.min(), ts.time.max()))
            iplot += 1
        os.chdir(self.workdir)

    def plotZonal(self):
        """
        Plot surface zonal flow profiles.
        """
        iplot = 1
        for datadir in self.dataliste:
            os.chdir(self.workdir + '/' + datadir)
            if self.ave:
                gr = MagicGraph(ivar=1, ave=True)
            else:
                gr = MagicGraph(ivar=1)
            ax = plt.subplot(self.nrow, self.ncol, iplot)
            vpm = gr.vphi.mean(axis=0)
            theta = np.linspace(-90., 90, gr.ntheta)

            ax.plot(vpm[:, 1], theta)

            roequat = vpm[gr.ntheta/2, 0]*gr.ek*(1.-gr.radratio)
            print('{:7.3e} {:7.3e}'.format(gr.ra, roequat))
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(10)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(10)
            ax.set_title('Ra = {:.1e}'.format(gr.ra), fontsize=10)
            ax.set_xlim(1.1*vpm[:,0].min(), 1.1*vpm[:,0].max())
            ax.set_ylim(theta.min(), theta.max())
            ax.axvline(0., color='k', linestyle='--')


            iplot += 1
        os.chdir(self.workdir)

    def plotSurf(self):
        """
        Plot radial cuts in (phi, theta)  planes using the Hammer projection.
        """
        cmap = plt.get_cmap(self.cm)

        iplot = 1
        for datadir in self.dataliste:
            os.chdir(self.workdir + '/' + datadir)
            print(datadir)
            try:
                if self.ave:
                    gr = MagicGraph(ivar=1, ave=True)
                else:
                    gr = MagicGraph(ivar=1)

                rad = self.r/(1-gr.radratio) # as we give a normalised radius
                ind = np.nonzero(np.where(abs(gr.radius-rad) \
                                == min(abs(gr.radius-rad)), 1, 0))
                indPlot = ind[0][0]

                if self.field in ('Vs', 'vs'):
                    vr = gr.vr
                    vt = gr.vtheta
                    thlin = np.linspace(0., np.pi, gr.ntheta)
                    th3D = np.zeros_like(vr)
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    data = vr * np.sin(th3D) + vt * np.cos(th3D)
                    label = 'Vs'
                elif self.field in ('Vz', 'vz'):
                    vr = gr.vr
                    vt = gr.vtheta
                    thlin = np.linspace(0., np.pi, gr.ntheta)
                    th3D = np.zeros_like(vr)
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    data = vr * np.cos(th3D) - vt * np.sin(th3D)
                    label = 'Vz'
                else:
                    data, data_ic, label = selectField(gr, self.field)
            except AttributeError:
                continue

            data = symmetrize(data, gr.minc)

            phi2, th2 = np.mgrid[-np.pi:np.pi:gr.nphi*1j,
                                np.pi/2.:-np.pi/2.:gr.ntheta*1j]
            xx, yy = hammer2cart(th2, phi2)


            ax = plt.subplot(self.nrow, self.ncol, iplot, frameon=False)
            if self.cut != 1:
                self.normed = False
                vmin = - max(abs(data[..., indPlot].max()), abs(data[..., indPlot].min()))
                vmin = self.cut*vmin
                vmax = -vmin
                cs = np.linspace(vmin, vmax, self.levels)
                im = ax.contourf(xx, yy, data[..., indPlot], cs, extend='both',
                                  cmap=cmap, aa=True)
            else:
                cs = self.levels
                im = ax.contourf(xx, yy, data[..., indPlot], cs,
                                  cmap=cmap, aa=True)
            rad = gr.radius[indPlot] * (1. - gr.radratio)
            ax.set_title('{}, r/ro={:.3f}, Ra={:.1e}'.format(label, rad, gr.ra),
                    fontsize=10)
            ax.axis('off')


            if self.field not in ['entropy', 's', 'S'] and self.normed is True:
                im.set_clim(-max(abs(data[..., indPlot].max()),
                                 abs(data[..., indPlot].min())),
                             max(abs(data[..., indPlot].max()),
                                 abs(data[..., indPlot].min())))
            iplot += 1
        os.chdir(self.workdir)

    def plotEquat(self):
        """
        Plot equatorial cuts in (phi, r)  planes.
        """
        cmap = plt.get_cmap(self.cm)
        iplot = 1
        for datadir in self.dataliste:
            os.chdir(self.workdir + '/' + datadir)
            print(datadir)
            try:
                if self.ave:
                    gr = MagicGraph(ivar=1, ave=True)
                else:
                    gr = MagicGraph(ivar=1)
                if self.field in ('Vs', 'vs'):
                    vr = gr.vr
                    vt = gr.vtheta
                    thlin = np.linspace(0., np.pi, gr.ntheta)
                    th3D = np.zeros_like(vr)
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    data = vr * np.sin(th3D) + vt * np.cos(th3D)
                    label = 'Vs'
                elif self.field in ('Vz', 'vz'):
                    vr = gr.vr
                    vt = gr.vtheta
                    thlin = np.linspace(0., np.pi, gr.ntheta)
                    th3D = np.zeros_like(vr)
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    data = vr * np.cos(th3D) - vt * np.sin(th3D)
                    label = 'Vz'
                elif self.field in ('vortz'):
                    philoc = np.linspace(0., 2.*np.pi/gr.minc, gr.npI)
                    rrloc, pphiloc = np.meshgrid(gr.radius, philoc)
                    dr = rderavg(rrloc*gr.vphi[:,gr.ntheta/2,:], spectral=False,
                                 eta=gr.radratio, exclude=True)
                    equator = 1./rrloc*(dr - phideravg(gr.vr[:, gr.ntheta/2, :], gr.minc))
                    if labTex:
                        label = r'$\omega_z$'
                    else:
                        label = 'omega'
                else:
                    data, data_ic, label = selectField(gr, self.field)
            except AttributeError:
                continue


            label += ' Ra = {:.1e}'.format(gr.ra)

            if self.field not in ('vortz'):
                equator = data[:, gr.ntheta//2,:]

            equator = symmetrize(equator, gr.minc)

            phi = np.linspace(0., 2.*np.pi, gr.nphi)
            rr, pphi = np.meshgrid(gr.radius, phi)
            xx = rr * np.cos(pphi)
            yy = rr * np.sin(pphi)

            ax = plt.subplot(self.nrow, self.ncol, iplot, frameon=False)
            if self.bw:
                im = ax.contour(xx, yy, equator, self.levels, colors='k',
                                linewidths=0.5)
            else:
                if self.cut != 1:
                    self.normed = False
                    vmin = - max(abs(equator.max()), abs(equator.min()))
                    vmin = self.cut*vmin
                    vmax = -vmin
                    cs = np.linspace(vmin, vmax, self.levels)
                    im = ax.contourf(xx, yy, equator, cs, extend='both',
                                      cmap=cmap)
                else:
                    cs = self.levels
                    im = ax.contourf(xx, yy, equator, cs, cmap=cmap)
            ax.plot(gr.radius[0] * np.cos(phi), gr.radius[0]*np.sin(phi), 'k-')
            ax.plot(gr.radius[-1] * np.cos(phi), gr.radius[-1]*np.sin(phi), 'k-')

            # Variable conductivity
            if hasattr(gr, 'nVarCond'):
                if gr.nVarCond == 2:
                    radi = gr.con_radratio * gr.radius[0]
                    ax.plot(radi*np.cos(phi), radi*np.sin(phi), 'k--')

            #if hasattr(gr, 'cmbHflux'):
                #tit1 = r"${\cal Q}_{cmb} = {:.1f}$".format(gr.cmbHflux)
                #if gr.strat >= 1:
                    #tit1 = r"$N_\rho = {:.0f}$".format(gr.strat)
                #else:
                    #tit1 = r"$N_\rho = 10^{-2}$"
            tit1 = datadir

            ax.text(0.5, 0.5, tit1, fontsize=14,
                              horizontalalignment='center',
                              verticalalignment='center',
                              transform = ax.transAxes)

            #ax.set_title(label, fontsize=10)
            ax.axis('off')
            #fig.colorbar(im)

            if self.field not in ['entropy', 's', 'S'] and self.normed is True:
                im.set_clim(-max(abs(equator.max()), abs(equator.min())),
                             max(abs(equator.max()), abs(equator.min())))
            iplot += 1
        os.chdir(self.workdir)

    def plotAvg(self):
        """
        Plot azimutal averages in (theta, r) planes.
        """
        cmap = plt.get_cmap(self.cm)

        iplot = 1
        for datadir in self.dataliste:
            os.chdir(self.workdir + '/' + datadir)
            print(datadir)
            try:
                if self.ave:
                    gr = MagicGraph(ivar=1, ave=True)
                else:
                    gr = MagicGraph(ivar=1)
                if self.field in ('Vs', 'vs'):
                    vr = gr.vr
                    vt = gr.vtheta
                    thlin = np.linspace(0., np.pi, gr.ntheta)
                    th3D = np.zeros_like(vr)
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    data = vr * np.sin(th3D) + vt * np.cos(th3D)
                    label = 'Vs'
                elif self.field in ('Vz', 'vz'):
                    vr = gr.vr
                    vt = gr.vtheta
                    thlin = np.linspace(0., np.pi, gr.ntheta)
                    th3D = np.zeros_like(vr)
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    data = vr * np.cos(th3D) - vt * np.sin(th3D)
                    label = 'Vz'
                elif self.field in ('Cr', 'cr'):
                    vr = gr.vr
                    vt = gr.vtheta
                    vp = gr.vphi.copy()
                    vp = gr.vphi- gr.vphi.mean(axis=0) # convective vp
                    thlin = np.linspace(0., np.pi, gr.ntheta)
                    th3D = np.zeros_like(vr)
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    vs = vr * np.sin(th3D) + vt * np.cos(th3D)
                    data =  vs * vp
                    denom = np.sqrt(np.mean(vs**2, axis=0)* np.mean(vp**2, axis=0))
                    label = r'$\langle v_s v_\phi\rangle$'
                else:
                    data, data_ic, label = selectField(gr, self.field)
            except AttributeError:
                continue

            #label += ' Ra = {:.1e}'.format(gr.ra)
            label = 'Ra = {:.1e}'.format(gr.ra)

            if self.field not in ('Cr', 'cr', 'ra', 'ratio', 'Cz', 'cz'):
                phiavg = data.mean(axis=0)
            else:
                ro = gr.radius[0]
                ri = gr.radius[-1]
                fac = 2./(np.pi*(ro**2-ri**2))
                facOTC = ro**2.*(np.pi-2.*np.arcsin(gr.radratio))/2. \
                         -ri**2*np.sqrt(1.-gr.radratio**2)/gr.radratio
                facOTC = 1./facOTC
                facITC = ri**2*np.sqrt(1.-gr.radratio**2)/gr.radratio \
                         +(ro**2-ri**2)* np.arcsin(gr.radratio) \
                         -ri**2/2.*(np.pi - 2.*np.arcsin(gr.radratio))
                facITC = 1./facITC
                phiavg = data.mean(axis=0)

                TC = np.array([], dtype=data.dtype)
                outTC = np.array([], dtype=data.dtype)
                inTC = np.array([], dtype=data.dtype)
                integ = np.array([], dtype=data.dtype)
                for k, th in enumerate(gr.colatitude):
                    rr = gr.radius[::-1]
                    dat = phiavg[k, ::-1] * rr
                    corr = intcheb(dat, gr.nr-1, ri, ro)
                    TC = np.append(TC, corr)
                    if th >= np.arcsin(gr.radratio)  and \
                       th <= np.pi - np.arcsin(gr.radratio):
                        # Outside tangent cylinder
                        val = trapz(dat[rr >= ri/np.sin(th)],
                                    rr[rr >= ri/np.sin(th)])
                        outTC = np.append(outTC, val)
                        integ = np.append(integ, th)
                        # Inside tangent cylinder
                        val = trapz(dat[rr < ri/np.sin(th)],
                                    rr[rr < ri/np.sin(th)])
                        inTC = np.append(inTC, val)
                    else:
                        val= intcheb(dat, gr.nr-1, ri, ro)
                        inTC = np.append(inTC, val)

                mask = np.where(denom == 0, 1, 0)
                phiavg /= (denom+mask)


            th = np.linspace(0., np.pi, gr.ntheta)
            rr, tth = np.meshgrid(gr.radius, th)
            xx = rr * np.sin(tth)
            yy = rr * np.cos(tth)

            titmax = phiavg.max()
            titmin = phiavg.min()

            ax = plt.subplot(self.nrow, self.ncol, iplot, frameon=False)
            if self.bw:
                im = ax.contour(xx, yy, phiavg, self.levels, colors='k',
                                linewidths=0.5)
            else:
                if self.cut != 1:
                    self.normed = False
                    vmin = - max(abs(phiavg.max()), abs(phiavg.min()))
                    vmin = self.cut*vmin
                    vmax = -vmin
                    cs = np.linspace(vmin, vmax, self.levels)
                    im = ax.contourf(xx, yy, phiavg, cs, extend='both',
                                      cmap=cmap)
                else:
                    cs = self.levels
                    im = ax.contourf(xx, yy, phiavg, cs, cmap=cmap)
            ax.plot(gr.radius[0]*np.sin(th), gr.radius[0]*np.cos(th),
                   'k-')
            ax.plot(gr.radius[-1]*np.sin(th), gr.radius[-1]*np.cos(th),
                   'k-')

            ax.plot([0., 0], [gr.radius[-1], gr.radius[0]], 'k-')
            ax.plot([0., 0], [-gr.radius[-1], -gr.radius[0]], 'k-')

            # Variable conductivity
            if hasattr(gr, 'nVarCond'):
                if gr.nVarCond == 2:
                    radi = gr.con_radratio * gr.radius[0]
                    ax.plot(radi*np.sin(th), radi*np.cos(th), 'k--')

            ax.set_title(label, fontsize=12)
            ax.axis('off')
            #fig.colorbar(im)

            """
            if gr.strat >= 1:
                tit1 = r"$N_\rho = {:.0f}$".format(gr.strat)
            else:
                tit1 = r"$N_\rho = 10^{-2}$"
            """
            #plt.title(tit1, fontsize=12)

            """
            if int(titmin) == 0:
                tit1 = r'$+{}$'.format(titmax) +'\n'+r'${:.1f}$'.format(titmin)
            else:
                tit1 = r'$+{}$'.format(titmax) +'\n'+r'${}$'.format(titmin)
            ax.text(0., 0.5, tit1, fontsize=12,
                              horizontalalignment='left',
                              verticalalignment='center',
                              transform = ax.transAxes)
            tit2 = r'$+{:.1e}/-{:.1e}$'.format(titmax, titmin)
            """
            #ax.text(0.9, 0.05, tit2, fontsize=12,
                              #horizontalalignment='left',
                              #verticalalignment='center',
                              #transform = ax.transAxes)
            if self.field not in ['entropy', 's', 'S'] and self.normed is True:
                im.set_clim(-max(abs(phiavg.max()), abs(phiavg.min())),
                             max(abs(phiavg.max()), abs(phiavg.min())))

            iplot += 1
        os.chdir(self.workdir)

    def plotSlice(self):
        cmap = plt.get_cmap(self.cm)

        iplot = 1
        for datadir in self.dataliste:
            os.chdir(self.workdir + '/' + datadir)
            print(datadir)
            try:
                if self.ave:
                    gr = MagicGraph(ivar=1, ave=True)
                else:
                    gr = MagicGraph(ivar=1)
                if self.field in ('Vs', 'vs'):
                    vr = gr.vr
                    vt = gr.vtheta
                    thlin = np.linspace(0., np.pi, gr.ntheta)
                    th3D = np.zeros((gr.npI, gr.ntheta, gr.nr), dtype=vr.dtype)
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    data = vr * np.sin(th3D) + vt * np.cos(th3D)
                    label = 'Vs'
                elif self.field in ('Vz', 'vz'):
                    vr = gr.vr
                    vt = gr.vtheta
                    thlin = np.linspace(0., np.pi, gr.ntheta)
                    th3D = np.zeros((gr.npI, gr.ntheta, gr.nr), dtype=vr.dtype)
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    data = vr * np.cos(th3D) - vt * np.sin(th3D)
                    label = 'Vz'
                else:
                    data, data_ic, label = selectField(gr, self.field)
            except AttributeError:
                continue

            #label += ' Ra = {:.1e}'.format(gr.ra)
            label = 'Ra = {:.1e}'.format(gr.ra)

            phiavg = data[0, ...]

            th = np.linspace(0., np.pi, gr.ntheta)
            rr, tth = np.meshgrid(gr.radius, th)
            xx = rr * np.sin(tth)
            yy = rr * np.cos(tth)

            titmax = phiavg.max()
            titmin = phiavg.min()

            # liste2
            """
            if gr.strat == 4:
                vmax = 0.6*phiavg.max()
                vmin = -vmax
            elif gr.strat == 5:
                vmax = -0.5*phiavg.min()
                vmin = -vmax
            else:
                vmax = -0.8*phiavg.min()
                vmin = -vmax
            """

            ax = plt.subplot(self.nrow, self.ncol, iplot, frameon=False)
            if self.cut != 1:
                self.normed = False
                vmin = - max(abs(phiavg.max()), abs(phiavg.min()))
                vmin = self.cut*vmin
                vmax = -vmin
                cs = np.linspace(vmin, vmax, self.levels)
                im = ax.contourf(xx, yy, phiavg, cs, extend='both', cmap=cmap)
            else:
                cs = self.levels
                im = ax.contourf(xx, yy, phiavg, cs, cmap=cmap)
            ax.plot(gr.radius[0]*np.sin(th), gr.radius[0]*np.cos(th),
                   'k-')
            ax.plot(gr.radius[-1]*np.sin(th), gr.radius[-1]*np.cos(th),
                   'k-')

            ax.plot([0., 0], [gr.radius[-1], gr.radius[0]], 'k-')
            ax.plot([0., 0], [-gr.radius[-1], -gr.radius[0]], 'k-')

            # Variable conductivity
            if hasattr(gr, 'nVarCond'):
                if gr.nVarCond == 2:
                    radi = gr.con_radratio * gr.radius[0]
                    ax.plot(radi*np.sin(th), radi*np.cos(th), 'k--')

            ax.set_title(label, fontsize=12)
            ax.axis('off')
            #fig.colorbar(im)

            if gr.strat >= 1:
                tit1 = r"$N_\rho = {:.0f}$".format(gr.strat)
            else:
                tit1 = r"$N_\rho = 10^{-2}$"
            #plt.title(tit1, fontsize=12)

            if int(titmin) == 0:
                tit1 = r'$+{}$'.format(titmax) +'\n'+r'${:.1f}$'.format(titmin)
            else:
                tit1 = r'$+{}$'.format(titmax) +'\n'+r'${}$'.format(titmin)
            ax.text(0., 0.5, tit1, fontsize=12,
                              horizontalalignment='left',
                              verticalalignment='center',
                              transform = ax.transAxes)
            #ax.text(0.9, 0.05, tit2, fontsize=12,
                              #horizontalalignment='left',
                              #verticalalignment='center',
                              #transform = ax.transAxes)
            if self.field not in ['entropy', 's', 'S'] and self.normed is True:
                im.set_clim(-max(abs(phiavg.max()), abs(phiavg.min())),
                             max(abs(phiavg.max()), abs(phiavg.min())))

            iplot += 1
        os.chdir(self.workdir)

if __name__ == '__main__':
    CompSims(file='liste', field='ts', ncol=5, type=None)
