# -*- coding: utf-8 -*-
from magic import *
from scipy.integrate import trapz
from matplotlib.ticker import ScalarFormatter
import pylab as P
import numpy as N
import os


class CompSims:

    def __init__(self, file='liste', field='ts', ncol=4, cm='RdYlBu_r', dpi=96,
                 normed=True, levels=16, type=None,
                 r=0.9, bw=False, ave=False, cut=1):
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

        P.ioff()
        if type == 'avg' or type == 'slice':
            P.figure(figsize=(self.ncol*1.5, self.nrow*3), dpi=dpi)
        elif type == 'equat':
            P.figure(figsize=(self.ncol*2.5, self.nrow*2.5), dpi=dpi)
        elif type == 'surf':
            P.figure(figsize=(self.ncol*3, self.nrow*1.7), dpi=dpi)
        else:
            P.figure(figsize=(self.ncol*3, self.nrow*3), dpi=dpi)
        if self.nrow == 1:
            if type == 'surf':
                P.subplots_adjust(left=0.05, right=0.98, top=0.85, bottom=0.05)
            else:
                P.subplots_adjust(left=0.05, right=0.98, top=0.92, bottom=0.05)
        else:
            if type == 'surf':
                P.subplots_adjust(left=0.05, right=0.98, top=0.92, bottom=0.05)
            else:
                P.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.05)

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
        P.show()
        P.ion()

    def plotTs(self):
        iplot = 1
        #myyfmt = ScalarFormatter(useOffset=True)
        #myyfmt.set_powerlimits((1,1))
        for datadir in self.dataliste:
            os.chdir(self.workdir + '/' + datadir)
            print datadir
            ts = MagicTs(field='e_kin', iplot=False, all=True)

            ax = P.subplot(self.nrow, self.ncol, iplot)
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
            ax.set_title('Ra = %.1e' % ts.ra, fontsize=10)
            ax.set_xlim((ts.time.min(), ts.time.max()))
            iplot += 1
        os.chdir(self.workdir)

    def plotEmag(self):
        iplot = 1
        #myyfmt = ScalarFormatter(useOffset=True)
        #myyfmt.set_powerlimits((1,1))
        for datadir in self.dataliste:
            os.chdir(self.workdir + '/' + datadir)
            print datadir
            ts = MagicTs(field='e_mag_oc', iplot=False, all=True)

            ax = P.subplot(self.nrow, self.ncol, iplot)
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
            ax.set_title('Ra = %.1e, Pm= %.1f' % (ts.ra, ts.prmag), fontsize=10)
            ax.set_xlim((ts.time.min(), ts.time.max()))
            iplot += 1
        os.chdir(self.workdir)

    def plotFlux(self):
        iplot = 1
        for datadir in self.dataliste:
            os.chdir(self.workdir + '/' + datadir)
            print datadir
            ts = MagicTs(field='misc', iplot=False, all=True)

            ax = P.subplot(self.nrow, self.ncol, iplot)
            ax.plot(ts.time, ts.botnuss, 'b-')
            ax.plot(ts.time, ts.topnuss, 'g-')

            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(10)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(10)
            ax.set_title('Ra = %.1e' % ts.ra, fontsize=10)
            ax.set_xlim((ts.time.min(), ts.time.max()))
            iplot += 1
        os.chdir(self.workdir)

    def plotZonal(self):
        iplot = 1
        for datadir in self.dataliste:
            os.chdir(self.workdir + '/' + datadir)
            if self.ave:
                gr = MagicGraph(ivar=1, ave=True)
            else:
                gr = MagicGraph(ivar=1)
            ax = P.subplot(self.nrow, self.ncol, iplot)
            vpm = gr.vphi.mean(axis=0)
            theta = N.linspace(-90., 90, gr.ntheta)

            ax.plot(vpm[:, 1], theta)

            roequat = vpm[gr.ntheta/2, 0]*gr.ek*(1.-gr.radratio)
            print '%7.3e %7.3e' % (gr.ra, roequat)
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(10)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(10)
            ax.set_title('Ra = %.1e' % gr.ra, fontsize=10)
            ax.set_xlim(1.1*vpm[:,0].min(), 1.1*vpm[:,0].max())
            ax.set_ylim(theta.min(), theta.max())
            ax.axvline(0., color='k', linestyle='--')


            iplot += 1
        os.chdir(self.workdir)

    def plotSurf(self):
        cmap = P.get_cmap(self.cm)

        iplot = 1
        for datadir in self.dataliste:
            os.chdir(self.workdir + '/' + datadir)
            print datadir
            try:
                if self.ave:
                    gr = MagicGraph(ivar=1, ave=True)
                else:
                    gr = MagicGraph(ivar=1)

                rad = self.r/(1-gr.radratio) # as we give a normalised radius
                ind = N.nonzero(N.where(abs(gr.radius-rad) \
                                == min(abs(gr.radius-rad)), 1, 0))
                indPlot = ind[0][0]

                if self.field in ('Vs', 'vs'):
                    vr = gr.vr
                    vt = gr.vtheta
                    thlin = N.linspace(0., N.pi, gr.ntheta)
                    th3D = N.zeros((gr.nphi, gr.ntheta, gr.nr), dtype='f')
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    data = vr * N.sin(th3D) + vt * N.cos(th3D)
                    label = 'Vs'
                elif self.field in ('Vz', 'vz'):
                    vr = gr.vr
                    vt = gr.vtheta
                    thlin = N.linspace(0., N.pi, gr.ntheta)
                    th3D = N.zeros((gr.nphi, gr.ntheta, gr.nr), dtype='f')
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    data = vr * N.cos(th3D) - vt * N.sin(th3D)
                    label = 'Vz'
                else:
                    data, label = selectField(gr, self.field)
            except AttributeError:
                continue

            data = symmetrize(data, gr.minc)

            phi2, th2 = N.mgrid[-N.pi:N.pi:gr.nphi*1j,
                                N.pi/2.:-N.pi/2.:gr.ntheta*1j]
            xx, yy = hammer2cart(th2, phi2)


            ax = P.subplot(self.nrow, self.ncol, iplot, frameon=False)
            if self.cut != 1:
                self.normed = False
                vmin = - max(abs(data[..., indPlot].max()), abs(data[..., indPlot].min()))
                vmin = self.cut*vmin
                vmax = -vmin
                cs = N.linspace(vmin, vmax, self.levels)
                im = ax.contourf(xx, yy, data[..., indPlot], cs, extend='both',
                                  cmap=cmap, aa=True)
            else:
                cs = self.levels
                im = ax.contourf(xx, yy, data[..., indPlot], cs, 
                                  cmap=cmap, aa=True)
            rad = gr.radius[indPlot] * (1. - gr.radratio)
            P.title('%s, $r/r_o$=%.3f, Ra=%.1e' % (label, rad, gr.ra), 
                    fontsize=10)
            P.axis('off')


            if self.field not in ['entropy', 's', 'S'] and self.normed is True:
                im.set_clim(-max(abs(data[..., indPlot].max()), 
                                 abs(data[..., indPlot].min())), 
                             max(abs(data[..., indPlot].max()),
                                 abs(data[..., indPlot].min())))
            iplot += 1
        os.chdir(self.workdir)

    def plotEquat(self):
        cmap = P.get_cmap(self.cm)
        iplot = 1
        for datadir in self.dataliste:
            os.chdir(self.workdir + '/' + datadir)
            print datadir
            try:
                if self.ave:
                    gr = MagicGraph(ivar=1, ave=True)
                else:
                    gr = MagicGraph(ivar=1)
                if self.field in ('Vs', 'vs'):
                    vr = gr.vr
                    vt = gr.vtheta
                    thlin = N.linspace(0., N.pi, gr.ntheta)
                    th3D = N.zeros((gr.npI, gr.ntheta, gr.nr), dtype='f')
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    data = vr * N.sin(th3D) + vt * N.cos(th3D)
                    label = 'Vs'
                elif self.field in ('Vz', 'vz'):
                    vr = gr.vr
                    vt = gr.vtheta
                    thlin = N.linspace(0., N.pi, gr.ntheta)
                    th3D = N.zeros((gr.npI, gr.ntheta, gr.nr), dtype='f')
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    data = vr * N.cos(th3D) - vt * N.sin(th3D)
                    label = 'Vz'
                elif self.field in ('vortz'):
                    philoc = N.linspace(0., 2.*N.pi/gr.minc, gr.npI)
                    rrloc, pphiloc = N.meshgrid(gr.radius, philoc)
                    dr = rderavg(rrloc*gr.vphi[:,gr.ntheta/2,:], spectral=False,
                                 eta=gr.radratio, exclude=True)
                    equator = 1./rrloc*(dr - phideravg(gr.vr[:, gr.ntheta/2, :]))
                    label = r'$\omega_z$'
                else:
                    data, label = selectField(gr, self.field)
            except AttributeError:
                continue


            label += ' Ra = %.1e' % gr.ra

            if self.field not in ('vortz'):
                equator = data[:, gr.ntheta/2,:]

            equator = symmetrize(equator, gr.minc)

            phi = N.linspace(0., 2.*N.pi, gr.nphi)
            rr, pphi = N.meshgrid(gr.radius, phi)
            xx = rr * N.cos(pphi)
            yy = rr * N.sin(pphi)

            ax = P.subplot(self.nrow,self.ncol,iplot, frameon=False)
            if self.bw:
                im = ax.contour(xx, yy, equator, self.levels, colors='k',
                                linewidths=0.5)
            else:
                if self.cut != 1:    
                    self.normed = False
                    vmin = - max(abs(equator.max()), abs(equator.min()))
                    vmin = self.cut*vmin
                    vmax = -vmin
                    cs = N.linspace(vmin, vmax, self.levels)
                    im = ax.contourf(xx, yy, equator, cs, extend='both',
                                      cmap=cmap)
                else:
                    cs = self.levels
                    im = ax.contourf(xx, yy, equator, cs, cmap=cmap)
            ax.plot(gr.radius[0] * N.cos(phi), gr.radius[0]*N.sin(phi), 'k-')
            ax.plot(gr.radius[-1] * N.cos(phi), gr.radius[-1]*N.sin(phi), 'k-')

            # Variable conductivity
            if hasattr(gr, 'nVarCond'):
                if gr.nVarCond == 2:
                    radi = gr.con_radratio * gr.radius[0]
                    ax.plot(radi*N.cos(phi), radi*N.sin(phi), 'k--')

            if hasattr(gr, 'strat'):
                if gr.strat >= 1:
                    tit1 = r"$N_\rho = %.0f$"  % gr.strat
                else:
                    tit1 = r"$N_\rho = 10^{-2}$"
            else:
                tit1 = datadir

            ax.text(0.5, 0.5, tit1, fontsize=18,
                              horizontalalignment='center',
                              verticalalignment='center',
                              transform = ax.transAxes)

            #P.title(label, fontsize=10)
            P.axis('off')
            #fig.colorbar(im)

            if self.field not in ['entropy', 's', 'S'] and self.normed is True:
                im.set_clim(-max(abs(equator.max()), abs(equator.min())), 
                             max(abs(equator.max()), abs(equator.min())))
            iplot += 1
        os.chdir(self.workdir)

    def plotAvg(self):
        """
        Plot the azimutal average of a given field.
        """
        cmap = P.get_cmap(self.cm)

        iplot = 1
        for datadir in self.dataliste:
            os.chdir(self.workdir + '/' + datadir)
            print datadir
            try:
                if self.ave:
                    gr = MagicGraph(ivar=1, ave=True)
                else:
                    gr = MagicGraph(ivar=1)
                if self.field in ('Vs', 'vs'):
                    vr = gr.vr
                    vt = gr.vtheta
                    thlin = N.linspace(0., N.pi, gr.ntheta)
                    th3D = N.zeros((gr.npI, gr.ntheta, gr.nr), dtype='f')
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    data = vr * N.sin(th3D) + vt * N.cos(th3D)
                    label = 'Vs'
                elif self.field in ('Vz', 'vz'):
                    vr = gr.vr
                    vt = gr.vtheta
                    thlin = N.linspace(0., N.pi, gr.ntheta)
                    th3D = N.zeros((gr.npI, gr.ntheta, gr.nr), dtype='f')
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    data = vr * N.cos(th3D) - vt * N.sin(th3D)
                    label = 'Vz'
                elif self.field in ('Cr', 'cr'):
                    vr = gr.vr
                    vt = gr.vtheta
                    vp = gr.vphi.copy()
                    vp = gr.vphi- gr.vphi.mean(axis=0) # convective vp
                    thlin = N.linspace(0., N.pi, gr.ntheta)
                    th3D = N.zeros((gr.npI, gr.ntheta, gr.nr), dtype='f')
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    vs = vr * N.sin(th3D) + vt * N.cos(th3D)
                    data =  vs * vp
                    denom = N.sqrt(N.mean(vs**2, axis=0)* N.mean(vp**2, axis=0))
                    label = r'$\langle v_s v_\phi\rangle$'
                else:
                    data, label = selectField(gr, self.field)
            except AttributeError:
                continue

            #label += ' Ra = %.1e' % gr.ra
            label = 'Ra = %.1e' % gr.ra
            
            if self.field not in ('Cr', 'cr', 'ra', 'ratio', 'Cz', 'cz'):
                phiavg = data.mean(axis=0)
            else:
                ro = gr.radius[0]
                ri = gr.radius[-1]
                fac = 2./(N.pi*(ro**2-ri**2))
                facOTC = ro**2.*(N.pi-2.*N.arcsin(gr.radratio))/2. \
                         -ri**2*N.sqrt(1.-gr.radratio**2)/gr.radratio
                facOTC = 1./facOTC
                facITC = ri**2*N.sqrt(1.-gr.radratio**2)/gr.radratio \
                         +(ro**2-ri**2)* N.arcsin(gr.radratio) \
                         -ri**2/2.*(N.pi - 2.*N.arcsin(gr.radratio))
                facITC = 1./facITC
                phiavg = data.mean(axis=0)

                TC = N.array([], dtype='f')
                outTC = N.array([], dtype='f')
                inTC = N.array([], dtype='f')
                integ = N.array([], dtype='f')
                for k, th in enumerate(gr.colatitude):
                    rr = gr.radius[::-1]
                    dat = phiavg[k, ::-1] * rr
                    corr = intcheb(dat, gr.nr-1, ri, ro)
                    TC = N.append(TC, corr)
                    if th >= N.arcsin(gr.radratio)  and \
                       th <= N.pi - N.arcsin(gr.radratio):
                        # Outside tangent cylinder
                        val = trapz(dat[rr >= ri/N.sin(th)], 
                                    rr[rr >= ri/N.sin(th)])
                        outTC = N.append(outTC, val)
                        integ = N.append(integ, th)
                        # Inside tangent cylinder
                        val = trapz(dat[rr < ri/N.sin(th)], 
                                    rr[rr < ri/N.sin(th)])
                        inTC = N.append(inTC, val)
                    else:
                        val= intcheb(dat, gr.nr-1, ri, ro)
                        inTC = N.append(inTC, val)

                mask = N.where(denom == 0, 1, 0)
                phiavg /= (denom+mask)


            th = N.linspace(0., N.pi, gr.ntheta)
            rr, tth = N.meshgrid(gr.radius, th)
            xx = rr * N.sin(tth)
            yy = rr * N.cos(tth)

            titmax = phiavg.max()
            titmin = phiavg.min()

            ax = P.subplot(self.nrow, self.ncol, iplot, frameon=False)
            if self.bw:
                im = ax.contour(xx, yy, phiavg, self.levels, colors='k', 
                                linewidths=0.5)
            else:
                if self.cut != 1:    
                    self.normed = False
                    vmin = - max(abs(phiavg.max()), abs(phiavg.min()))
                    vmin = self.cut*vmin
                    vmax = -vmin
                    cs = N.linspace(vmin, vmax, self.levels)
                    im = ax.contourf(xx, yy, phiavg, cs, extend='both',
                                      cmap=cmap)
                else:
                    cs = self.levels
                    im = ax.contourf(xx, yy, phiavg, cs, cmap=cmap)
            ax.plot(gr.radius[0]*N.sin(th), gr.radius[0]*N.cos(th),
                   'k-')
            ax.plot(gr.radius[-1]*N.sin(th), gr.radius[-1]*N.cos(th),
                   'k-')

            P.plot([0., 0], [gr.radius[-1], gr.radius[0]], 'k-')
            P.plot([0., 0], [-gr.radius[-1], -gr.radius[0]], 'k-')
            
            # Variable conductivity
            if hasattr(gr, 'nVarCond'):
                if gr.nVarCond == 2:
                    radi = gr.con_radratio * gr.radius[0]
                    ax.plot(radi*N.sin(th), radi*N.cos(th), 'k--')

            P.title(label, fontsize=12)
            P.axis('off')
            #fig.colorbar(im)

            """
            if gr.strat >= 1:
                tit1 = r"$N_\rho = %.0f$"  % gr.strat
            else:
                tit1 = r"$N_\rho = 10^{-2}$"
            """
            #P.title(tit1, fontsize=12)

            """
            if int(titmin) == 0:
                tit1 = r'$+%i$' % titmax +'\n'+r'$%.1f$' % titmin
            else:
                tit1 = r'$+%i$' % titmax +'\n'+r'$%i$' % titmin
            ax.text(0., 0.5, tit1, fontsize=12,
                              horizontalalignment='left',
                              verticalalignment='center',
                              transform = ax.transAxes)
            tit2 = r'$+%.1e/-%.1e$' % (titmax, titmin)
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
        cmap = P.get_cmap(self.cm)

        iplot = 1
        for datadir in self.dataliste:
            os.chdir(self.workdir + '/' + datadir)
            print datadir
            try:
                if self.ave:
                    gr = MagicGraph(ivar=1, ave=True)
                else:
                    gr = MagicGraph(ivar=1)
                if self.field in ('Vs', 'vs'):
                    vr = gr.vr
                    vt = gr.vtheta
                    thlin = N.linspace(0., N.pi, gr.ntheta)
                    th3D = N.zeros((gr.npI, gr.ntheta, gr.nr), dtype='f')
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    data = vr * N.sin(th3D) + vt * N.cos(th3D)
                    label = 'Vs'
                elif self.field in ('Vz', 'vz'):
                    vr = gr.vr
                    vt = gr.vtheta
                    thlin = N.linspace(0., N.pi, gr.ntheta)
                    th3D = N.zeros((gr.npI, gr.ntheta, gr.nr), dtype='f')
                    for i in range(gr.ntheta):
                        th3D[:, i, :] = thlin[i]
                    data = vr * N.cos(th3D) - vt * N.sin(th3D)
                    label = 'Vz'
                else:
                    data, label = selectField(gr, self.field)
            except AttributeError:
                continue

            #label += ' Ra = %.1e' % gr.ra
            label = 'Ra = %.1e' % gr.ra
            
            phiavg = data[0, ...]

            th = N.linspace(0., N.pi, gr.ntheta)
            rr, tth = N.meshgrid(gr.radius, th)
            xx = rr * N.sin(tth)
            yy = rr * N.cos(tth)

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

            ax = P.subplot(self.nrow, self.ncol, iplot, frameon=False)
            if self.cut != 1:    
                self.normed = False
                vmin = - max(abs(phiavg.max()), abs(phiavg.min()))
                vmin = self.cut*vmin
                vmax = -vmin
                cs = N.linspace(vmin, vmax, self.levels)
                im = ax.contourf(xx, yy, phiavg, cs, extend='both', cmap=cmap)
            else:
                cs = self.levels
                im = ax.contourf(xx, yy, phiavg, cs, cmap=cmap)
            ax.plot(gr.radius[0]*N.sin(th), gr.radius[0]*N.cos(th),
                   'k-')
            ax.plot(gr.radius[-1]*N.sin(th), gr.radius[-1]*N.cos(th),
                   'k-')

            P.plot([0., 0], [gr.radius[-1], gr.radius[0]], 'k-')
            P.plot([0., 0], [-gr.radius[-1], -gr.radius[0]], 'k-')
            
            # Variable conductivity
            if hasattr(gr, 'nVarCond'):
                if gr.nVarCond == 2:
                    radi = gr.con_radratio * gr.radius[0]
                    ax.plot(radi*N.sin(th), radi*N.cos(th), 'k--')

            P.title(label, fontsize=12)
            P.axis('off')
            #fig.colorbar(im)

            if gr.strat >= 1:
                tit1 = r"$N_\rho = %.0f$"  % gr.strat
            else:
                tit1 = r"$N_\rho = 10^{-2}$"
            #P.title(tit1, fontsize=12)

            if int(titmin) == 0:
                tit1 = r'$+%i$' % titmax +'\n'+r'$%.1f$' % titmin
            else:
                tit1 = r'$+%i$' % titmax +'\n'+r'$%i$' % titmin
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
