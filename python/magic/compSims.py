# -*- coding: utf-8 -*-
from magic import *
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import numpy as np
import os
try:
    from scipy.integrate import trapz
except:
    from scipy.integrate import trapezoid as trapz

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

    def __init__(self, file='listRuns', field='ts', ncol=4, cm=None, dpi=96,
                 normed=True, levels=65, type=None, fullPath=False,
                 r=0.9, bw=False, ave=False, cut=1):
        """
        :param file: the input file that contains the list of directories that one
                     wants to analyse
         or self.field == 'heat':type file: str
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
        :param fullPath: set to True if the full path is specified in the input file
        :type fullPath: bool
        :param dpi: dot per inch when saving PNGs
        :type dpi: int
        :param bw: when set to True, display grey-scaled contour levels
        :type bw: bool
        :param cut: adjust the contour extrema to max(abs(data))*cut
        :type cut: float
        """
        self.list_of_runs = []
        self.workdir = os.getcwd()
        self.fullPath = fullPath
        self.field = field
        self.cm = cm
        self.normed = normed
        self.cut = cut
        self.levels = levels
        self.r = r
        self.bw = bw # for black and white outputs
        self.ave = ave # for G_ave.TAG files
        with open(file, 'r') as f:
            for line in f.readlines():
                self.list_of_runs.append(line.strip())
        self.ncol = ncol
        self.nplot = len(self.list_of_runs)
        if (self.nplot % self.ncol != 0):
            self.nrow = self.nplot // self.ncol + 1
        else:
            self.nrow = self.nplot // self.ncol

        if type == 'avg' or type == 'slice':
            self.fig = plt.figure(figsize=(self.ncol*1.5, self.nrow*3), dpi=dpi)
        elif type == 'equat':
            self.fig = plt.figure(figsize=(self.ncol*2.5, self.nrow*2.5), dpi=dpi)
        elif type == 'surf':
            self.fig = plt.figure(figsize=(self.ncol*3, self.nrow*1.7), dpi=dpi)
        else:
            self.fig = plt.figure(figsize=(self.ncol*3, self.nrow*3), dpi=dpi)
        if self.nrow == 1:
            if type == 'surf':
                self.fig.subplots_adjust(left=0.05, right=0.98, top=0.85, bottom=0.05)
            else:
                self.fig.subplots_adjust(left=0.05, right=0.98, top=0.92, bottom=0.05)
        else:
            if type == 'surf':
                self.fig.subplots_adjust(left=0.05, right=0.98, top=0.92, bottom=0.05)
            else:
                self.fig.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.05)

        if self.field == 'ts':
            self.plotTs()
        elif self.field == 'e_mag':
            self.plotEmag()
        elif self.field == 'flux' or self.field == 'heat':
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

        self.fig.tight_layout()

    def plotTs(self):
        """
        Plot time-series of the kinetic energy
        """
        iplot = 1
        for datadir in self.list_of_runs:
            if not self.fullPath:
                os.chdir(os.path.join(self.workdir, datadir))
            else:
                os.chdir(datadir)
            print(datadir)
            ts = MagicTs(field='e_kin', iplot=False, all=True)

            ax = self.fig.add_subplot(self.nrow, self.ncol, iplot)
            if hasattr(ts, 'time'):
                ax.semilogy(ts.time, ts.ekin_pol, c='C0')
                ax.semilogy(ts.time, ts.ekin_tor, c='C1')
                ax.semilogy(ts.time, ts.ekin_pol_axi, c='C0', ls='--')
                ax.semilogy(ts.time, ts.ekin_tor_axi, c='C1', ls='--')
                ax.semilogy(ts.time, ts.ekin_tot, 'k-')

                ax.set_title('Ra = {:.1e}'.format(ts.ra), fontsize=10)
                ax.set_xlim((ts.time.min(), ts.time.max()))
            iplot += 1
        os.chdir(self.workdir)

    def plotEmag(self):
        """
        Plot time-series of the magnetic energy
        """
        iplot = 1
        for datadir in self.list_of_runs:
            if not self.fullPath:
                os.chdir(os.path.join(self.workdir, datadir))
            else:
                os.chdir(datadir)
            print(datadir)
            ts = MagicTs(field='e_mag_oc', iplot=False, all=True)

            ax = self.fig.add_subplot(self.nrow, self.ncol, iplot)
            if hasattr(ts, 'time'):
                ax.semilogy(ts.time, ts.emagoc_pol, c='C0')
                ax.semilogy(ts.time, ts.emagoc_tor, c='C1')
                ax.semilogy(ts.time, ts.emagoc_pol_axi, c='C0', ls='--')
                ax.semilogy(ts.time, ts.emagoc_tor_axi, c='C1', ls='--')
                ax.semilogy(ts.time, ts.emag_tot, 'k-')

                ax.set_title('Ra = {:.1e}, Pm= {:.1f}'.format(ts.ra, ts.prmag),
                             fontsize=10)
                ax.set_xlim((ts.time.min(), ts.time.max()))
            iplot += 1
        os.chdir(self.workdir)

    def plotFlux(self):
        """
        Plot time-series of the top and bottom Nusselt numbers
        """
        iplot = 1
        for datadir in self.list_of_runs:
            if not self.fullPath:
                os.chdir(os.path.join(self.workdir, datadir))
            else:
                os.chdir(datadir)
            print(datadir)
            ts = MagicTs(field='heat', iplot=False, all=True)

            ax = self.fig.add_subplot(self.nrow, self.ncol, iplot)
            if hasattr(ts, 'time'):
                if ts.botnuss.max() == 1:
                    ax.plot(ts.time, ts.deltaTnuss)
                else:
                    ax.plot(ts.time, ts.botnuss)
                    ax.plot(ts.time, ts.topnuss)

                ax.set_title('Ra = {:.1e}'.format(ts.ra), fontsize=10)
                ax.set_xlim((ts.time.min(), ts.time.max()))
            iplot += 1
        os.chdir(self.workdir)

    def plotZonal(self):
        """
        Plot surface zonal flow profiles.
        """
        iplot = 1
        for datadir in self.list_of_runs:
            if not self.fullPath:
                os.chdir(os.path.join(self.workdir, datadir))
            else:
                os.chdir(datadir)
            if self.ave:
                gr = MagicGraph(ivar=1, ave=True)
            else:
                gr = MagicGraph(ivar=1)
            ax = self.fig.add_subplot(self.nrow, self.ncol, iplot)
            vpm = gr.vphi.mean(axis=0)
            theta = np.linspace(-90., 90, gr.ntheta)

            ax.plot(vpm[:, 1], theta)

            roequat = vpm[gr.ntheta//2, 0]*gr.ek*(1.-gr.radratio)
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
        if self.cm is None:
            self.cm = default_cmap(self.field)
        cmap = plt.get_cmap(self.cm)

        iplot = 1
        for datadir in self.list_of_runs:
            if not self.fullPath:
                os.chdir(os.path.join(self.workdir, datadir))
            else:
                os.chdir(datadir)
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

            data = symmetrize(data[..., indPlot], gr.minc)

            ax = self.fig.add_subplot(self.nrow, self.ncol, iplot, frameon=False)
            radialContour(data, fig=self.fig, ax=ax, levels=self.levels,
                          cm=cmap, cbar=False)

            iplot += 1
        os.chdir(self.workdir)

    def plotEquat(self):
        """
        Plot equatorial cuts in (phi, r)  planes.
        """
        if self.cm is None:
            self.cm = default_cmap(self.field)
        cmap = plt.get_cmap(self.cm)

        iplot = 1
        for datadir in self.list_of_runs:
            if not self.fullPath:
                os.chdir(os.path.join(self.workdir, datadir))
            else:
                os.chdir(datadir)
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
                    dr = rderavg(rrloc*gr.vphi[:,gr.ntheta//2,:],
                                 gr.radius, exclude=True)
                    equator = 1./rrloc*(dr - phideravg(gr.vr[:, gr.ntheta//2, :],
                                                       gr.minc))
                    if labTex:
                        label = r'$\omega_z$'
                    else:
                        label = 'omega'
                else:
                    data, data_ic, label = selectField(gr, self.field)
            except AttributeError:
                continue

            if self.field not in ('vortz'):
                equator = data[:, gr.ntheta//2,:]

            equator = symmetrize(equator, gr.minc)

            ax = self.fig.add_subplot(self.nrow, self.ncol, iplot, frameon=False)
            equatContour(equator, gr.radius, levels=self.levels, cm=cmap,
                         cbar=False, fig=self.fig, ax=ax)

            iplot += 1
        os.chdir(self.workdir)

    def plotAvg(self):
        """
        Plot azimutal averages in (theta, r) planes.
        """
        if self.cm is None:
            self.cm = default_cmap(self.field)
        cmap = plt.get_cmap(self.cm)

        iplot = 1
        for datadir in self.list_of_runs:
            if not self.fullPath:
                os.chdir(os.path.join(self.workdir, datadir))
            else:
                os.chdir(datadir)
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

            ax = self.fig.add_subplot(self.nrow, self.ncol, iplot, frameon=False)

            merContour(phiavg, gr.radius, fig=self.fig, ax=ax, cbar=False,
                       cm=cmap, levels=self.levels)

            iplot += 1
        os.chdir(self.workdir)

    def plotSlice(self):
        if self.cm is None:
            self.cm = default_cmap(self.field)
        cmap = plt.get_cmap(self.cm)

        iplot = 1
        for datadir in self.list_of_runs:
            if not self.fullPath:
                os.chdir(os.path.join(self.workdir, datadir))
            else:
                os.chdir(datadir)
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

            phiavg = data[0, ...]

            ax = self.fig.add_subplot(self.nrow, self.ncol, iplot, frameon=False)

            merContour(phiavg, gr.radius, levels=self.levels, cm=cmap,
                       fig=self.fig, ax=ax, cbar=False)

            iplot += 1
        os.chdir(self.workdir)

if __name__ == '__main__':
    CompSims(file='listRuns', field='ts', ncol=5, type=None)
