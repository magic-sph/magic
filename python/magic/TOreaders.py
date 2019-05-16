# -*- coding: utf-8 -*-
import glob
import re
import os
import numpy as np
import matplotlib.pyplot as plt
from .npfile import *
from .log import MagicSetup
from .libmagic import scanDir


class TOMovie:
    """
    This class allows to read and display the :ref:`TO_mov.TAG <secTO_movieFile>`
    generated when :ref:`l_TOmovie=.true. <varl_TOmovie>` is True.

    >>> # This will allow you to pick up one TO_mov files among the existing ones
    >>> t = TOMovie()

    >>> # Read TO_mov.N0m2, time-averaged it and display it with 65 contour levels
    >>> t = TOMovie(file='TO_mov.N0m2', avg=True, levels=65, cm='seismic')
    """

    def __init__(self, file=None, iplot=True, cm='RdYlBu_r',
                 cut=0.8, levels=16, avg=True, precision=np.float32):
        """
        :param file: the filename of the TO_mov file
        :type file: str
        :param cmap: the name of the color map
        :type cmap: str
        :param levels: the number of contour levels
        :type levels: int
        :param cut: adjust the contour extrema to max(abs(data))*cut
        :type cut: float
        :param iplot: a boolean to specify if one wants to plot or not the
                      results
        :type iplot: bool
        :param avg: time average of the different forces
        :type avg: bool
        :param precision: precision of the input file, np.float32 for single
                          precision, np.float64 for double precision
        :type precision: str
        """

        if file is None:
            dat = glob.glob('TO_mov.*')
            str1 = 'Which TO movie do you want ?\n'
            for k, movie in enumerate(dat):
                str1 += ' %i) %s\n' % (k+1, movie)
            index = int(input(str1))
            try:
                filename = dat[index-1]
            except IndexError:
                print('Non valid index: %s has been chosen instead' % dat[0])
                filename = dat[0]

        else:
            filename = file

        mot = re.compile(r'TO_mov\.(.*)')
        end = mot.findall(filename)[0]

        # DETERMINE THE NUMBER OF LINES BY READING THE LOG FILE
        logfile = open('log.%s' % end, 'r')
        mot = re.compile(r' ! WRITING TO MOVIE FRAME NO\s*(\d*).*')
        for line in logfile.readlines():
            if mot.match(line):
                nlines = int(mot.findall(line)[0])
        logfile.close()

        self.nvar = nlines

        # READ the movie file
        infile = npfile(filename, endian='B')
        # HEADER
        version = infile.fort_read('|S64')
        n_type, n_surface, const, n_fields = infile.fort_read(precision)
        movtype = infile.fort_read(precision)
        n_fields = int(n_fields)
        self.movtype = np.asarray(movtype)
        n_surface = int(n_surface)

        # RUN PARAMETERS
        runid = infile.fort_read('|S64')
        n_r_mov_tot, n_r_max, n_theta_max, n_phi_tot, self.minc, self.ra, \
            self.ek, self.pr, self.prmag, \
            self.radratio, self.tScale = infile.fort_read(precision)
        n_r_mov_tot = int(n_r_mov_tot)
        self.n_r_max = int(n_r_max)
        self.n_theta_max = int(n_theta_max)
        self.n_phi_tot = int(n_phi_tot)

        # GRID
        self.radius = infile.fort_read(precision)
        self.radius = self.radius[:self.n_r_max]  # remove inner core
        self.theta = infile.fort_read(precision)
        self.phi = infile.fort_read(precision)

        surftype = 'phi_constant'
        shape = (self.n_r_max, self.n_theta_max)

        self.time = np.zeros(self.nvar, precision)
        self.asVphi = np.zeros((self.nvar, self.n_theta_max, self.n_r_max),
                               precision)
        self.rey = np.zeros_like(self.asVphi)
        self.adv = np.zeros_like(self.asVphi)
        self.visc = np.zeros_like(self.asVphi)
        self.lorentz = np.zeros_like(self.asVphi)
        self.coriolis = np.zeros_like(self.asVphi)
        self.dtVp = np.zeros_like(self.asVphi)

        # READ the data
        for k in range(self.nvar):
            n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                movieDipLon, movieDipStrength, movieDipStrengthGeo \
                = infile.fort_read(precision)
            self.time[k] = t_movieS
            self.asVphi[k, ...] = infile.fort_read(precision, shape=shape).T
            self.rey[k, ...] = infile.fort_read(precision, shape=shape).T
            self.adv[k, ...] = infile.fort_read(precision, shape=shape).T
            self.visc[k, ...] = infile.fort_read(precision, shape=shape).T
            self.lorentz[k, ...] = infile.fort_read(precision, shape=shape).T
            self.coriolis[k, ...] = infile.fort_read(precision, shape=shape).T
            self.dtVp[k, ...] = infile.fort_read(precision, shape=shape).T

        if iplot:
            cmap = plt.get_cmap(cm)
            self.plot(cut, levels, avg, cmap)

    def plot(self, cut=0.8, levs=16, avg=True, cmap='RdYlBu_r'):
        """
        Plotting function

        :param cut: adjust the contour extrema to max(abs(data))*cut
        :type cut: float
        :param levs: number of contour levels
        :type levs: int
        :param avg: when set to True, quantities are time-averaged
        :type avg: bool
        :param cmap: name of the colormap
        :type cmap: str
        """
        th = np.linspace(np.pi/2., -np.pi/2., self.n_theta_max)
        rr, tth = np.meshgrid(self.radius, th)
        xx = rr * np.cos(tth)
        yy = rr * np.sin(tth)
        xxout = rr.max() * np.cos(th)
        yyout = rr.max() * np.sin(th)
        xxin = rr.min() * np.cos(th)
        yyin = rr.min() * np.sin(th)
        fig = plt.figure(figsize=(20, 5))
        fig.subplots_adjust(top=0.99, right=0.99, bottom=0.01, left=0.01,
                            wspace=0.01)

        if avg:
            cor = self.coriolis.mean(axis=0)
            asVp = self.asVphi.mean(axis=0)
            rey = self.rey.mean(axis=0)
            adv = self.adv.mean(axis=0)
            lor = self.lorentz.mean(axis=0)
            vis = self.visc.mean(axis=0)
            dt = self.dtVp[:-1].mean(axis=0)

            vmin = - max(abs(cor.max()), abs(cor.min()))
            vmin = cut * vmin
            vmax = -vmin
            cs = np.linspace(vmin, vmax, levs)

            vmin = - max(abs(asVp.max()), abs(asVp.min()))
            vmin = cut * vmin
            vmax = -vmin
            csVp = np.linspace(vmin, vmax, levs)

            ax = fig.add_subplot(181)
            ax.axis('off')
            im = ax.contourf(xx, yy, asVp, csVp, cmap=cmap, extend='both')
            ax.plot(xxout, yyout, 'k-', lw=1.5)
            ax.plot(xxin, yyin, 'k-', lw=1.5)
            ax.text(0.05, 0., 'uphi', fontsize=20, horizontalalignment='left',
                    verticalalignment='center')
            ax = fig.add_subplot(182)
            ax.axis('off')
            im = ax.contourf(xx, yy, adv, cs, cmap=cmap, extend='both')
            ax.plot(xxout, yyout, 'k-', lw=1.5)
            ax.plot(xxin, yyin, 'k-', lw=1.5)
            ax.text(0.05, 0., 'Adv', fontsize=20, horizontalalignment='left',
                    verticalalignment='center')
            ax = fig.add_subplot(183)
            ax.axis('off')
            im = ax.contourf(xx, yy, rey, cs, cmap=cmap, extend='both')
            ax.plot(xxout, yyout, 'k-', lw=1.5)
            ax.plot(xxin, yyin, 'k-', lw=1.5)
            ax.text(0.05, 0., 'Rey', fontsize=20, horizontalalignment='left',
                    verticalalignment='center')
            ax = fig.add_subplot(184)
            ax.axis('off')
            im = ax.contourf(xx, yy, vis, cs, cmap=cmap, extend='both')
            ax.plot(xxout, yyout, 'k-', lw=1.5)
            ax.plot(xxin, yyin, 'k-', lw=1.5)
            ax.text(0.05, 0., 'Visc', fontsize=20, horizontalalignment='left',
                    verticalalignment='center')
            ax = fig.add_subplot(185)
            ax.axis('off')
            im = ax.contourf(xx, yy, lor, cs, cmap=cmap, extend='both')
            ax.plot(xxout, yyout, 'k-', lw=1.5)
            ax.plot(xxin, yyin, 'k-', lw=1.5)
            ax.text(0.05, 0., 'Lo.', fontsize=20, horizontalalignment='left',
                    verticalalignment='center')
            ax = fig.add_subplot(186)
            ax.axis('off')
            im = ax.contourf(xx, yy, cor, cs, cmap=cmap, extend='both')
            ax.plot(xxout, yyout, 'k-', lw=1.5)
            ax.plot(xxin, yyin, 'k-', lw=1.5)
            ax.text(0.05, 0., 'Cor.', fontsize=20, horizontalalignment='left',
                    verticalalignment='center')
            ax = fig.add_subplot(187)
            ax.axis('off')
            balance = adv+cor+vis+lor+rey
            im = ax.contourf(xx, yy, balance, cs, cmap=cmap, extend='both')
            ax.plot(xxout, yyout, 'k-', lw=1.5)
            ax.plot(xxin, yyin, 'k-', lw=1.5)
            ax.text(0.05, 0., 'sum', fontsize=20, horizontalalignment='left',
                    verticalalignment='center')
            ax = fig.add_subplot(188)
            ax.axis('off')
            im = ax.contourf(xx, yy, dt, cs, cmap=cmap, extend='both')
            ax.plot(xxout, yyout, 'k-', lw=1.5)
            ax.plot(xxin, yyin, 'k-', lw=1.5)
            ax.text(0.01, 0., 'dvp/dt', fontsize=20,
                    horizontalalignment='left', verticalalignment='center')

        else:
            plt.ion()

            vmin = - max(abs(self.coriolis.max()), abs(self.coriolis.min()))
            vmin = cut * vmin
            vmax = -vmin
            cs = np.linspace(vmin, vmax, levs)

            vmin = - max(abs(self.asVphi.max()), abs(self.asVphi.min()))
            vmin = cut * vmin
            vmax = -vmin
            cs1 = np.linspace(vmin, vmax, levs)

            for k in range(self.nvar-1):  # avoid last dvp/dt which is wrong
                bal = self.asVphi[k, ...]+self.adv[k, ...]+self.rey[k, ...] + \
                      self.visc[k, ...]+self.lorentz[k, ...] + \
                      self.coriolis[k, ...]
                if k == 0:
                    ax1 = fig.add_subplot(181)
                    ax1.axis('off')
                    im = ax1.contourf(xx, yy, self.asVphi[k, ...], cs1,
                                      cmap=cmap, extend='both')
                    ax1.plot(xxout, yyout, 'k-', lw=1.5)
                    ax1.plot(xxin, yyin, 'k-', lw=1.5)
                    ax1.text(0.05, 0., 'uphi', fontsize=20,
                             horizontalalignment='left',
                             verticalalignment='center')
                    ax2 = fig.add_subplot(182)
                    ax2.axis('off')
                    im = ax2.contourf(xx, yy, self.adv[k, ...], cs,
                                      cmap=cmap, extend='both')
                    ax2.plot(xxout, yyout, 'k-', lw=1.5)
                    ax2.plot(xxin, yyin, 'k-', lw=1.5)
                    ax2.text(0.05, 0., 'Adv', fontsize=20,
                             horizontalalignment='left',
                             verticalalignment='center')
                    ax3 = fig.add_subplot(183)
                    ax3.axis('off')
                    im = ax3.contourf(xx, yy, self.rey[k, ...], cs,
                                      cmap=cmap, extend='both')
                    ax3.plot(xxout, yyout, 'k-', lw=1.5)
                    ax3.plot(xxin, yyin, 'k-', lw=1.5)
                    ax3.text(0.05, 0., 'Rey', fontsize=20,
                             horizontalalignment='left',
                             verticalalignment='center')
                    ax4 = fig.add_subplot(184)
                    ax4.axis('off')
                    im = ax4.contourf(xx, yy, self.visc[k, ...], cs,
                                      cmap=cmap, extend='both')
                    ax4.plot(xxout, yyout, 'k-', lw=1.5)
                    ax4.plot(xxin, yyin, 'k-', lw=1.5)
                    ax4.text(0.05, 0., 'Visc', fontsize=20,
                             horizontalalignment='left',
                             verticalalignment='center')
                    ax5 = fig.add_subplot(185)
                    ax5.axis('off')
                    im = ax5.contourf(xx, yy, self.lorentz[k, ...], cs,
                                      cmap=cmap, extend='both')
                    ax5.plot(xxout, yyout, 'k-', lw=1.5)
                    ax5.plot(xxin, yyin, 'k-', lw=1.5)
                    ax5.text(0.05, 0., 'Lo.', fontsize=20,
                             horizontalalignment='left',
                             verticalalignment='center')
                    ax6 = fig.add_subplot(186)
                    ax6.axis('off')
                    im = ax6.contourf(xx, yy, self.coriolis[k, ...], cs,
                                      cmap=cmap, extend='both')
                    ax6.plot(xxout, yyout, 'k-', lw=1.5)
                    ax6.plot(xxin, yyin, 'k-', lw=1.5)
                    ax6.text(0.05, 0., 'Cor.', fontsize=20,
                             horizontalalignment='left',
                             verticalalignment='center')
                    ax7 = fig.add_subplot(187)
                    ax7.axis('off')
                    im = ax7.contourf(xx, yy, bal, cs,
                                      cmap=cmap, extend='both')
                    ax7.plot(xxout, yyout, 'k-', lw=1.5)
                    ax7.plot(xxin, yyin, 'k-', lw=1.5)
                    ax7.text(0.05, 0., 'sum', fontsize=20,
                             horizontalalignment='left',
                             verticalalignment='center')
                    ax8 = fig.add_subplot(188)
                    ax8.axis('off')
                    im = ax8.contourf(xx, yy, self.dtVp[k, ...], cs,
                                      cmap=cmap, extend='both')
                    ax8.plot(xxout, yyout, 'k-', lw=1.5)
                    ax8.plot(xxin, yyin, 'k-', lw=1.5)
                    ax8.text(0.05, 0., 'dvp/dt', fontsize=20,
                             horizontalalignment='left',
                             verticalalignment='center')

                    man = plt.get_current_fig_manager()
                    man.canvas.draw()
                else:
                    plt.cla()
                    ax1.contourf(xx, yy, self.asVphi[k, ...], cs1,
                                 cmap=cmap, extend='both')
                    ax2.contourf(xx, yy, self.adv[k, ...], cs,
                                 cmap=cmap, extend='both')
                    ax3.contourf(xx, yy, self.rey[k, ...], cs,
                                 cmap=cmap, extend='both')
                    ax4.contourf(xx, yy, self.visc[k, ...], cs,
                                 cmap=cmap, extend='both')
                    ax5.contourf(xx, yy, self.lorentz[k, ...], cs,
                                 cmap=cmap, extend='both')
                    ax6.contourf(xx, yy, self.coriolis[k, ...], cs,
                                 cmap=cmap, extend='both')
                    ax7.contourf(xx, yy, bal, cs,
                                 cmap=cmap, extend='both')
                    ax8.contourf(xx, yy, self.dtVp[k, ...], cs,
                                 cmap=cmap, extend='both')

                    plt.axis('off')
                    man.canvas.draw()


class MagicTOZ(MagicSetup):
    """
    This class can be used to read the TOZ.TAG files produced by the TO outputs

    >>> # read the content of TOZ_1.tag
    >>> # where tag is the most recent file in the current directory
    >>> toz = MagicTOZ(itoz=1)
    >>> # read the content of TOZ_ave.test
    >>> toz = MagicTOZ(tag='test', ave=True)
    """

    def __init__(self, datadir='.', itoz=None, tag=None, precision=np.float32,
                 ave=False):
        """
        :param datadir: current working directory
        :type datadir: str
        :param tag: file suffix (tag), if not specified the most recent one in
                    the current directory is chosen
        :type tag: str
        :param precision: single or double precision
        :type precision: str
        :param iplot: display the output plot when set to True (default is
                      True)
        :type iplot: bool
        :param ave: plot a time-averaged TOZ file when set to True
        :type ave: bool
        :param itoz: the number of the TOZ file you want to plot
        :type itoz: int
        """

        if ave:
            self.name = 'TOZ_ave'
            n_fields = 9
        else:
            self.name = 'TOZ_'
            n_fields = 8

        if tag is not None:
            if itoz is not None:
                file = '%s%i.%s' % (self.name, itoz, tag)
                filename = os.path.join(datadir, file)
            else:
                pattern = os.path.join(datadir, '%s*%s' % (self.name, tag))
                files = scanDir(pattern)
                if len(files) != 0:
                    filename = files[-1]
                else:
                    print('No such tag... try again')
                    return

            if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
        else:
            if itoz is not None:
                pattern = os.path.join(datadir, '%s%i*' % (self.name, itoz))
                files = scanDir(pattern)
                filename = files[-1]
            else:
                pattern = os.path.join(datadir, '%s*' % self.name)
                files = scanDir(pattern)
                filename = files[-1]
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists(os.path.join(datadir, 'log.%s' % ending)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % ending)


        infile = npfile(filename, endian='B')

        if ave:
            self.nSmax, self.omega_ic, self.omega_oc = infile.fort_read(precision)
        else:
            self.time, self.nSmax, self.omega_ic, \
                       self.omega_oc = infile.fort_read(precision)
        self.nSmax = int(self.nSmax)

        self.cylRad = infile.fort_read(precision)

        self.zall = np.zeros((2*self.nSmax, self.nSmax), dtype=precision)
        self.vp = np.zeros((2*self.nSmax, self.nSmax), dtype=precision)
        self.dvp = np.zeros((2*self.nSmax, self.nSmax), dtype=precision)
        self.Rstr = np.zeros((2*self.nSmax, self.nSmax), dtype=precision)
        self.Astr = np.zeros((2*self.nSmax, self.nSmax), dtype=precision)
        self.LF = np.zeros((2*self.nSmax, self.nSmax), dtype=precision)
        self.Cor = np.zeros((2*self.nSmax, self.nSmax), dtype=precision)
        self.str = np.zeros((2*self.nSmax, self.nSmax), dtype=precision)

        if n_fields == 9:
            self.CL = np.zeros((2*self.nSmax, self.nSmax), dtype=precision)

        for nS in range(self.nSmax):
            nZmaxNS = int(infile.fort_read(precision))
            data = infile.fort_read(precision)
            data = data.reshape((n_fields, nZmaxNS))
            self.zall[:nZmaxNS, nS] = data[0, :]
            self.vp[:nZmaxNS, nS] = data[1, :]
            self.dvp[:nZmaxNS, nS] = data[2, :]
            self.Rstr[:nZmaxNS, nS] = data[3, :]
            self.Astr[:nZmaxNS, nS] = data[4, :]
            self.LF[:nZmaxNS, nS] = data[5, :]
            self.str[:nZmaxNS, nS] = data[6, :]
            self.Cor[:nZmaxNS, nS] = data[7, :]
            if n_fields == 9:
                self.CL[:nZmaxNS, nS] = data[8, :]

        infile.close()

        S = np.zeros_like(self.zall)
        for i in range(2*self.nSmax):
            S[i, :] = self.cylRad


class MagicTOHemi(MagicSetup):
    """
    This class can be used to read and display z-integrated quantities
    produced by the TO outputs. Those are basically the TO[s|n]hn.TAG files

    >>> to = MagicTOHemi(hemi='n', iplot=True) # For the Northern hemisphere
    """

    def __init__(self, datadir='.', hemi='n', tag=None, precision=np.float32,
                 iplot=False):
        """
        :param datadir: current working directory
        :type datadir: str
        :param hemi: Northern or Southern hemisphere ('n' or 's')
        :type hemi: str
        :param tag: file suffix (tag), if not specified the most recent one in
                    the current directory is chosen
        :type tag: str
        :param precision: single or double precision
        :type precision: str
        :param iplot: display the output plot when set to True (default is
                      True)
        :type iplot: bool
        """

        if tag is not None:
            pattern = os.path.join(datadir, 'TO%shs.%s' % (hemi, tag))
            files = scanDir(pattern)
            if len(files) != 0:
                filename = files[-1]
            else:
                print('No such tag... try again')
                return

            if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
        else:
            pattern = os.path.join(datadir, 'TO%shs.*' % hemi)
            files = scanDir(pattern)
            filename = files[-1]
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists(os.path.join(datadir, 'log.%s' % ending)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % ending)

        if not os.path.exists(filename):
            print('No such file')
            return

        infile = npfile(filename, endian='B')

        self.nSmax = int(infile.fort_read(precision))
        self.cylRad = infile.fort_read(precision)

        ind = 0
        while 1:
            try:
                data = infile.fort_read(precision)
                if ind == 0:
                    self.time = np.r_[data[0]]
                else:
                    self.time = np.append(self.time, data[0])

                data = data[1:]
                data = data.reshape((17, self.nSmax))

                if ind == 0:
                    self.vp = data[0, :]
                    self.dvp = data[1, :]
                    self.ddvp = data[2, :]
                    self.vpr = data[3, :]
                    self.rstr = data[4, :]
                    self.astr = data[5, :]
                    self.LF = data[6, :]
                    self.viscstr = data[7, :]
                    self.tay = data[8, :]
                else:
                    self.vp = np.vstack((self.vp, data[0, :]))
                    self.dvp = np.vstack((self.dvp, data[1, :]))
                    self.ddvp = np.vstack((self.ddvp, data[2, :]))
                    self.vpr = np.vstack((self.vpr, data[3, :]))
                    self.rstr = np.vstack((self.rstr, data[4, :]))
                    self.astr = np.vstack((self.astr, data[5, :]))
                    self.LF = np.vstack((self.LF, data[6, :]))
                    self.viscstr = np.vstack((self.viscstr, data[7, :]))
                    self.tay = np.vstack((self.tay, data[8, :]))

                ind += 1
            except TypeError:
                break
        infile.close()

        if iplot:
            self.plot()

    def plot(self):
        """
        Plotting function
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.cylRad, self.vp.mean(axis=0))
        ax.set_xlabel('Cylindrical radius')
        ax.axhline(linestyle='--', linewidth=1.5, color='k')
        ax.set_ylabel('Vphi')
        ax.set_xlim(self.cylRad[0], self.cylRad[-1])

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.cylRad, self.rstr.mean(axis=0), label='Reynolds stress')
        ax.plot(self.cylRad, self.astr.mean(axis=0),
                label='Axi. Reynolds stress')
        ax.plot(self.cylRad, self.LF.mean(axis=0), label='Lorentz force')
        ax.plot(self.cylRad, self.viscstr.mean(axis=0), label='Viscous stress')
        ax.set_xlabel('Cylindrical radius')
        ax.set_ylabel('Integrated forces')
        ax.legend(loc='upper left', frameon=False)
        ax.set_xlim(self.cylRad[0], self.cylRad[-1])


class MagicTaySphere(MagicSetup):
    """
    This class can be used to read and display quantities
    produced by the TO outputs. Those are basically the TaySphere.TAG files

    >>> to = MagicTaySphere(iplot=True) # For the Northern hemisphere
    >>> print(to.time, to.e_kin)
    """

    def __init__(self, datadir='.', tag=None, precision=np.float32,
                 iplot=False):
        """
        :param datadir: current working directory
        :type datadir: str
        :param tag: file suffix (tag), if not specified the most recent one in
                    the current directory is chosen
        :type tag: str
        :param precision: single or double precision
        :type precision: str
        :param iplot: display the output plot when set to True (default is
                      True)
        :type iplot: bool
        """

        if tag is not None:
            pattern = os.path.join(datadir, 'TaySphere.%s' % tag)
            files = scanDir(pattern)
            if len(files) != 0:
                filename = files[-1]
            else:
                print('No such tag... try again')
                return

            if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
        else:
            pattern = os.path.join(datadir, 'TaySphere.*')
            files = scanDir(pattern)
            filename = files[-1]
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists(os.path.join(datadir, 'log.%s' % ending)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % ending)

        if not os.path.exists(filename):
            print('No such file')
            return

        infile = npfile(filename, endian='B')

        self.n_r_max = int(infile.fort_read(precision))
        self.radius = infile.fort_read(precision)

        ind = 0
        while 1:
            try:
                data = infile.fort_read(precision)
                if ind == 0:
                    self.time = np.r_[data[0]]
                    self.vpRMS = np.r_[data[1]]
                    self.vgRMS = np.r_[data[2]]
                    self.tayRMS = np.r_[data[3]]
                    self.taySRMS = np.r_[data[4]]
                    self.tayRRMS = np.r_[data[5]]
                    self.tayVRMS = np.r_[data[6]]
                    self.e_kin = np.r_[data[7]]
                else:
                    self.time = np.append(self.time, data[0])
                    self.vpRMS = np.append(self.vgRMS, data[1])
                    self.vgRMS = np.append(self.vgRMS, data[2])
                    self.tayRMS = np.append(self.tayRMS, data[3])
                    self.taySRMS = np.append(self.taySRMS, data[4])
                    self.tayRRMS = np.append(self.tayRRMS, data[5])
                    self.tayVRMS = np.append(self.tayVRMS, data[6])
                    self.e_kin = np.append(self.e_kin, data[7])

                data = data[8:]
                data = data.reshape((8, self.n_r_max))

                if ind == 0:
                    self.vp = data[0, :]
                    self.dvp = data[1, :]
                    self.rstr = data[2, :]
                    self.astr = data[3, :]
                    self.LF = data[4, :]
                    self.viscstr = data[5, :]
                    self.cor = data[6, :]
                    self.tay = data[7, :]
                else:
                    self.vp = np.vstack((self.vp, data[0, :]))
                    self.dvp = np.vstack((self.dvp, data[1, :]))
                    self.rstr = np.vstack((self.rstr, data[2, :]))
                    self.astr = np.vstack((self.astr, data[3, :]))
                    self.LF = np.vstack((self.LF, data[4, :]))
                    self.viscstr = np.vstack((self.viscstr, data[5, :]))
                    self.cor = np.vstack((self.cor, data[6, :]))
                    self.tay = np.vstack((self.tay, data[7, :]))

                ind += 1
            except TypeError:
                break
        infile.close()

        if iplot:
            self.plot()

    def plot(self):
        """
        Plotting function
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.radius, self.vp.mean(axis=0))
        ax.set_xlabel('Radius')
        ax.axhline(linestyle='--', linewidth=1.5, color='k')
        ax.set_ylabel('Vphi')
        ax.set_xlim(self.radius[0], self.radius[-1])

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.radius, self.rstr.mean(axis=0), label='Reynolds stress')
        ax.plot(self.radius, self.astr.mean(axis=0),
                label='Axi. Reynolds stress')
        ax.plot(self.radius, self.LF.mean(axis=0), label='Lorentz force')
        ax.plot(self.radius, self.viscstr.mean(axis=0), label='Viscous stress')
        ax.set_xlabel('Radius')
        ax.set_ylabel('Forces')
        ax.legend(loc='upper left', frameon=False)
        ax.set_xlim(self.radius[0], self.radius[-1])


class MagicPV(MagicSetup):
    """
    This class can be used to read and display quantities
    produced by the PV outputs. Those are basically the PVZ.TAG and the Vcy.TAG
    files

    >>> # To plot the content of PVZ.test
    >>> pv = MagicPV(field='PVZ', tag='test', iplot=True)
    """

    def __init__(self, field='PVZ', datadir='.', tag=None, iplot=False,
                 precision=np.float32):
        """
        :param field: file prefix (either 'Vcy' or 'PVZ')
        :type field: str
        :param datadir: current working directory
        :type datadir: str
        :param tag: file suffix (tag), if not specified the most recent one in
                    the current directory is chosen
        :type tag: str
        :param precision: single or double precision
        :type precision: str
        :param iplot: display the output plot when set to True (default is
                      True)
        :type iplot: bool
        """

        if tag is not None:
            pattern = os.path.join(datadir, '%s.%s' % (field, tag))
            files = scanDir(pattern)
            if len(files) != 0:
                filename = files[-1]
            else:
                print('No such tag... try again')
                return

            if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
        else:
            pattern = os.path.join(datadir, '%s.*' % field)
            files = scanDir(pattern)
            filename = files[-1]
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists(os.path.join(datadir, 'log.%s' % ending)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % ending)

        if not os.path.exists(filename):
            print('No such file')
            return

        infile = npfile(filename, endian='B')

        if field == 'Vcy':
            self.time, self.nSmax, self.nZmax, self.n_phi_max, self.omega_ic, \
                self.omega_ma, self.radratio, self.minc = \
                infile.fort_read(precision)
            self.n_phi_max = int(self.n_phi_max)
        elif field == 'PVZ':
            self.time, self.nSmax, self.nZmax, self.omega_ic, self.omega_ma = \
                                 infile.fort_read(precision)
        self.nSmax = int(self.nSmax)
        self.nZmax = int(self.nZmax)

        self.cyl_rad = infile.fort_read(precision)
        self.z = infile.fort_read(precision)

        if field == 'Vcy':
            self.vs = np.zeros((self.nZmax*self.n_phi_max, self.nSmax),
                               dtype=precision)
            self.vphi = np.zeros_like(self.vs)
            self.vz = np.zeros_like(self.vs)
            self.vort = np.zeros_like(self.vs)
            self.dtvort = np.zeros_like(self.vs)
            for i in range(self.nSmax):
                nZmax = int(infile.fort_read(precision))
                start = (self.nZmax-nZmax)/2 * self.n_phi_max
                stop = start+nZmax*self.n_phi_max
                self.vs[:nZmax*self.n_phi_max, i] = infile.fort_read(precision)
                self.vphi[start:stop, i] = infile.fort_read(precision)
                self.vz[:nZmax*self.n_phi_max, i] = infile.fort_read(precision)
                self.vort[:nZmax*self.n_phi_max, i] = \
                    infile.fort_read(precision)
                self.dtvort[:nZmax*self.n_phi_max, i] = \
                    infile.fort_read(precision)

            self.vs = self.vs.reshape((self.nZmax, self.n_phi_max, self.nSmax))
            self.vz = self.vz.reshape((self.nZmax, self.n_phi_max, self.nSmax))
            self.vphi = self.vphi.reshape((self.nZmax, self.n_phi_max,
                                           self.nSmax))
            self.vort = self.vort.reshape((self.nZmax, self.n_phi_max,
                                           self.nSmax))
            self.dtvort = self.dtvort.reshape((self.nZmax, self.n_phi_max,
                                               self.nSmax))

            self.vphi = np.transpose(self.vphi, axes=(1, 0, 2))
            self.vs = np.transpose(self.vs, axes=(1, 0, 2))
            self.vz = np.transpose(self.vz, axes=(1, 0, 2))
            self.vort = np.transpose(self.vort, axes=(1, 0, 2))
            self.dtvort = np.transpose(self.dtvort, axes=(1, 0, 2))

        elif field == 'PVZ':
            self.omS = np.zeros((self.nZmax, self.nSmax), dtype=precision)
            for i in range(self.nSmax):
                self.omS[:, i] = infile.fort_read(precision)

        if iplot:
            self.plot(field)

    def plot(self, field):
        """
        Plotting routine

        :param field: file prefix (either 'Vcy' or 'PVZ')
        :type field: str
        """
        if field == 'PVZ':
            S, Z = np.meshgrid(self.cyl_rad, self.z)
            fig = plt.figure(figsize=(4.5, 8))
            ax = fig.add_subplot(111)
            vmax = abs(self.omS).max()
            vmin = -vmax
            cs = np.linspace(vmin, vmax, 17)
            ax.contourf(S, Z, self.omS, cs, cmap=plt.get_cmap('RdYlBu'),
                        extend='both')
            theta = np.linspace(np.pi/2, -np.pi/2, 32)
            rcmb = self.z.max()
            ax.plot(rcmb*np.cos(theta), rcmb*np.sin(theta), 'k-', lw=1.5)
        elif field == 'Vcy':
            S, Z = np.meshgrid(self.cyl_rad, self.z)
            rcmb = 1./(1.-self.radratio)
            ricb = self.radratio/(1.-self.radratio)
            for field in [self.vphi]:
                dat = field.mean(axis=0)
                fig = plt.figure(figsize=(4.5, 8))
                ax = fig.add_subplot(111)
                vmax = abs(dat).max()
                vmin = -vmax
                cs = np.linspace(vmin, vmax, 17)
                ax.contourf(S, Z, dat, cs, cmap=plt.get_cmap('RdYlBu'),
                            extend='both')
                theta = np.linspace(np.pi/2, -np.pi/2, 32)
                ax.plot(rcmb*np.cos(theta), rcmb*np.sin(theta), 'k-', lw=1.5)
                ax.plot(ricb*np.cos(theta), ricb*np.sin(theta), 'k-', lw=1.5)

            phi = np.linspace(0., 2*np.pi, self.n_phi_max)
            S, pphi = np.meshgrid(self.cyl_rad, phi)
            xx = S*np.cos(pphi)
            yy = S*np.sin(pphi)
            for field in [self.vs]:
                dat = field[:, self.nZmax/2, :]
                fig = plt.figure(figsize=(6, 6))
                fig.subplots_adjust(top=0.99, right=0.99, bottom=0.01,
                                    left=0.01)
                ax = fig.add_subplot(111, frameon=False)

                vmax = abs(dat).max()
                vmin = -vmax
                cs = np.linspace(vmin, vmax, 17)
                ax.contourf(xx, yy, dat, cs, cmap=plt.get_cmap('RdYlBu'),
                            extend='both')
                theta = np.linspace(0., 2*np.pi, 32)
                ax.plot(rcmb*np.cos(theta), rcmb*np.sin(theta), 'k-', lw=1.5)
                ax.plot(ricb*np.cos(theta), ricb*np.sin(theta), 'k-', lw=1.5)
                ax.axis('off')


if __name__ == '__main__':
    MagicTOZ(tag='start', itoz=1, ave=False, verbose=True)

    MagicTOHemi(hemi='n', iplot=True)

    file = 'TO_mov.test'
    TOMovie(file=file)
    plt.show()
