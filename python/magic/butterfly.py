# -*- coding: utf-8 -*-
import glob
import re
import os
import copy
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as S
from magic.plotlib import cut
from magic.setup import labTex
from scipy.integrate import simps
from .npfile import *


class Butterfly:
    """
    This class can be used to display the time evolution of the magnetic field for
    various latitudes (i.e. the well-known butterfly diagrams). These diagrams are
    usually constructed using MagIC's :ref:`movie files <secMovieFile>`: either
    radial cuts (like Br_CMB_mov.TAG) or azimuthal-average (like AB_mov.TAG).

    >>> # Read Br_CMB_mov.ccondAnelN3MagRa2e7Pm2ggg
    >>> t1 = Butterfly(file='Br_CMB_mov.ccondAnelN3MagRa2e7Pm2ggg', step=1,
                       iplot=False)
    >>> # Plot it
    >>> t1.plot(levels=33, cm='seismic', cut=0.6)
    """

    def __init__(self, file=None, step=1, iplot=True, rad=0.8,
                 lastvar=None, nvar='all', levels=20, cm='RdYlBu_r',
                 precision=np.float32, cut=0.8):
        """
        :param file: when specified, the constructor reads this file, otherwise
                     a list with the possible options is displayed
        :type file: str
        :param rad: radial level (normalised to the outer boundary radius)
        :type rad: float
        :param iplot: display/hide the plots (default is True)
        :type iplot: bool
        :param nvar: the number of time steps (lines) of the movie file one
                     wants to plot starting from the last line
        :type nvar: int
        :param lastvar: the number of the last time step to be read
        :type lastvar: int
        :param step: the stepping between two lines
        :type step: int
        :param levels: the number of contour levels
        :type levels: int
        :param cm: the name of the color map
        :type cm: str
        :param cut: adjust the contour extrema to max(abs(data))*cut
        :type cut: float
        :param precision: precision of the input file, np.float32 for single
                          precision, np.float64 for double precision
        :type precision: bool
        """

        self.precision = precision

        if file is None:
            dat = glob.glob('*_mov.*')
            str1 = 'Which movie do you want ?\n'
            for k, movie in enumerate(dat):
                str1 += ' {}) {}\n'.format(k+1, movie)
            index = input(str1)
            try:
                filename = dat[int(index)-1]
            except IndexError:
                print('Non valid index: {} has been chosen instead'.format(dat[0]))
                filename = dat[0]

        else:
            filename = file
        mot = re.compile(r'.*_mov\.(.*)')
        end = mot.findall(filename)[0]

        # DETERMINE THE NUMBER OF LINES BY READING THE LOG FILE
        logfile = open('log.{}'.format(end), 'r')
        mot = re.compile(r'  ! WRITING MOVIE FRAME NO\s*(\d*).*')
        for line in logfile.readlines():
            if mot.match(line):
                 nlines = int(mot.findall(line)[0])
        logfile.close()
        if lastvar is None:
            self.var2 = nlines
        else:
            self.var2 = lastvar
        if str(nvar) == 'all':
            self.nvar = nlines
            self.var2 = nlines
        else:
            self.nvar = nvar

        # READ the movie file
        infile = npfile(filename, endian='B')
        # HEADER
        version = infile.fort_read('|S64')
        n_type, n_surface, const, n_fields = infile.fort_read(self.precision)
        movtype = infile.fort_read(self.precision)
        self.movtype = int(movtype)
        n_surface = int(n_surface)

        # RUN PARAMETERS
        runid = infile.fort_read('|S64')
        n_r_mov_tot, n_r_max, n_theta_max, n_phi_tot, minc, self.ra, \
             self.ek, self.pr, self.prmag, self.radratio, tScale = infile.fort_read(self.precision)
        self.n_r_max = int(n_r_max)
        self.n_theta_max = int(n_theta_max)
        self.n_phi_tot = int(n_phi_tot)

        # GRID
        self.rad = rad
        self.radius = infile.fort_read(self.precision)
        ind = np.nonzero(np.where(abs(self.radius-self.rad) \
                        == min(abs(self.radius-self.rad)), 1, 0))
        self.indPlot = ind[0][0]
        self.theta = infile.fort_read(self.precision)
        self.phi = infile.fort_read(self.precision)

        if n_surface == 0:
            self.surftype = '3d volume'
            shape = (self.n_r_max, self.n_theta_max, self.n_phi_tot)
        elif n_surface == 1:
            self.surftype = 'r_constant'
            shape = (self.n_theta_max, self.n_phi_tot)
            self.data = np.zeros((self.n_theta_max, self.nvar), self.precision)
        elif n_surface == 2:
            self.surftype = 'theta_constant'
            shape = (self.n_r_max, self.n_phi_tot)
        elif n_surface == 3:
            self.surftype = 'phi_constant'
            shape = (self.n_r_max, self.n_theta_max)
            if self.movtype in [8, 9]:
                shape = (int(n_r_mov_tot)+2, self.n_theta_max)
            else:
                shape = (self.n_r_max, self.n_theta_max)
            self.data = np.zeros((self.n_theta_max, self.nvar), self.precision)

        self.cmap = plt.get_cmap(cm)
        self.time = np.zeros(self.nvar, self.precision)

        for i in range(self.var2-self.nvar):
            n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                 movieDipLon, movieDipStrength, \
                 movieDipStrengthGeo = infile.fort_read(self.precision)
            data = infile.fort_read(self.precision, shape=shape)
        for k in range(self.nvar):
            n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                 movieDipLon, movieDipStrength, \
                 movieDipStrengthGeo = infile.fort_read(self.precision)
            if k % step == 0:
                self.time[k] = t_movieS
                #print(k+self.var2-self.nvar)
                if self.surftype == 'r_constant':
                    data = infile.fort_read(self.precision, shape=shape)
                    self.data[:, k] = data.mean(axis=1)
                elif self.surftype == 'phi_constant':
                    data = infile.fort_read(self.precision, shape=shape)
                    self.data[:, k] = data[self.indPlot, :]

            else: # Nevertheless read
                data = infile.fort_read(self.precision, shape=shape)


        if step != 1:
            self.time = self.time[::step]
            self.data = self.data[:, ::step]

        if iplot:
            self.plot(levels=levels,cm=self.cmap,cut=cut)

    def __add__(self, new):
        """
        Overload of the addition operator

        >>> # Read 2 files
        >>> b1 = Butterfly(file='AB_mov.test1', iplot=False)
        >>> b2 = Butterfly(file='AB_mov.test2', iplot=False)
        >>> # Stack them and display the whole thing
        >>> b = b1+b2
        >>> b.plot(levels=33, contour=True, cut=0.8, cm='seismic')
        """
        out = copy.deepcopy(new)
        out.time = np.concatenate((self.time, new.time), axis=0)
        out.data = np.concatenate((self.data, new.data), axis=1)
        return out

    def plot(self, levels=12, contour=False, renorm=False, cut=0.5, mesh=3,
             cm='RdYlBu_R'):
        """
        Plotting function

        :param cm: name of the colormap
        :type cm: str
        :param levels: the number of contour levels (only used when iplot=True and
                       contour=True)
        :type levels: int
        :param contour: when set to True, display contour levels (pylab.contourf),
                        when set to False, display an image (pylab.imshow)
        :type contour: bool
        :param renorm: when set to True, it re-bins the time series in case of
                       irregularly time-spaced data
        :type renorm: bool
        :param mesh: when renorm=True, factor of regriding: NewTime = mesh*OldTime
        :type mesh: int
        :param cut: adjust the contour extrema to max(abs(data))*cut
        :type cut: float
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        vmax = max(abs(self.data.max()), abs(self.data.min()))
        if contour:
            im = ax.contourf(self.time, self.theta*180/np.pi, self.data[::-1,:],
                            levels, cmap=cm, aa=True)
        else:
            extent = self.time.min(), self.time.max(), -90, 90
            if renorm:
                nx = mesh*len(self.time)
                x = np.linspace(self.time.min(), self.time.max(), nx)
                data = np.zeros((len(self.theta), nx), self.precision)
                for i in range(self.data.shape[0]):
                    tckp = S.splrep(self.time, self.data[i, :])
                    data[i, :] = S.splev(x, tckp)
                im = ax.imshow(data, origin='upper', aspect='auto',
                           cmap=cm, extent=extent, interpolation='bicubic')
                #im = ax.contourf(zi[:-1, :],  cmap=cm)
                #vmax = max(abs(zi.max()), abs(zi.min()))
            else:
                im = ax.imshow(self.data, origin='upper', aspect='auto',
                           cmap=cm, extent=extent, interpolation='bicubic')
        im.set_clim(-cut*vmax, cut*vmax)
        ax.set_xlabel('Time')
        ax.set_ylabel('Latitude')
        if self.surftype == 'phi_constant':
            if labTex:
                txt = r'$r = {:.1f}\ r_o$'.format(self.rad)
            else:
                txt = r'r = {:.1f} ro'.format(self.rad)
            ax.text(0.8, 0.07, txt, transform =ax.transAxes)
            #ax.text(0.05, 0.9, 'a)', transform =ax.transAxes)
        ax.set_yticks([-90,-45,0,45,90])
        if labTex:
            ax.set_yticklabels(['$-90^\circ$', '$-45^\circ$', '$0^\circ$',
                                '$45^\circ$', '$90^\circ$'])
        else:
            ax.set_yticklabels(['-90', '-45', '0', '45', '90'])

    def fourier2D(self, renorm=False):
        """
        This function allows to conduct some basic Fourier analysis on the
        data. It displays two figures: the first one is a contour levels
        in the (Frequency, Latitude) plane, the second one is integrated over
        latitudes (thus a simple, power vs Frequency plot)

        >>> # Load the data without plotting
        >>> b1 = Butterfly(file='AB_mov.test1', iplot=False)
        >>> # Fourier analysis
        >>> b1.fourier2D()

        :param renorm: when set to True, it rebins the time series in case of
                       irregularly spaced data
        :type renorm: bool
        """
        if renorm:
            nx = 3*len(self.time)
            x = np.linspace(self.time.min(), self.time.max(), nx)
            data = np.zeros((len(self.theta), nx), self.precision)
            for i in range(self.data.shape[0]):
                tckp = S.splrep(self.time, self.data[i, :])
                data[i, :] = S.splev(x, tckp)
            self.data = data
            self.time = x
            print("####################")
            print("Warning time and data have been replaced by the extrapolated values !!!")
            print("####################")
        nt = self.time.shape[0]
        w1 = np.fft.fft(self.data, axis=1)
        self.amp = np.abs(w1[:, 1:nt//2+1])
        dw = 2.*np.pi/(self.time[-1]-self.time[0])
        w = dw*np.arange(nt)
        self.omega = w[1:nt//2+1]
        self.amp1D = np.zeros_like(self.omega)
        for i in range(len(self.omega)):
            self.amp1D[i] = simps(self.amp[:, i], self.theta)

        fig = plt.figure()
        ax = fig.add_subplot(211)
        extent = self.omega.min(), self.omega.max(), -90, 90
        ax.imshow(self.amp, extent=extent, aspect='auto', origin='upper')
        ax.set_xlabel('Frequency')
        ax.set_ylabel('Latitude')

        ax = fig.add_subplot(212)
        ax.semilogy(self.omega, self.amp1D)
        ax.set_xlabel('Frequency')
        ax.set_ylabel('Spectrum')
        ax.set_xlim(self.omega.min(), self.omega.max())
        print('Fourier frequency:{:.2f}'.format(self.omega[self.amp1D==self.amp1D.max()]))


if __name__ == '__main__':
    t1 = Butterfly(file='Br_CMB_mov.ccondAnelN3MagRa2e7Pm2ggg', step=1,
                   iplot=False)
    t1.plot()
    plt.show()
