# -*- coding: utf-8 -*-
import glob
import re
import os
import copy
import numpy as N
import matplotlib.pyplot as P
import scipy.interpolate as S
from magic.libmagic import cut, hammer2cart
from magic.setup import labTex
from scipy.integrate import simps
from .npfile import *


class Butterfly:

    def __init__(self, file=None, avg=False, step=1, iplot=True, rad=0.8,
                 lastvar=None, nvar='all', levels=20, cm='RdYlBu_r',
                 precision='Float32'):
        """
        :param nvar: the number of lines of the movie file we want to plot
                     starting from the last line
        :param lastvar: the rank of the last line to be read
        :param step: the stepping between two lines             
        :param levels: the number of contour levels
        :param cm: the name of the color map
        :param png: if png=True, write the png outputs
        :param precision: precision of the input file, Float32 for single precision,
                          Float64 for double precision
        """

        self.precision = precision

        if file == None:
            dat = glob.glob('*_mov.*')
            str1 = 'Which movie do you want ?\n'
            for k, movie in enumerate(dat):
                str1 += ' %i) %s\n' % (k+1, movie)
            index = input(str1)
            try:
                filename = dat[index-1]
            except IndexError:
                print('Non valid index: %s has been chosen instead' % dat[0])
                filename = dat[0]

        else:
            filename = file
        mot = re.compile(r'.*_mov\.(.*)')
        end = mot.findall(filename)[0]

        # DETERMINE THE NUMBER OF LINES BY READING THE LOG FILE
        logfile = open('log.%s' % end, 'r')
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
        ind = N.nonzero(N.where(abs(self.radius-self.rad) \
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
            self.data = N.zeros((self.n_theta_max, self.nvar), self.precision)
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
            self.data = N.zeros((self.n_theta_max, self.nvar), self.precision)

        self.cmap = P.get_cmap(cm)
        self.time = N.zeros(self.nvar, self.precision)

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
                    #self.data[:, k] = data.mean(axis=1)
                    self.data[:, k] = data[:, self.n_phi_tot/2.]
                elif self.surftype == 'phi_constant':
                    data = infile.fort_read(self.precision, shape=shape)
                    self.data[:, k] = data[self.indPlot, :]

            else: # Nevertheless read
                data = infile.fort_read(self.precision, shape=shape)


        if step != 1:
            self.time = self.time[::step]
            self.data = self.data[:, ::step]

        if iplot:
            self.plot()

    def __add__(self, new):
        out = copy.deepcopy(new)
        out.time = N.concatenate((self.time, new.time), axis=0)
        out.data = N.concatenate((self.data, new.data), axis=1)
        return out

    def plot(self, levels=12, contour=False, renorm=False, cut=0.5, mesh=3):
        """
        :param levels: the number of contour levels (only used in contour=True)
        :param contour: contour levels instead of image
        :param renorm: rebin the time series in case of irregularly spaced data
        :param cut: levels cutting
        """
        fig = P.figure()
        ax = fig.add_subplot(111)
        vmax = max(abs(self.data.max()), abs(self.data.min()))
        if contour:
            im = ax.contourf(self.time, self.theta*180/N.pi, self.data[::-1,:], 
                            levels, cmap=self.cmap, aa=True)
        else:
            extent = self.time.min(), self.time.max(), -90, 90
            if renorm:
                nx = mesh*len(self.time)
                x = N.linspace(self.time.min(), self.time.max(), nx)
                data = N.zeros((len(self.theta), nx), self.precision)
                for i in range(self.data.shape[0]):
                    tckp = S.splrep(self.time, self.data[i, :])
                    data[i, :] = S.splev(x, tckp)
                im = ax.imshow(data, origin='upper', aspect='auto', 
                           cmap=self.cmap, extent=extent, interpolation='bicubic')
                #im = ax.contourf(zi[:-1, :],  cmap=self.cmap)
                #vmax = max(abs(zi.max()), abs(zi.min()))
            else:
                im = ax.imshow(self.data, origin='upper', aspect='auto', 
                           cmap=self.cmap, extent=extent, interpolation='bicubic')
        im.set_clim(-cut*vmax, cut*vmax)
        ax.set_xlabel('Time')
        ax.set_ylabel('Latitude')
        if self.surftype == 'phi_constant':
            if labTex:
                txt = r'$r = %.1f\ r_o$' % self.rad
            else:
                txt = r'r = %.1f ro' % self.rad
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
        :param renorm: rebin the time series in case of irregularly spaced data
        """
        if renorm:
            nx = 3*len(self.time)
            x = N.linspace(self.time.min(), self.time.max(), nx)
            data = N.zeros((len(self.theta), nx), self.precision)
            for i in range(self.data.shape[0]):
                tckp = S.splrep(self.time, self.data[i, :])
                data[i, :] = S.splev(x, tckp)
            self.data = data
            self.time = x
            print("####################")
            print("Warning time and data have been replaced by the extrapolated values !!!")
            print("####################")
        nt = self.time.shape[0]
        w1 = N.fft.fft(self.data, axis=1)
        self.amp = N.abs(w1[:, 1:nt/2+1])
        dw = 2.*N.pi/(self.time[-1]-self.time[0])
        w = dw*N.arange(nt)
        self.omega = w[1:nt/2+1]
        self.amp1D = N.zeros_like(self.omega)
        for i in range(len(self.omega)):
            self.amp1D[i] = simps(self.amp[:, i], self.theta)

        fig = P.figure()
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
        print('Fourier frequency:%.2f' % (self.omega[self.amp1D==self.amp1D.max()]))


if __name__ == '__main__':
    t1 = Butterfly(file='Br_CMB_mov.ccondAnelN3MagRa2e7Pm2ggg', step=1, 
                   iplot=False)
    t1.plot()
    P.show()

