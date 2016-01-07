# -*- coding: utf-8 -*-
import glob
import re
import os
import copy
import numpy as N
import matplotlib.pyplot as P
from .npfile import *
from magic.libmagic import symmetrize, hammer2cart



class Movie:
    """
    This class allows to read the :ref:`movie files <secMovieFile>` generated
    by the MagIC code.

    >>> m = Movie()
    >>> # This returns a list of the available movies in the directory
    >>> # and lets you decide which one you want to read

    >>> # Reads and display AV_mov.test
    >>> m = Movie(filed='AV_mov.test') 
    >>> print(m.data) # access to the data

    >>> # Read three movie files (no display)
    >>> m1 = Movie(file='AV_mov.testa', iplot=True)
    >>> m2 = Movie(file='AV_mov.testb', iplot=True)
    >>> m3 = Movie(file='AV_mov.testc', iplot=True)
    >>> # Stack them together
    >>> m = m1+m2+m3
    >>> # Display
    >>> m.plot(levels=33, cm='seismic', cut=0.5)

    >>> # Store the outputs in movie/img_#.png
    >>> # Only from the timesteps 280 to 380
    >>> m = Movie(file='AB_mov.test', png=True, nvar=100, lastvar=380)
    """

    def __init__(self, file=None, iplot=True, step=1, png=False,
                 lastvar=None, nvar='all', levels=12, cm='RdYlBu_r', cut=0.5,
                 bgcolor=None, fluct=False, normed=False, avg=False, 
                 std=False, dpi=80, normRad=False, precision='Float32',
                 deminc=True):
        """
        :param nvar: the number of timesteps of the movie file we want to plot
                     starting from the last line
        :type nvar: int
        :param png: if png=True, write the png files instead of display
        :type png: bool
        :param iplot: if iplot=True, display otherwise just read
        :type iplot: bool
        :param lastvar: the number of the last timesteps to be read
        :type lastvar: int
        :param step: the stepping between two timesteps
        :type step: int
        :param levels: the number of contour levels
        :type levels: int
        :param cm: the name of the color map
        :type cm: str
        :param fluct: if fluct=True, substract the axisymmetric part
        :type fluct: bool
        :param normed: the colormap is rescaled every timestep when set to True,
                       otherwise it is calculated from the global extrema
        :type normed: bool
        :param avg: if avg=True, time-average is displayed
        :type avg: bool
        :param std: if std=True, standard deviation is displayed
        :type std: bool
        :param dpi: dot per inch when saving PNGs
        :type dpi: int
        :param normRad: if normRad=True, then we normalise for each radial level
        :type normRad: bool
        :param precision: precision of the input file, Float32 for single precision,
                          Float64 for double precision
        :type precision: str
        :param cut: adjust the contour extrema to max(abs(data))*cut
        :type cut: float
        :param bgcolor: background color of the figure
        :type bgcolor: str
        :param deminc: a logical to indicate if one wants do get rid of the
                       possible azimuthal symmetry
        :type deminc: bool
        """

        if avg or std:
            iplot = False
        if file == None:
            dat = glob.glob('*[Mm]ov.*')
            str1 = 'Which movie do you want ?\n'
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
        mot = re.compile(r'.*[Mm]ov\.(.*)')
        end = mot.findall(filename)[0]

        # DETERMINE THE NUMBER OF LINES BY READING THE LOG FILE
        logfile = open('log.%s' % end, 'r')
        mot = re.compile(r'  ! WRITING MOVIE FRAME NO\s*(\d*).*')
        mot2 = re.compile(r' ! WRITING TO MOVIE FRAME NO\s*(\d*).*')
        for line in logfile.readlines():
            if mot.match(line):
                nlines = int(mot.findall(line)[0])
            elif mot2.match(line):
                nlines = int(mot2.findall(line)[0])
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
        n_type, n_surface, const, n_fields = infile.fort_read(precision)
        movtype = infile.fort_read(precision)
        n_fields = int(n_fields)
        if n_fields > 1:
            print('!!! Warning: several fields in the movie file !!!')
            print('!!! Only the last one will be displayed       !!!')
            print('!!! For TO_mov.TAG, use TOMovie(...) instead  !!!')
        self.movtype = int(movtype[0])
        n_surface = int(n_surface)

        # RUN PARAMETERS
        runid = infile.fort_read('|S64')
        n_r_mov_tot, n_r_max, n_theta_max, n_phi_tot, self.minc, self.ra, \
             self.ek, self.pr, self.prmag, self.radratio, self.tScale = infile.fort_read(precision)
        self.minc = int(self.minc)
        n_r_mov_tot = int(n_r_mov_tot)
        self.n_r_max = int(n_r_max)
        self.n_theta_max = int(n_theta_max)
        self.n_phi_tot = int(n_phi_tot)

        # GRID
        self.radius = infile.fort_read(precision)
        self.radius = self.radius[:self.n_r_max] # remove inner core
        self.theta = infile.fort_read(precision)
        self.phi = infile.fort_read(precision)

        if n_surface == 0:
            self.surftype = '3d volume'
            shape = (n_r_mov_tot+2, self.n_theta_max, self.n_phi_tot)
        elif n_surface == 1:
            self.surftype = 'r_constant'
            shape = (self.n_theta_max, self.n_phi_tot)
            self.data = N.zeros((self.nvar, self.n_phi_tot, self.n_theta_max), precision)
        elif n_surface == 2:
            self.surftype = 'theta_constant'
            if self.movtype in [1, 2, 3]: # read inner core
                shape = (n_r_mov_tot+2, self.n_phi_tot)
            else:
                shape = (self.n_r_max, self.n_phi_tot)
            self.data = N.zeros((self.nvar, self.n_phi_tot, self.n_r_max), precision)
        elif n_surface == 3:
            self.surftype = 'phi_constant'
            if self.movtype in [1, 2, 3]: # read inner core
                shape = (n_r_mov_tot+2, 2*self.n_theta_max)
            elif self.movtype in [8, 9]:
                shape = (n_r_mov_tot+2, self.n_theta_max)
            elif self.movtype in [4, 5, 6, 7, 16, 17, 18, 47, 54]:
                shape = (self.n_r_max, 2*self.n_theta_max)
            elif self.movtype in [10, 11, 12, 19, 92]:
                shape = (self.n_r_max, self.n_theta_max)
            # Inner core is not stored here
            self.data = N.zeros((self.nvar, self.n_theta_max, self.n_r_max), precision)

        self.time = N.zeros(self.nvar, precision)

        # READ the data

        # If one skip the beginning, nevertheless read but do not store
        for i in range(self.var2-self.nvar):
            n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                                   movieDipLon, movieDipStrength, \
                            movieDipStrengthGeo = infile.fort_read(precision)
            for ll in range(n_fields):
                dat = infile.fort_read(precision, shape=shape)
        # then read the remaining requested nvar lines
        for k in range(self.nvar):
            n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                                   movieDipLon, movieDipStrength, \
                            movieDipStrengthGeo = infile.fort_read(precision)
            self.time[k] = t_movieS
            for ll in range(n_fields):
                dat = infile.fort_read(precision, shape=shape)
                if n_surface == 2:
                    if self.movtype in [1, 2, 3]:
                        dat = dat[:self.n_r_max, :].T
                        self.data[k, ...] = dat
                    else:
                        self.data[k, ...] = dat.T
                elif n_surface == 3:
                    if self.movtype in [1, 2, 3]:
                        dat = dat[:self.n_r_max, :self.n_theta_max].T
                        self.data[k, :, ::2] = dat[:, :n_r_max/2+1]
                        self.data[k, :, 1::2] = dat[:, n_r_max/2+1:]
                    elif self.movtype in [8, 9]:
                        dat = dat[:self.n_r_max, :].T
                        self.data[k, ...] = dat
                    elif self.movtype in [4, 5, 6, 7, 16, 17, 18, 47, 54, 91]:
                        dat = dat[:, :self.n_theta_max].T
                        self.data[k, :, ::2] = dat[:, :n_r_max/2+1]
                        self.data[k, :, 1::2] = dat[:, n_r_max/2+1:]
                    elif self.movtype in [10, 11, 12, 19, 92]:
                        self.data[k, ...] = dat.T
                else:
                    self.data[k, ...] = dat.T
            if fluct:
                self.data[k, ...] = self.data[k, ...]-self.data[k, ...].mean(axis=0)

        infile.close()

        if normRad:
            norm = N.sqrt(N.mean(self.data**2, axis=1))
            norm = norm.mean(axis=0)
            self.data[:, :, norm!=0.] /= norm[norm!=0.]

        if iplot:
            cmap = P.get_cmap(cm)
            self.plot(cut, levels, cmap, png, step, normed, dpi, bgcolor, deminc)
        if avg or std:
            cmap = P.get_cmap(cm)
            self.avgStd(std, cut, levels, cmap)

    def __add__(self, new):
        """
        Built-in function to sum two movies

        .. note:: So far this function only works for two movies with the same 
                  grid sizes. At some point, we might introduce grid extrapolation 
                  to allow any summation/
        """
        out = copy.deepcopy(new)
        out.time = N.concatenate((self.time, new.time), axis=0)
        out.data = N.concatenate((self.data, new.data), axis=0)
        out.nvar = self.nvar+new.nvar
        out.var2 = out.nvar
        return out

    def avgStd(self, std=False, cut=0.5, levels=12, cmap='RdYlBu_r'):
        """
        Plot time-average or standard deviation

        :param std: the standard deviation is computed instead the average
                    when std is True
        :type std: bool
        :param levels: number of contour levels
        :type levels: int
        :param cmap: name of the colormap
        :type cmap: str
        :param cut: adjust the contour extrema to max(abs(data))*cut
        :type cut: float
        """
        if std:
            avg = self.data.std(axis=0)
        else:
            avg = self.data.mean(axis=0)
        vmin = - max(abs(avg.max()), abs(avg.min()))
        vmin = cut * vmin
        vmax = -vmin
        cs = N.linspace(vmin, vmax, levels)

        if self.surftype == 'phi_constant':
            th = N.linspace(N.pi/2., -N.pi/2., self.n_theta_max)
            rr, tth = N.meshgrid(self.radius, th)
            xx = rr * N.cos(tth)
            yy = rr * N.sin(tth)
            xxout = rr.max() * N.cos(th)
            yyout = rr.max() * N.sin(th)
            xxin = rr.min() * N.cos(th)
            yyin = rr.min() * N.sin(th)
            fig = P.figure(figsize=(4, 8))
        elif self.surftype == 'r_constant':
            th = N.linspace(N.pi/2., -N.pi/2., self.n_theta_max)
            phi = N.linspace(-N.pi, N.pi, self.n_phi_tot)
            ttheta, pphi = N.meshgrid(th, phi)
            xx, yy = hammer2cart(ttheta, pphi)
            xxout, yyout = hammer2cart(th, -N.pi)
            xxin, yyin = hammer2cart(th, N.pi)
            fig = P.figure(figsize=(8, 4))
        elif self.surftype == 'theta_constant':
            phi = N.linspace(0., 2.*N.pi, self.n_phi_tot)
            rr, pphi = N.meshgrid(self.radius, phi)
            xx = rr * N.cos(pphi)
            yy = rr * N.sin(pphi)
            xxout = rr.max() * N.cos(pphi)
            yyout = rr.max() * N.sin(pphi)
            xxin = rr.min() * N.cos(pphi)
            yyin = rr.min() * N.sin(pphi)
            fig = P.figure(figsize=(6, 6))
        elif self.surftype == '3d volume':
            self.data = self.data[..., 0]
            th = N.linspace(N.pi/2., -N.pi/2., self.n_theta_max)
            phi = N.linspace(-N.pi, N.pi, self.n_phi_tot)
            ttheta, pphi = N.meshgrid(th, phi)
            xx, yy = hammer2cart(ttheta, pphi)
            xxout, yyout = hammer2cart(th, -N.pi)
            xxin, yyin = hammer2cart(th, N.pi)
            fig = P.figure(figsize=(8, 4))

        fig.subplots_adjust(top=0.99, right=0.99, bottom=0.01, left=0.01)
        ax = fig.add_subplot(111, frameon=False)
        im = ax.contourf(xx, yy, avg, cs, cmap=cmap, extend='both')
        ax.plot(xxout, yyout, 'k-', lw=1.5)
        ax.plot(xxin, yyin, 'k-', lw=1.5)
        ax.axis('off')

    def plot(self, cut=0.5, levels=12, cmap='RdYlBu_r', png=False, step=1, 
             normed=False, dpi=80, bgcolor=None, deminc=True):
        """
        Plotting function (it can also write the png files)

        :param levels: number of contour levels
        :type levels: int
        :param cmap: name of the colormap
        :type cmap: str
        :param cut: adjust the contour extrema to max(abs(data))*cut
        :type cut: float
        :param png: save the movie as a series of png files when
                    set to True
        :type png: bool
        :param dpi: dot per inch when saving PNGs
        :type dpi: int
        :param bgcolor: background color of the figure
        :type bgcolor: str
        :param normed: the colormap is rescaled every timestep when set to True,
                       otherwise it is calculated from the global extrema
        :type normed: bool
        :param step: the stepping between two timesteps
        :type step: int
        :param deminc: a logical to indicate if one wants do get rid of the
                       possible azimuthal symmetry
        :type deminc: bool
        """

        if png:
            P.ioff()
            if not os.path.exists('movie'):
                os.mkdir('movie')
        else:
            P.ion()

        if not normed:
            vmin = - max(abs(self.data.max()), abs(self.data.min()))
            vmin = cut * vmin
            vmax = -vmin
            #vmin, vmax = self.data.min(), self.data.max()
            cs = N.linspace(vmin, vmax, levels)

        if self.surftype == 'phi_constant':
            #if self.movtype in [1, 7]:
                #th = N.linspace(0., 2.*N.pi, 2*self.n_theta_max)
            #else:
            th = N.linspace(N.pi/2., -N.pi/2., self.n_theta_max)
            rr, tth = N.meshgrid(self.radius, th)
            xx = rr * N.cos(tth)
            yy = rr * N.sin(tth)
            xxout = rr.max() * N.cos(th)
            yyout = rr.max() * N.sin(th)
            xxin = rr.min() * N.cos(th)
            yyin = rr.min() * N.sin(th)
            fig = P.figure(figsize=(4, 8))
        elif self.surftype == 'r_constant':
            th = N.linspace(N.pi/2., -N.pi/2., self.n_theta_max)
            if deminc:
                phi = N.linspace(-N.pi, N.pi, self.n_phi_tot*self.minc+1)
                xxout, yyout = hammer2cart(th, -N.pi)
                xxin, yyin = hammer2cart(th, N.pi)
            else:
                phi = N.linspace(-N.pi/self.minc, N.pi/self.minc, self.n_phi_tot)
                xxout, yyout = hammer2cart(th, -N.pi/self.minc)
                xxin, yyin = hammer2cart(th, N.pi/self.minc)
            ttheta, pphi = N.meshgrid(th, phi)
            xx, yy = hammer2cart(ttheta, pphi)
            fig = P.figure(figsize=(8, 4))
        elif self.surftype == 'theta_constant':
            if deminc:
                phi = N.linspace(0., 2.*N.pi, self.n_phi_tot*self.minc+1)
            else:
                phi = N.linspace(0., 2.*N.pi/self.minc, self.n_phi_tot)
            rr, pphi = N.meshgrid(self.radius, phi)
            xx = rr * N.cos(pphi)
            yy = rr * N.sin(pphi)
            xxout = rr.max() * N.cos(pphi)
            yyout = rr.max() * N.sin(pphi)
            xxin = rr.min() * N.cos(pphi)
            yyin = rr.min() * N.sin(pphi)
            fig = P.figure(figsize=(6, 6))
        elif self.surftype == '3d volume':
            self.data = self.data[..., 0]
            th = N.linspace(N.pi/2., -N.pi/2., self.n_theta_max)
            phi = N.linspace(-N.pi, N.pi, self.n_phi_tot)
            ttheta, pphi = N.meshgrid(th, phi)
            xx, yy = hammer2cart(ttheta, pphi)
            xxout, yyout = hammer2cart(th, -N.pi)
            xxin, yyin = hammer2cart(th, N.pi)
            fig = P.figure(figsize=(8, 4))

        fig.subplots_adjust(top=0.99, right=0.99, bottom=0.01, left=0.01)
        ax = fig.add_subplot(111, frameon=False)

        for k in range(self.nvar):
            if k == 0:
                if normed:
                    vmin = - max(abs(self.data[k, ...].max()), abs(self.data[k, ...].min()))
                    vmin = cut * vmin
                    vmax = -vmin
                    cs = N.linspace(vmin, vmax, levels)
                if self.surftype in ['r_constant', 'theta_constant']:
                    if deminc:
                        dat = symmetrize(self.data[k, ...], self.minc)
                    else:
                        dat = self.data[k, ...]
                    im = ax.contourf(xx, yy, dat, cs, cmap=cmap, extend='both')
                else:
                    im = ax.contourf(xx, yy, self.data[k, ...], cs, cmap=cmap, extend='both')
                ax.plot(xxout, yyout, 'k-', lw=1.5)
                ax.plot(xxin, yyin, 'k-', lw=1.5)
                ax.axis('off')
                man = P.get_current_fig_manager()
                man.canvas.draw()
            if k != 0 and k % step == 0:
                if not png:
                    print(k+self.var2-self.nvar)
                P.cla()
                if normed:
                    vmin = - max(abs(self.data[k, ...].max()), abs(self.data[k, ...].min()))
                    vmin = cut * vmin
                    vmax = -vmin
                    cs = N.linspace(vmin, vmax, levels)
                if self.surftype in ['r_constant', 'theta_constant']:
                    if deminc:
                        dat = symmetrize(self.data[k, ...], self.minc)
                    else:
                        dat = self.data[k, ...]
                    im = ax.contourf(xx, yy, dat, cs, cmap=cmap, extend='both')
                else:
                    im = ax.contourf(xx, yy, self.data[k, ...], cs, cmap=cmap, extend='both')
                ax.plot(xxout, yyout, 'k-', lw=1.5)
                ax.plot(xxin, yyin, 'k-', lw=1.5)
                ax.axis('off')
                man.canvas.draw()
            if png:
                filename = 'movie/img%05d.png' % k
                print('write %s' % filename)
                #st = 'echo %i' % ivar + ' > movie/imgmax'
                if bgcolor is not None:
                    fig.savefig(filename, facecolor=bgcolor, dpi=dpi)
                else:
                    fig.savefig(filename, dpi=dpi)



if __name__ == '__main__':
    Movie(step=1)
    P.show()
