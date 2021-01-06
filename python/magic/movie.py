# -*- coding: utf-8 -*-
import glob
import re
import os
import copy
import numpy as np
import matplotlib.pyplot as plt
from .npfile import *
from magic.libmagic import symmetrize, chebgrid
from magic.plotlib import hammer2cart

def getNlines(file_name, endian='B', precision=np.float32):
    """
    This function determines the number of lines of a binary file.

    :param file_name: name of the input file
    :type file_name: str
    :param endian: endianness of the file ('B' or 'l')
    :type endian: str
    :param precision: precision of the data contained in the input file
                      (np.float32 or np.float64)
    :type endian: str
    :returns: the number of lines
    :rtype: int
    """

    f = npfile(file_name, endian=endian)
    n_lines = 0
    while 1:
        try:
            f.fort_read(precision)
            n_lines += 1
        except TypeError:
            break

    f.close()
    return n_lines


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
    >>> m1 = Movie(file='AV_mov.testa', iplot=False)
    >>> m2 = Movie(file='AV_mov.testb', iplot=False)
    >>> m3 = Movie(file='AV_mov.testc', iplot=False)
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
                 std=False, dpi=80, normRad=False, precision=np.float32,
                 deminc=True, ifield=0):
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
        :param precision: precision of the input file, np.float32 for single
                          precision, np.float64 for double precision
        :type precision: str
        :param cut: adjust the contour extrema to max(abs(data))*cut
        :type cut: float
        :param bgcolor: background color of the figure
        :type bgcolor: str
        :param deminc: a logical to indicate if one wants do get rid of the
                       possible azimuthal symmetry
        :type deminc: bool
        :param ifield: in case of a multiple-field movie file, you can change
                       the default field displayed using the parameter ifield
        :type ifield: int
        """

        if avg or std:
            iplot = False
        if file is None:
            dat = glob.glob('*[Mm]ov.*')
            str1 = 'Which movie do you want ?\n'
            for k, movie in enumerate(dat):
                str1 += ' {}) {}\n'.format(k+1, movie)
            index = int(input(str1))
            try:
                filename = dat[index-1]
            except IndexError:
                print('Non valid index: {} has been chosen instead'.format(dat[0]))
                filename = dat[0]

        else:
            filename = file
        mot = re.compile(r'.*[Mm]ov\.(.*)')
        end = mot.findall(filename)[0]

        # Read the movie file
        infile = npfile(filename, endian='B')
        # Header
        version = infile.fort_read('|S64')
        n_type, n_surface, const, n_fields = infile.fort_read(precision)
        movtype = infile.fort_read(precision)
        self.n_fields = int(n_fields)
        if self.n_fields > 1:
            print('!!! Warning: several fields in the movie file !!!')
            print('!!! {} fields !!!'.format(self.n_fields))
            print('!!! The one displayed is controlled by the    !!!')
            print('!!! input parameter ifield (=0 by default)    !!!')
        self.movtype = int(movtype[ifield])
        n_surface = int(n_surface)

        # Run parameters
        runid = infile.fort_read('|S64')
        n_r_mov_tot, n_r_max, n_theta_max, n_phi_tot, self.minc, self.ra, \
            self.ek, self.pr, self.prmag, self.radratio, self.tScale =   \
            infile.fort_read(precision)
        self.minc = int(self.minc)
        n_r_mov_tot = int(n_r_mov_tot)
        self.n_r_max = int(n_r_max)
        self.n_r_ic_max = n_r_mov_tot-self.n_r_max
        self.n_theta_max = int(n_theta_max)
        self.n_phi_tot = int(n_phi_tot)

        # Grid
        self.radius = infile.fort_read(precision)
        self.radius_ic = np.zeros((self.n_r_ic_max+2), precision)
        self.radius_ic[:-1] = self.radius[self.n_r_max-1:]

        self.radius = self.radius[:self.n_r_max]  # remove inner core
        # Overwrite radius to ensure double-precision of the
        # grid (useful for Cheb der)
        rout = 1./(1.-self.radratio)
        rin = self.radratio/(1.-self.radratio)
        self.radius *= rout
        self.radius_ic *= rout
        # self.radius = chebgrid(self.n_r_max-1, rout, rin)
        self.theta = infile.fort_read(precision)
        self.phi = infile.fort_read(precision)

        # Determine the number of lines by reading the log.TAG file
        logfile = open('log.{}'.format(end), 'r')
        mot = re.compile(r'  ! WRITING MOVIE FRAME NO\s*(\d*).*')
        mot2 = re.compile(r' ! WRITING TO MOVIE FRAME NO\s*(\d*).*')
        nlines = 0
        for line in logfile.readlines():
            if mot.match(line):
                nlines = int(mot.findall(line)[0])
            elif mot2.match(line):
                nlines = int(mot2.findall(line)[0])
        logfile.close()

        # In case no 'nlines' can be determined from the log file:
        if nlines == 0:
            nlines = getNlines(filename, endian='B', precision=precision)
            nlines -= 8  # Remove 8 lines of header
            nlines /= (self.n_fields+1)

        if lastvar is None:
            self.var2 = nlines
        else:
            self.var2 = lastvar
        if str(nvar) == 'all':
            self.nvar = nlines
            self.var2 = nlines
        else:
            self.nvar = nvar

        if n_surface == 0:
            self.surftype = '3d volume'
            if self.movtype in [1, 2, 3]:
                shape = (n_r_mov_tot+2, self.n_theta_max, self.n_phi_tot)
            else:
                shape = (self.n_r_max, self.n_theta_max, self.n_phi_tot)
            self.data = np.zeros((self.n_fields, self.nvar, self.n_phi_tot,
                                  self.n_theta_max, self.n_r_max), precision)
            self.data_ic = np.zeros((self.n_fields, self.nvar, self.n_phi_tot,
                                    self.n_theta_max, self.n_r_ic_max+2),
                                    precision)
        elif n_surface == 1 or n_surface == -1:
            self.surftype = 'r_constant'
            shape = (self.n_theta_max, self.n_phi_tot)
            self.data = np.zeros((self.n_fields, self.nvar, self.n_phi_tot,
                                  self.n_theta_max), precision)
        elif n_surface == 2:
            self.surftype = 'theta_constant'
            if self.movtype in [1, 2, 3, 14]:  # read inner core
                shape = (n_r_mov_tot+2, self.n_phi_tot)
            else:
                shape = (self.n_r_max, self.n_phi_tot)
            self.data = np.zeros((self.n_fields, self.nvar, self.n_phi_tot,
                                  self.n_r_max), precision)
            self.data_ic = np.zeros((self.n_fields, self.nvar, self.n_phi_tot,
                                    self.n_r_ic_max+2), precision)
        elif n_surface == 3:
            self.surftype = 'phi_constant'
            if self.movtype in [1, 2, 3, 14]:  # read inner core
                shape = (n_r_mov_tot+2, self.n_theta_max)
                self.n_theta_plot = 2*self.n_theta_max
            elif self.movtype in [8, 9]:
                shape = (n_r_mov_tot+2, self.n_theta_max)
                self.n_theta_plot = self.n_theta_max
            elif self.movtype in [4, 5, 6, 7, 16, 17, 18, 47, 54, 109]:
                shape = (self.n_r_max, self.n_theta_max)
                self.n_theta_plot = 2*self.n_theta_max
            elif self.movtype in [10, 11, 12, 19, 92, 94, 95, 110]:
                shape = (self.n_r_max, self.n_theta_max)
                self.n_theta_plot = self.n_theta_max
            # Inner core is not stored here
            self.data = np.zeros((self.n_fields, self.nvar, self.n_theta_plot,
                                  self.n_r_max), precision)
            self.data_ic = np.zeros((self.n_fields, self.nvar,
                                     self.n_theta_plot, self.n_r_ic_max+2),
                                     precision)

        self.time = np.zeros(self.nvar, precision)

        # Read the data

        # If one skip the beginning, nevertheless read but do not store
        for i in range(self.var2-self.nvar):
            n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                movieDipLon, movieDipStrength, movieDipStrengthGeo = \
                infile.fort_read(precision)
            for ll in range(self.n_fields):
                dat = infile.fort_read(precision)
        # then read the remaining requested nvar lines
        for k in range(self.nvar):
            n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                movieDipLon, movieDipStrength, movieDipStrengthGeo = \
                infile.fort_read(precision)
            self.time[k] = t_movieS
            for ll in range(self.n_fields):
                dat = infile.fort_read(precision)
                if n_surface == 0:
                    dat = dat.reshape(shape)
                    if self.movtype in [1, 2, 3]:
                        datic = dat[self.n_r_max:, ...].T
                        dat = dat[:self.n_r_max, ...].T
                        self.data[ll, k, ...] = dat
                        self.data_ic[ll, k, ...] = datic
                    else:
                        self.data[ll, k, ...] = dat.T
                elif n_surface == 2:
                    dat = dat.reshape(shape)
                    if self.movtype in [1, 2, 3, 14]:
                        datic = dat[self.n_r_max:, :].T
                        dat = dat[:self.n_r_max, :].T
                        self.data[ll, k, ...] = dat
                        self.data_ic[ll, k, ...] = datic
                    else:
                        self.data[ll, k, ...] = dat.T
                elif n_surface == 3:
                    if self.movtype in [1, 2, 3, 14]:
                        len1 = (self.n_r_max*self.n_theta_max*2)
                        datoc = dat[:len1]
                        datic = dat[len1:]
                        datoc0 = datoc[:len(datoc)//2].reshape(self.n_r_max,
                                                               self.n_theta_max)
                        datoc1 = datoc[len(datoc)//2:].reshape(self.n_r_max,
                                                               self.n_theta_max)
                        dat = np.hstack((datoc0, datoc1))
                        datic0 = datic[:len(datic)//2].reshape(self.n_r_ic_max+2,
                                                               self.n_theta_max)
                        datic1 = datic[len(datic)//2:].reshape(self.n_r_ic_max+2,
                                                               self.n_theta_max)
                        dat = np.hstack((datoc0, datoc1))
                        datic = np.hstack((datic0, datic1))
                        self.data[ll, k, ...] = dat.T
                        self.data_ic[ll, k, ...] = datic.T
                    elif self.movtype in [8, 9]:
                        dat = dat.reshape(shape)
                        datic = dat[self.n_r_max:, :].T
                        dat = dat[:self.n_r_max, :].T
                        self.data_ic[ll, k, ...] = datic
                        self.data[ll, k, ...] = dat
                    elif self.movtype in [4, 5, 6, 7, 16, 17, 18, 47, 54, 91, 109]:
                        dat0 = dat[:len(dat)//2].reshape(shape)
                        dat1 = dat[len(dat)//2:].reshape(shape)
                        dat = np.hstack((dat0, dat1))
                        self.data[ll, k, ...] = dat.T
                    elif self.movtype in [10, 11, 12, 19, 92, 94, 95, 110]:
                        dat = dat.reshape(shape)
                        self.data[ll, k, ...] = dat.T
                else:
                    dat = dat.reshape(shape)
                    self.data[ll, k, ...] = dat.T
                if fluct:
                    self.data[ll, k, ...] -= self.data[ll, k, ...].mean(axis=0)

        infile.close()

        if normRad:
            for ll in range(self.n_fields):
                norm = np.sqrt(np.mean(self.data[ll, ...]**2, axis=1))
                norm = norm.mean(axis=0)
                self.data[ll, :, :, norm != 0.] /= norm[norm != 0.]

        if iplot:
            cmap = plt.get_cmap(cm)
            self.plot(ifield, cut, levels, cmap, png, step, normed, dpi,
                      bgcolor, deminc)
        if avg or std:
            cmap = plt.get_cmap(cm)
            self.avgStd(ifield, std, cut, levels, cmap)

    def __add__(self, new):
        """
        Built-in function to sum two movies

        .. note:: So far this function only works for two movies with the same
                  grid sizes. At some point, we might introduce grid
                  extrapolation to allow any summation.
        """
        out = copy.deepcopy(new)
        if new.time[0] == self.time[-1]:
            out.time = np.concatenate((self.time, new.time[1:]), axis=0)
            out.data = np.concatenate((self.data, new.data[:, 1:, ...]),
                                      axis=1)
            out.nvar = self.nvar+new.nvar-1
            out.var2 = out.nvar
        else:
            out.time = np.concatenate((self.time, new.time), axis=0)
            out.data = np.concatenate((self.data, new.data), axis=1)
            out.nvar = self.nvar+new.nvar
            out.var2 = out.nvar
        return out

    def avgStd(self, ifield=0, std=False, cut=0.5, levels=12, cmap='RdYlBu_r',
               ic=False):
        """
        Plot time-average or standard deviation

        :param ifield: in case of a multiple-field movie file, you can change
                       the default field displayed using the parameter ifield
        :type ifield: int
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
            avg = self.data[ifield, ...].std(axis=0)
            if ic:
                avg_ic = self.data_ic[ifield, ...].std(axis=0)
        else:
            avg = self.data[ifield, ...].mean(axis=0)
            if ic:
                avg_ic = self.data_ic[ifield, ...].mean(axis=0)
        vmin = - max(abs(avg.max()), abs(avg.min()))
        vmin = cut * vmin
        vmax = -vmin
        cs = np.linspace(vmin, vmax, levels)

        if self.surftype == 'phi_constant':
            if self.n_theta_plot == self.n_theta_max:
                th = np.linspace(np.pi/2., -np.pi/2., self.n_theta_max)
                fig = plt.figure(figsize=(4, 8))
                th0 = th
            else:
                th0 = np.linspace(np.pi/2., -np.pi/2., self.n_theta_max)
                th1 = np.linspace(np.pi/2., 3.*np.pi/2., self.n_theta_max)
                th = np.concatenate((th0, th1))
                fig = plt.figure(figsize=(6.5, 6))
                # Plotting trick using th0
                th0 = np.linspace(np.pi/2, np.pi/2+2.*np.pi, 2*self.n_theta_max)
            rr, tth = np.meshgrid(self.radius, th)
            xx = rr * np.cos(tth)
            yy = rr * np.sin(tth)
            xxout = rr.max() * np.cos(th0)
            yyout = rr.max() * np.sin(th0)
            xxin = rr.min() * np.cos(th0)
            yyin = rr.min() * np.sin(th0)
            if ic:
                rr, tth = np.meshgrid(self.radius_ic, th)
                xx_ic = rr * np.cos(tth)
                yy_ic = rr * np.sin(tth)
        elif self.surftype == 'r_constant':
            th = np.linspace(np.pi/2., -np.pi/2., self.n_theta_max)
            phi = np.linspace(-np.pi, np.pi, self.n_phi_tot)
            ttheta, pphi = np.meshgrid(th, phi)
            xx, yy = hammer2cart(ttheta, pphi)
            xxout, yyout = hammer2cart(th, -np.pi)
            xxin, yyin = hammer2cart(th, np.pi)
            fig = plt.figure(figsize=(8, 4))
        elif self.surftype == 'theta_constant':
            phi = np.linspace(0., 2.*np.pi, self.n_phi_tot)
            rr, pphi = np.meshgrid(self.radius, phi)
            xx = rr * np.cos(pphi)
            yy = rr * np.sin(pphi)
            xxout = rr.max() * np.cos(pphi)
            yyout = rr.max() * np.sin(pphi)
            xxin = rr.min() * np.cos(pphi)
            yyin = rr.min() * np.sin(pphi)
            fig = plt.figure(figsize=(6, 6))
        elif self.surftype == '3d volume':
            self.data = self.data[ifield, ..., 0]
            th = np.linspace(np.pi/2., -np.pi/2., self.n_theta_max)
            phi = np.linspace(-np.pi, np.pi, self.n_phi_tot)
            ttheta, pphi = np.meshgrid(th, phi)
            xx, yy = hammer2cart(ttheta, pphi)
            xxout, yyout = hammer2cart(th, -np.pi)
            xxin, yyin = hammer2cart(th, np.pi)
            fig = plt.figure(figsize=(8, 4))

        fig.subplots_adjust(top=0.99, right=0.99, bottom=0.01, left=0.01)
        ax = fig.add_subplot(111, frameon=False)
        ax.contourf(xx, yy, avg, cs, cmap=cmap, extend='both')
        if ic:
            ax.contourf(xx_ic, yy_ic, avg_ic, cs, cmap=cmap, extend='both')
        ax.plot(xxout, yyout, 'k-', lw=1.5)
        ax.plot(xxin, yyin, 'k-', lw=1.5)
        ax.axis('off')

    def plot(self, ifield=0, cut=0.5, levels=12, cmap='RdYlBu_r', png=False,
             step=1, normed=False, dpi=80, bgcolor=None, deminc=True,
             ic=False):
        """
        Plotting function (it can also write the png files)

        :param ifield: in case of a multiple-field movie file, you can change
                       the default field displayed using the parameter ifield
        :type ifield: int
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
            plt.ioff()
            if not os.path.exists('movie'):
                os.mkdir('movie')
        else:
            plt.ion()

        if not normed:
            vmin = - max(abs(self.data[ifield, ...].max()),
                         abs(self.data[ifield, ...].min()))
            vmin = cut * vmin
            vmax = -vmin
            # vmin, vmax = self.data.min(), self.data.max()
            cs = np.linspace(vmin, vmax, levels)

        if self.surftype == 'phi_constant':
            if self.n_theta_plot == self.n_theta_max:
                th = np.linspace(np.pi/2., -np.pi/2., self.n_theta_max)
                fig = plt.figure(figsize=(4, 8))
                th0 = th
            else:
                th0 = np.linspace(np.pi/2., -np.pi/2., self.n_theta_max)
                th1 = np.linspace(np.pi/2., 3.*np.pi/2., self.n_theta_max)
                th = np.concatenate((th0, th1))
                fig = plt.figure(figsize=(6.5, 6))
                # Plotting trick using th0
                th0 = np.linspace(np.pi/2, np.pi/2+2.*np.pi, 2*self.n_theta_max)
            rr, tth = np.meshgrid(self.radius, th)
            xx = rr * np.cos(tth)
            yy = rr * np.sin(tth)
            xxout = rr.max() * np.cos(th0)
            yyout = rr.max() * np.sin(th0)
            xxin = rr.min() * np.cos(th0)
            yyin = rr.min() * np.sin(th0)

            if ic:
                rr, tth = np.meshgrid(self.radius_ic, th)
                xx_ic = rr * np.cos(tth)
                yy_ic = rr * np.sin(tth)
        elif self.surftype == 'r_constant':
            th = np.linspace(np.pi/2., -np.pi/2., self.n_theta_max)
            if deminc:
                phi = np.linspace(-np.pi, np.pi, self.n_phi_tot*self.minc+1)
                xxout, yyout = hammer2cart(th, -np.pi)
                xxin, yyin = hammer2cart(th, np.pi)
            else:
                phi = np.linspace(-np.pi/self.minc, np.pi/self.minc,
                                  self.n_phi_tot)
                xxout, yyout = hammer2cart(th, -np.pi/self.minc)
                xxin, yyin = hammer2cart(th, np.pi/self.minc)
            ttheta, pphi = np.meshgrid(th, phi)
            xx, yy = hammer2cart(ttheta, pphi)
            fig = plt.figure(figsize=(8, 4))
        elif self.surftype == 'theta_constant':
            if deminc:
                phi = np.linspace(0., 2.*np.pi, self.n_phi_tot*self.minc+1)
            else:
                phi = np.linspace(0., 2.*np.pi/self.minc, self.n_phi_tot)
            rr, pphi = np.meshgrid(self.radius, phi)
            xx = rr * np.cos(pphi)
            yy = rr * np.sin(pphi)
            xxout = rr.max() * np.cos(pphi)
            yyout = rr.max() * np.sin(pphi)
            xxin = rr.min() * np.cos(pphi)
            yyin = rr.min() * np.sin(pphi)
            if ic:
                rr, pphi = np.meshgrid(self.radius_ic, phi)
                xx_ic = rr * np.cos(pphi)
                yy_ic = rr * np.sin(pphi)
            fig = plt.figure(figsize=(6, 6))
        elif self.surftype == '3d volume':
            self.data = self.data[ifield, ..., 0]
            th = np.linspace(np.pi/2., -np.pi/2., self.n_theta_max)
            phi = np.linspace(-np.pi, np.pi, self.n_phi_tot)
            ttheta, pphi = np.meshgrid(th, phi)
            xx, yy = hammer2cart(ttheta, pphi)
            xxout, yyout = hammer2cart(th, -np.pi)
            xxin, yyin = hammer2cart(th, np.pi)
            fig = plt.figure(figsize=(8, 4))

        fig.subplots_adjust(top=0.99, right=0.99, bottom=0.01, left=0.01)
        ax = fig.add_subplot(111, frameon=False)

        for k in range(self.nvar):
            if k == 0:
                if normed:
                    vmin = - max(abs(self.data[ifield, k, ...].max()),
                                 abs(self.data[ifield, k, ...].min()))
                    vmin = cut * vmin
                    vmax = -vmin
                    cs = np.linspace(vmin, vmax, levels)
                if self.surftype in ['r_constant', 'theta_constant']:
                    if deminc:
                        dat = symmetrize(self.data[ifield, k, ...], self.minc)
                        if ic:
                            datic = symmetrize(self.data_ic[ifield, k, ...],
                                               self.minc)
                    else:
                        dat = self.data[ifield, k, ...]
                        if ic:
                            datic = self.data_ic[ifield, k, ...]
                    im = ax.contourf(xx, yy, dat, cs, cmap=cmap, extend='both')
                    if ic:
                        ax.contourf(xx_ic, yy_ic, datic, cs, cmap=cmap,
                                    extend='both')
                else:
                    im = ax.contourf(xx, yy, self.data[ifield, k, ...], cs,
                                     cmap=cmap, extend='both')
                    if ic:
                        im_ic = ax.contourf(xx_ic, yy_ic,
                                            self.data_ic[ifield, k, ...], cs,
                                            cmap=cmap, extend='both')
                ax.plot(xxout, yyout, 'k-', lw=1.5)
                ax.plot(xxin, yyin, 'k-', lw=1.5)
                ax.axis('off')
                man = plt.get_current_fig_manager()
                man.canvas.draw()
            if k != 0 and k % step == 0:
                if not png:
                    print(k+self.var2-self.nvar)
                plt.cla()
                if normed:
                    vmin = - max(abs(self.data[ifield, k, ...].max()),
                                 abs(self.data[ifield, k, ...].min()))
                    vmin = cut * vmin
                    vmax = -vmin
                    cs = np.linspace(vmin, vmax, levels)
                if self.surftype in ['r_constant', 'theta_constant']:
                    if deminc:
                        dat = symmetrize(self.data[ifield, k, ...], self.minc)
                        if ic:
                            datic = symmetrize(self.data_ic[ifield, k, ...],
                                               self.minc)
                    else:
                        dat = self.data[ifield, k, ...]
                        if ic:
                            datic = self.data_ic[ifield, k, ...]
                    ax.contourf(xx, yy, dat, cs, cmap=cmap, extend='both')
                    if ic:
                        ax.contourf(xx_ic, yy_ic, datic, cs, cmap=cmap,
                                    extend='both')
                else:
                    ax.contourf(xx, yy, self.data[ifield, k, ...], cs,
                                cmap=cmap, extend='both')
                    if ic:
                        ax.contourf(xx_ic, yy_ic, self.data_ic[ifield, k, ...],
                                    cs, cmap=cmap, extend='both')
                ax.plot(xxout, yyout, 'k-', lw=1.5)
                ax.plot(xxin, yyin, 'k-', lw=1.5)
                ax.axis('off')
                man.canvas.draw()
            if png:
                filename = 'movie/img{:05d}.png'.format(k)
                print('write {}'.format(filename))
                # st = 'echo {}'.format(ivar) + ' > movie/imgmax'
                if bgcolor is not None:
                    fig.savefig(filename, facecolor=bgcolor, dpi=dpi)
                else:
                    fig.savefig(filename, dpi=dpi)

    def timeLongitude(self, ifield=0, removeMean=True, lat0=0., levels=12,
                      cm='RdYlBu_r', deminc=True):
        """
        Plot the time-longitude diagram (input latitude can be chosen)

        :param ifield: in case of a multiple-field movie file, you can change
                       the default field displayed using the parameter ifield
        :type ifield: int
        :param lat0: value of the latitude
        :type lat0: float
        :param levels: number of contour levels
        :type levels: int
        :param cm: name of the colormap
        :type cm: str
        :param deminc: a logical to indicate if one wants do get rid of the
                       possible azimuthal symmetry
        :type deminc: bool
        :param removeMean: remove the time-averaged part when set to True
        :type removeMean: bool
        """

        if removeMean:
            datCut = self.data[ifield, ...]-self.data[ifield, ...].mean(axis=0)
        else:
            datCut = self.data[ifield, ...]


        th = np.linspace(np.pi/2., -np.pi/2., self.n_theta_max)
        lat0 *= np.pi/180.
        mask = np.where(abs(th-lat0) == abs(th-lat0).min(), 1, 0)
        idx = np.nonzero(mask)[0][0]

        th = np.linspace(np.pi/2., -np.pi/2., self.n_theta_max)
        if deminc:
            phi = np.linspace(-np.pi, np.pi, self.minc*self.n_phi_tot+1)
        else:
            phi = np.linspace(-np.pi/self.minc, np.pi/self.minc,
                              self.n_phi_tot)

        if deminc:
            dat = np.zeros((self.nvar, self.minc*self.n_phi_tot+1), np.float64)
        else:
            dat = np.zeros((self.nvar, self.n_phi_tot), np.float64)

        for k in range(self.nvar):
            if deminc:
                dat[k, :] = symmetrize(datCut[k, :, idx], self.minc)
            else:
                dat[k, :] = datCut[k, :, idx]


        fig = plt.figure()
        ax = fig.add_subplot(111)
        vmin = -max(abs(dat.max()), abs(dat.min()))
        vmax = -vmin
        cs = np.linspace(vmin, vmax, levels)
        ax.contourf(phi, self.time, dat, cs, cmap=plt.get_cmap(cm))

        ax.set_xlabel('Longitude')
        ax.set_ylabel('Time')

        m_max = self.n_phi_tot/3
        w2 = np.fft.fft2(dat)
        w2 = abs(w2[1:self.nvar/2+1, 0:m_max+1])

        dw = 2.*np.pi/(self.time[-1]-self.time[0])
        omega = dw*np.arange(self.nvar)
        omega = omega[1:self.nvar/2+1]
        ms = np.arange(m_max+1)

        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax1.contourf(ms, omega, np.log10(w2), 65, cmap=plt.get_cmap('jet'))
        ax1.set_yscale('log')
        ax1.set_xlabel(r'Azimuthal wavenumber')
        ax1.set_ylabel(r'Frequency')


if __name__ == '__main__':
    Movie(step=1)
    plt.show()
