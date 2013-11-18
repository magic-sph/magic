# -*- coding: utf-8 -*-
import glob
import re
import os
import copy
import numpy as N
import pylab as P
from npfile import *
from magic.libmagic import symmetrize, hammer2cart

__author__  = "$Author$"
__date__   = "$Date$"
__version__ = "$Revision$"



class Movie:

    def __init__(self, file=None, iplot=True, step=1, png=False,
                 lastvar=None, nvar='all', levels=12, cmap='RdYlBu_r', cut=0.5,
                 bgcolor=None, fluct=False, normed=False, avg=False, 
                 std=False, dpi=80):
        """
        :param nvar: the number of lines of the movie file we want to plot
                     starting from the last line
        :param png: if png=True, write the png files instead of display
        :param iplot: if iplot=True, display otherwise just read
        :param lastvar: the rank of the last line to be read
        :param step: the stepping between two lines             
        :param levels: the number of contour levels
        :param cmap: the name of the color map
        :param png: if png=True, write the png outputs
        :param fluct: if fluct=True, substract the axisymmetric part
        :param normed: if normed=True, the colormap is rescaled every timestep,
                       otherwise it is computed from the first line
        :param avg: if avg=True, time-average is displayed
        :param avg: if std=True, standard deviation is displayed
        :param dpi: dot per inch when saving PNGs
        """

        if avg or std:
            iplot = False
        if file == None:
            dat = glob.glob('*_mov.*')
            str1 = 'Which movie do you want ?\n'
            for k, movie in enumerate(dat):
                str1 += ' %i) %s\n' % (k+1, movie)
            index = input(str1)
            try:
                filename = dat[index-1]
            except IndexError:
                print 'Non valid index: %s has been chosen instead' % dat[0] 
                filename = dat[0]

        else:
            filename = file
        mot = re.compile(r'.*_mov\.(.*)')
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
        n_type, n_surface, const, n_fields = infile.fort_read('f')
        movtype = infile.fort_read('f')
        n_fields = int(n_fields)
        if n_fields > 1:
            print '!!! Warning: several fields in the movie file !!!'
            print '!!! Only the last one will be displayed       !!!'
            print '!!! For TO_mov.TAG, use TOMovie(...) instead  !!!'
        self.movtype = int(movtype[0])
        n_surface = int(n_surface)
        print 'movtype, nsurf', self.movtype, n_surface

        # RUN PARAMETERS
        runid = infile.fort_read('|S64')
        n_r_mov_tot, n_r_max, n_theta_max, n_phi_tot, self.minc, self.ra, \
             self.ek, self.pr, self.prmag, self.radratio, self.tScale = infile.fort_read('f')
        self.minc = int(self.minc)
        n_r_mov_tot = int(n_r_mov_tot)
        self.n_r_max = int(n_r_max)
        self.n_theta_max = int(n_theta_max)
        self.n_phi_tot = int(n_phi_tot)

        # GRID
        self.radius = infile.fort_read('f')
        self.radius = self.radius[:self.n_r_max] # remove inner core
        self.theta = infile.fort_read('f')
        self.phi = infile.fort_read('f')

        if n_surface == 0:
            self.surftype = '3d volume'
            shape = (n_r_mov_tot+2, self.n_theta_max, self.n_phi_tot)
        elif n_surface == 1:
            self.surftype = 'r_constant'
            shape = (self.n_theta_max, self.n_phi_tot)
            self.data = N.zeros((self.nvar, self.n_phi_tot, self.n_theta_max), 'f')
        elif n_surface == 2:
            self.surftype = 'theta_constant'
            if self.movtype in [1, 2, 3]: # read inner core
                shape = (n_r_mov_tot+2, self.n_phi_tot)
            else:
                shape = (self.n_r_max, self.n_phi_tot)
            self.data = N.zeros((self.nvar, self.n_phi_tot, self.n_r_max), 'f')
        elif n_surface == 3:
            self.surftype = 'phi_constant'
            if self.movtype in [1, 2, 3]: # read inner core
                shape = (n_r_mov_tot+2, 2*self.n_theta_max)
            elif self.movtype in [8, 9]:
                shape = (n_r_mov_tot+2, self.n_theta_max)
            elif self.movtype in [4, 5, 6, 7, 16, 17, 18, 47, 54, 91]:
                shape = (self.n_r_max, 2*self.n_theta_max)
            elif self.movtype in [10, 11, 12, 19]:
                shape = (self.n_r_max, self.n_theta_max)
            # Inner core is not stored here
            self.data = N.zeros((self.nvar, self.n_theta_max, self.n_r_max), 'f')

        self.time = N.zeros(self.nvar, 'f')

        # READ the data

        # If one skip the beginning, nevertheless read but do not store
        for i in range(self.var2-self.nvar):
            n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                                   movieDipLon, movieDipStrength, \
                            movieDipStrengthGeo = infile.fort_read('f')
            for ll in range(n_fields):
                dat = infile.fort_read('f', shape=shape)
        # then read the remaining requested nvar lines
        for k in range(self.nvar):
            n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                                   movieDipLon, movieDipStrength, \
                            movieDipStrengthGeo = infile.fort_read('f')
            self.time[k] = t_movieS
            for ll in range(n_fields):
                dat = infile.fort_read('f', shape=shape)
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
                    elif self.movtype in [10, 11, 12, 19]:
                        self.data[k, ...] = dat.T
                else:
                    self.data[k, ...] = dat.T
            if fluct:
                self.data[k, ...] = self.data[k, ...]-self.data[k, ...].mean(axis=0)

        if iplot:
            cmap = P.get_cmap(cmap)
            self.plot(cut, levels, cmap, png, step, normed, dpi, bgcolor)
        if avg or std:
            cmap = P.get_cmap(cmap)
            self.avgStd(std, cut, levels, cmap)

    def __add__(self, new):
        """
        Built-in function to sum two movies
        So far only works for same grid sizes: at some point we should introduce
        extrapolation to allow any summation
        """
        out = copy.deepcopy(new)
        out.time = N.concatenate((self.time, new.time), axis=0)
        out.data = N.concatenate((self.data, new.data), axis=0)
        out.nvar = out.nvar+new.nvar
        out.var2 = out.nvar
        return out

    def avgStd(self, std=False, cut=0.5, levels=12, cmap='RdYlBu_r'):
        """
        plot time-average or standard deviation

        :param std: if std=True standard deviation is computed instead avg
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
             normed=False, dpi=80, bgcolor=None):
        """
        plotting subroutine (can also write the png files)
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
            phi = N.linspace(-N.pi, N.pi, self.n_phi_tot*self.minc+1)
            ttheta, pphi = N.meshgrid(th, phi)
            xx, yy = hammer2cart(ttheta, pphi)
            xxout, yyout = hammer2cart(th, -N.pi)
            xxin, yyin = hammer2cart(th, N.pi)
            fig = P.figure(figsize=(8, 4))
        elif self.surftype == 'theta_constant':
            phi = N.linspace(0., 2.*N.pi, self.n_phi_tot*self.minc+1)
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
        ax = fig.add_subplot(111)

        for k in range(self.nvar):
            if k == 0:
                if normed:
                    vmin = - max(abs(self.data[k, ...].max()), abs(self.data[k, ...].min()))
                    vmin = cut * vmin
                    vmax = -vmin
                    cs = N.linspace(vmin, vmax, levels)
                if self.surftype in ['r_constant', 'theta_constant']:
                    im = ax.contourf(xx, yy, symmetrize(self.data[k, ...], self.minc),
                                     cs, cmap=cmap, extend='both')
                else:
                    im = ax.contourf(xx, yy, self.data[k, ...], cs, cmap=cmap, extend='both')
                ax.plot(xxout, yyout, 'k-', lw=1.5)
                ax.plot(xxin, yyin, 'k-', lw=1.5)
                man = P.get_current_fig_manager()
                man.canvas.draw()
            if k !=0 and k % step == 0:
                if not png:
                    print k+self.var2-self.nvar
                P.cla()
                if normed:
                    vmin = - max(abs(self.data[k, ...].max()), abs(self.data[k, ...].min()))
                    vmin = cut * vmin
                    vmax = -vmin
                    cs = N.linspace(vmin, vmax, levels)
                if self.surftype in ['r_constant', 'theta_constant']:
                    im = ax.contourf(xx, yy, symmetrize(self.data[k, ...], self.minc),
                                     cs, cmap=cmap, extend='both')
                else:
                    im = ax.contourf(xx, yy, self.data[k, ...], cs, cmap=cmap, extend='both')
                ax.plot(xxout, yyout, 'k-', lw=1.5)
                ax.plot(xxin, yyin, 'k-', lw=1.5)
                ax.axis('off')
                man.canvas.draw()
            if png:
                filename = 'movie/img%05d.png' % k
                print 'write %s' % filename
                #st = 'echo %i' % ivar + ' > movie/imgmax'
                if bgcolor is not None:
                    fig.savefig(filename, facecolor=bgcolor, dpi=dpi)
                else:
                    fig.savefig(filename, dpi=dpi)



if __name__ == '__main__':
    Movie(step=1)
    P.show()
