# -*- coding: utf-8 -*-
import glob
import re
import os
import numpy as N
import pylab as P
from npfile import *

__author__  = "$Author$"
__date__   = "$Date$"
__version__ = "$Revision$"


def hammer2cart(ttheta, pphi):
    xx = 2.*N.sqrt(2.) * N.cos(ttheta)*N.sin(pphi/2.)\
         /N.sqrt(1.+N.cos(ttheta)*N.cos(pphi/2.))
    yy = N.sqrt(2.) * N.sin(ttheta)\
         /N.sqrt(1.+N.cos(ttheta)*N.cos(pphi/2.))
    return xx, yy



class Movie:

    def __init__(self, file=None, avg=False, iplot=False, step=1, png=False,
                 lastvar=None, nvar='all', levels=12, cmap='RdYlBu_r', cut=0.5,
                 bgcolor=None, fluct=False, centered=True, normed=False):
        """
        :param nvar: the number of lines of the movie file we want to plot
                     starting from the last line
        :param lastvar: the rank of the last line to be read
        :param step: the stepping between two lines             
        :param levels: the number of contour levels
        :param cmap: the name of the color map
        :param png: if png=True, write the png outputs
        :param fluct: if fluct=True, substract the axisymmetric part
        :param centered: if centered=True, the colormap is centered around 0
        :param normed: if normed=True, the colormap is rescaled every timestep,
                       otherwise it is computed from the first line
        """
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
        print filename
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
        n_type, n_surface, const, n_fields = infile.fort_read('f')
        movtype = infile.fort_read('f')
        self.movtype = int(movtype)
        n_surface = int(n_surface)

        # RUN PARAMETERS
        runid = infile.fort_read('|S64')
        n_r_mov_tot, n_r_max, n_theta_max, n_phi_tot, self.minc, self.ra, \
             self.ek, self.pr, self.prmag, self.radratio, self.tScale = infile.fort_read('f')
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
            surftype = '3d volume'
            shape = (n_r_mov_tot+2, self.n_theta_max, self.n_phi_tot)
        elif n_surface == 1:
            surftype = 'r_constant'
            shape = (self.n_theta_max, self.n_phi_tot)
        elif n_surface == 2:
            surftype = 'theta_constant'
            shape = (self.n_r_max, self.n_phi_tot)
        elif n_surface == 3:
            surftype = 'phi_constant'
            if self.movtype in [8, 9]:
                shape = (n_r_mov_tot+2, self.n_theta_max)
            else:
                shape = (self.n_r_max, self.n_theta_max)

        cmap = P.get_cmap(cmap)
        self. time = N.zeros(self.nvar, 'f')
        if not png:
            P.ion()
            for i in range(self.var2-self.nvar):
                n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                     movieDipLon, movieDipStrength, \
                     movieDipStrengthGeo = infile.fort_read('f')
                self.data = infile.fort_read('f', shape=shape)
            for k in range(self.nvar):
                n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                     movieDipLon, movieDipStrength, \
                     movieDipStrengthGeo = infile.fort_read('f')
                self.time[k] = t_movieS
                if k == 0:
                    self.data = infile.fort_read('f', shape=shape).T
                    if fluct:
                        self.data = self.data-self.data.mean(axis=0)
                    if self.movtype in [8, 9]:
                        self.data = self.data[:, :self.n_r_max] # remove inner core
                    if surftype == 'phi_constant':
                        th = N.linspace(N.pi/2., -N.pi/2., self.n_theta_max)
                        rr, tth = N.meshgrid(self.radius, th)
                        xx = rr * N.cos(tth)
                        yy = rr * N.sin(tth)
                        xxout = rr.max() * N.cos(th)
                        yyout = rr.max() * N.sin(th)
                        xxin = rr.min() * N.cos(th)
                        yyin = rr.min() * N.sin(th)
                        fig = P.figure(figsize=(4, 8))
                        P.subplots_adjust(top=0.99, right=0.99, bottom=0.01,
                                          left=0.01)
                    elif surftype == 'r_constant':
                        th = N.linspace(N.pi/2., -N.pi/2., self.n_theta_max)
                        phi = N.linspace(-N.pi, N.pi, self.n_phi_tot)
                        ttheta, pphi = N.meshgrid(th, phi)
                        xx, yy = hammer2cart(ttheta, pphi)
                        #xxout, yyout = hammer2cart(N.pi/4., phi)
                        #xxin, yyin = hammer2cart(-N.pi/4., phi)
                        xxout, yyout = hammer2cart(th, -N.pi)
                        xxin, yyin = hammer2cart(th, N.pi)
                        fig = P.figure(figsize=(8, 4))
                        P.subplots_adjust(top=0.99, right=0.99, bottom=0.01,
                                          left=0.01)
                    elif surftype == 'theta_constant':
                        phi = N.linspace(0., 2.*N.pi, self.n_phi_tot)
                        rr, pphi = N.meshgrid(self.radius, phi)
                        xx = rr * N.cos(pphi)
                        yy = rr * N.sin(pphi)
                        xxout = rr.max() * N.cos(pphi)
                        yyout = rr.max() * N.sin(pphi)
                        xxin = rr.min() * N.cos(pphi)
                        yyin = rr.min() * N.sin(pphi)
                        fig = P.figure(figsize=(6, 6))
                        P.subplots_adjust(top=0.99, bottom=0.01, right=0.99,
                                          left=0.01)
                    elif surftype == '3d volume':
                        self.data = self.data[..., 0]
                        th = N.linspace(N.pi/2., -N.pi/2., self.n_theta_max)
                        phi = N.linspace(-N.pi, N.pi, self.n_phi_tot)
                        ttheta, pphi = N.meshgrid(th, phi)
                        xx, yy = hammer2cart(ttheta, pphi)
                        xxout, yyout = hammer2cart(th, -N.pi)
                        xxin, yyin = hammer2cart(th, N.pi)
                        fig = P.figure(figsize=(8, 4))
                        P.subplots_adjust(top=0.99, right=0.99, bottom=0.01,
                                          left=0.01)
                    ax = fig.add_subplot(111, frameon=False)
                    P.axis('off')
                    if centered:
                        vmin = - max(abs(self.data.max()),
                                     abs(self.data.min()))
                        vmin = cut * vmin
                        vmax = -vmin
                        cs = N.linspace(vmin, vmax, levels)
                    else:
                        vmin = self.data.min()
                        vmax = self.data.max()
                        cs = N.linspace(vmin, vmax, levels)
                    im = ax.contourf(xx, yy, self.data, cs, cmap=cmap,
                                     extend='both')
                    ax.plot(xxout, yyout, 'k-', lw=1.5)
                    ax.plot(xxin, yyin, 'k-', lw=1.5)
                    #P.colorbar(im)
                    man = P.get_current_fig_manager()
                    man.canvas.draw()
                elif k != 0 and k % step == 0:
                    print k+self.var2-self.nvar
                    P.cla()
                    self.data = infile.fort_read('f', shape=shape).T
                    if self.movtype in [8, 9]:
                        self.data = self.data[:, :self.n_r_max] # remove inner core
                    if surftype == '3d volume':
                        self.data = self.data[..., 0]
                    if fluct:
                        self.data = self.data-self.data.mean(axis=0)
                    if normed:
                        vmin = - max(abs(self.data.max()),
                                     abs(self.data.min()))
                        vmin = cut * vmin
                        vmax = -vmin
                        cs = N.linspace(vmin, vmax, levels)
                    im = ax.contourf(xx, yy, self.data, cs, cmap=cmap,
                                     extend='both')

                    ax.plot(xxout, yyout, 'k-', lw=1.5)
                    ax.plot(xxin, yyin, 'k-', lw=1.5)
                    #P.colorbar(im)
                    P.axis('off')
                    man.canvas.draw()
                else: # On lit quand meme
                    self.data = infile.fort_read('f', shape=shape)
        else:
            P.ioff()
            if not os.path.exists('movie'):
                os.mkdir('movie')
            for k in range(0, self.nvar, step):
                n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                     movieDipLon, movieDipStrength, \
                     movieDipStrengthGeo = infile.fort_read('f')
                self.time[k] = t_movieS
                if k == 0:
                    self.dat = infile.fort_read('f', shape=shape).T
                    if self.movtype in [8, 9]:
                        self.dat = self.dat[:, :self.n_r_max] # remove inner core
                    if surftype == 'phi_constant':
                        th = N.linspace(N.pi/2., -N.pi/2., self.n_theta_max)
                        rr, tth = N.meshgrid(self.radius, th)
                        xx = rr * N.cos(tth)
                        yy = rr * N.sin(tth)
                        xxout = rr.max() * N.cos(th)
                        yyout = rr.max() * N.sin(th)
                        xxin = rr.min() * N.cos(th)
                        yyin = rr.min() * N.sin(th)
                        fig = P.figure(figsize=(4, 8))
                        P.subplots_adjust(top=0.99, right=0.99, bottom=0.01,
                                          left=0.01)
                    elif surftype == 'r_constant':
                        th = N.linspace(N.pi/2., -N.pi/2., self.n_theta_max)
                        phi = N.linspace(-N.pi, N.pi, self.n_phi_tot)
                        ttheta, pphi = N.meshgrid(th, phi)
                        xx, yy = hammer2cart(ttheta, pphi)
                        xxout, yyout = hammer2cart(th, -N.pi)
                        xxin, yyin = hammer2cart(th, N.pi)

                        fig = P.figure(figsize=(8, 4))
                        P.subplots_adjust(top=0.99, right=0.99, bottom=0.01,
                                          left=0.01)
                    ax = fig.add_subplot(111, frameon=False)
                    P.axis('off')
                    if centered:
                        vmin = - max(abs(self.dat.max()),
                                     abs(self.dat.min()))
                        vmin = cut * vmin
                        vmax = -vmin
                        cs = N.linspace(vmin, vmax, levels)
                    else:
                        vmin = self.dat.min()
                        vmax = self.dat.max()
                        cs = N.linspace(vmin, vmax, levels)
                    im = ax.contourf(xx, yy, self.dat, cs, cmap=cmap,
                                     extend='both')
                    man = P.get_current_fig_manager()
                    man.canvas.draw()
                elif k != 0 and k % step == 0:
                    P.cla()
                    self.data = infile.fort_read('f', shape=shape).T
                    if self.movtype in [8, 9]:
                        self.data = self.data[:, :self.n_r_max] # remove inner core
                    im = ax.contourf(xx, yy, self.data, cs, cmap=cmap,
                                     extend='both')
                    ax.plot(xxout, yyout, 'k-', lw=1.5)
                    ax.plot(xxin, yyin, 'k-', lw=1.5)
                    P.axis('off')
                    man.canvas.draw()
                else: # On lit quand meme
                    self.data = infile.fort_read('f', shape=shape)
                filename = 'movie/img%05d.png' % k
                print 'write %s' % filename
                #st = 'echo %i' % ivar + ' > movie/imgmax'
                if bgcolor is not None:
                    P.savefig(filename, facecolor=bgcolor)
                else:
                    P.savefig(filename)

        #if len(self.dat.shape) > 2:
            #self.dat = self.data

        # Choose which type of representation
        #if iplot:
            #if surftype == 'r_constant':
                #self.surf()
            #elif surftype == 'theta_constant':
                #self.equat()
            #elif surftype == 'phi_constant':
                #self.avg()



if __name__ == '__main__':
    Movie(step=1)
    P.show()
