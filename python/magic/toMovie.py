# -*- coding: utf-8 -*-
import glob
import re
import os
import numpy as N
import pylab as P
from npfile import *

__author__  = "$Author: gastine $"
__date__   = "$Date: 2013-03-07 15:06:48 +0100 (jeu. 07 mars 2013) $"
__version__ = "$Revision: 456 $"


class TOMovie:

    def __init__(self, file=None, iplot=True, cm='RdYlBu_r',
                 cut=0.8, levels=16, avg=True):
        """
        :param file: the filename of the TO_mov file
        :param cmap: the name of the color map
        :param levels: the number of contour levels
        :param cut: a parameter to change the maxima of the contour levels
        :param iplot: a boolean to specify if one wants to plot or not the
                      results
        :param avg: time average of the different forces
        """
                 
        if file == None:
            dat = glob.glob('TO_mov.*')
            str1 = 'Which TO movie do you want ?\n'
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
        n_type, n_surface, const, n_fields = infile.fort_read('f')
        movtype = infile.fort_read('f')
        n_fields = int(n_fields)
        self.movtype = N.asarray(movtype)
        n_surface = int(n_surface)

        # RUN PARAMETERS
        runid = infile.fort_read('|S64')
        n_r_mov_tot, n_r_max, n_theta_max, n_phi_tot, self.minc, self.ra, \
             self.ek, self.pr, self.prmag, \
             self.radratio, self.tScale = infile.fort_read('f')
        n_r_mov_tot = int(n_r_mov_tot)
        self.n_r_max = int(n_r_max)
        self.n_theta_max = int(n_theta_max)
        self.n_phi_tot = int(n_phi_tot)

        # GRID
        self.radius = infile.fort_read('f')
        self.radius = self.radius[:self.n_r_max] # remove inner core
        self.theta = infile.fort_read('f')
        self.phi = infile.fort_read('f')

        surftype = 'phi_constant'
        shape = (self.n_r_max, self.n_theta_max)

        self.time = N.zeros(self.nvar, 'f')
        self.asVphi = N.zeros((self.nvar, self.n_theta_max,self.n_r_max), 'f')
        self.rey = N.zeros_like(self.asVphi)
        self.adv = N.zeros_like(self.asVphi)
        self.visc = N.zeros_like(self.asVphi)
        self.lorentz = N.zeros_like(self.asVphi)
        self.coriolis = N.zeros_like(self.asVphi)
        self.dtVp = N.zeros_like(self.asVphi)

        # READ the data
        for k in range(self.nvar):
            n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                 movieDipLon, movieDipStrength, \
                 movieDipStrengthGeo = infile.fort_read('f')
            self.time[k] = t_movieS
            self.asVphi[k, ...] = infile.fort_read('f', shape=shape).T
            self.rey[k, ...] = infile.fort_read('f', shape=shape).T
            self.adv[k, ...] = infile.fort_read('f', shape=shape).T
            self.visc[k, ...] = infile.fort_read('f', shape=shape).T
            self.lorentz[k, ...] = infile.fort_read('f', shape=shape).T
            self.coriolis[k, ...] = infile.fort_read('f', shape=shape).T
            self.dtVp[k, ...] = infile.fort_read('f', shape=shape).T

        if iplot:
            cmap = P.get_cmap(cm)
            self.plot(cm, cut, levels, avg)

    def plot(self, cm, cut, levs, avg):
        th = N.linspace(N.pi/2., -N.pi/2., self.n_theta_max)
        rr, tth = N.meshgrid(self.radius, th)
        xx = rr * N.cos(tth)
        yy = rr * N.sin(tth)
        xxout = rr.max() * N.cos(th)
        yyout = rr.max() * N.sin(th)
        xxin = rr.min() * N.cos(th)
        yyin = rr.min() * N.sin(th)
        fig = P.figure(figsize=(20, 5))
        P.subplots_adjust(top=0.99, right=0.99, bottom=0.01, left=0.01, wspace=0.01)

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
            cs = N.linspace(vmin, vmax, levs)

            ax = fig.add_subplot(181)
            P.axis('off')
            im = ax.contourf(xx, yy, asVp, levs, cmap=cm, extend='both')
            ax.plot(xxout, yyout, 'k-', lw=1.5)
            ax.plot(xxin, yyin, 'k-', lw=1.5)
            ax.text(0.05, 0., 'uphi', fontsize=20, horizontalalignment='left',
                                  verticalalignment='center')
            ax = fig.add_subplot(182)
            P.axis('off')
            im = ax.contourf(xx, yy, adv, cs, cmap=cm, extend='both')
            ax.plot(xxout, yyout, 'k-', lw=1.5)
            ax.plot(xxin, yyin, 'k-', lw=1.5)
            ax.text(0.05, 0., 'Adv', fontsize=20, horizontalalignment='left',
                                  verticalalignment='center')
            ax = fig.add_subplot(183)
            P.axis('off')
            im = ax.contourf(xx, yy, rey, cs, cmap=cm, extend='both')
            ax.plot(xxout, yyout, 'k-', lw=1.5)
            ax.plot(xxin, yyin, 'k-', lw=1.5)
            ax.text(0.05, 0., 'Rey', fontsize=20, horizontalalignment='left',
                                  verticalalignment='center')
            ax = fig.add_subplot(184)
            P.axis('off')
            im = ax.contourf(xx, yy, vis, cs, cmap=cm, extend='both')
            ax.plot(xxout, yyout, 'k-', lw=1.5)
            ax.plot(xxin, yyin, 'k-', lw=1.5)
            ax.text(0.05, 0., 'Visc', fontsize=20, horizontalalignment='left',
                                  verticalalignment='center')
            ax = fig.add_subplot(185)
            P.axis('off')
            im = ax.contourf(xx, yy, lor, cs, cmap=cm, extend='both')
            ax.plot(xxout, yyout, 'k-', lw=1.5)
            ax.plot(xxin, yyin, 'k-', lw=1.5)
            ax.text(0.05, 0., 'Lo.', fontsize=20, horizontalalignment='left',
                                  verticalalignment='center')
            ax = fig.add_subplot(186)
            P.axis('off')
            im = ax.contourf(xx, yy, cor, cs, cmap=cm, extend='both')
            ax.plot(xxout, yyout, 'k-', lw=1.5)
            ax.plot(xxin, yyin, 'k-', lw=1.5)
            ax.text(0.05, 0., 'Cor.', fontsize=20, horizontalalignment='left',
                                  verticalalignment='center')
            ax = fig.add_subplot(187)
            P.axis('off')
            balance = adv+cor+vis+lor+rey
            im = ax.contourf(xx, yy, balance, cs, cmap=cm, extend='both')
            ax.plot(xxout, yyout, 'k-', lw=1.5)
            ax.plot(xxin, yyin, 'k-', lw=1.5)
            ax.text(0.05, 0., 'sum', fontsize=20, horizontalalignment='left',
                                  verticalalignment='center')
            ax = fig.add_subplot(188)
            P.axis('off')
            im = ax.contourf(xx, yy, dt, cs, cmap=cm, extend='both')
            ax.plot(xxout, yyout, 'k-', lw=1.5)
            ax.plot(xxin, yyin, 'k-', lw=1.5)
            ax.text(0.01, 0., 'dvp/dt', fontsize=20, horizontalalignment='left',
                                  verticalalignment='center')

        else:
            P.ion()

            vmin = - max(abs(self.coriolis.max()), abs(self.coriolis.min()))
            vmin = cut * vmin
            vmax = -vmin
            cs = N.linspace(vmin, vmax, levs)

            vmin = - max(abs(self.asVphi.max()), abs(self.asVphi.min()))
            vmin = cut * vmin
            vmax = -vmin
            cs1 = N.linspace(vmin, vmax, levs)

            print self.nvar
            for k in range(self.nvar-1): # avoid last dvp/dt which is wrong
                bal = self.asVphi[k, ...]+self.adv[k, ...]+self.rey[k, ...]+\
                      self.visc[k, ...]+self.lorentz[k, ...]+self.coriolis[k, ...]
                if k == 0:
                    ax1 = fig.add_subplot(181)
                    P.axis('off')
                    im = ax1.contourf(xx, yy, self.asVphi[k, ...], cs1, 
                                     cmap=cm, extend='both')
                    ax1.plot(xxout, yyout, 'k-', lw=1.5)
                    ax1.plot(xxin, yyin, 'k-', lw=1.5)
                    ax1.text(0.05, 0., 'uphi', fontsize=20, horizontalalignment='left',
                                  verticalalignment='center')
                    ax2 = fig.add_subplot(182)
                    P.axis('off')
                    im = ax2.contourf(xx, yy, self.adv[k, ...], cs, 
                                     cmap=cm, extend='both')
                    ax2.plot(xxout, yyout, 'k-', lw=1.5)
                    ax2.plot(xxin, yyin, 'k-', lw=1.5)
                    ax2.text(0.05, 0., 'Adv', fontsize=20, horizontalalignment='left',
                                  verticalalignment='center')
                    ax3 = fig.add_subplot(183)
                    P.axis('off')
                    im = ax3.contourf(xx, yy, self.rey[k, ...], cs, 
                                     cmap=cm, extend='both')
                    ax3.plot(xxout, yyout, 'k-', lw=1.5)
                    ax3.plot(xxin, yyin, 'k-', lw=1.5)
                    ax3.text(0.05, 0., 'Rey', fontsize=20, horizontalalignment='left',
                                  verticalalignment='center')
                    ax4 = fig.add_subplot(184)
                    P.axis('off')
                    im = ax4.contourf(xx, yy, self.visc[k, ...], cs, 
                                     cmap=cm, extend='both')
                    ax4.plot(xxout, yyout, 'k-', lw=1.5)
                    ax4.plot(xxin, yyin, 'k-', lw=1.5)
                    ax4.text(0.05, 0., 'Visc', fontsize=20, horizontalalignment='left',
                                  verticalalignment='center')
                    ax5 = fig.add_subplot(185)
                    P.axis('off')
                    im = ax5.contourf(xx, yy, self.lorentz[k, ...], cs, 
                                     cmap=cm, extend='both')
                    ax5.plot(xxout, yyout, 'k-', lw=1.5)
                    ax5.plot(xxin, yyin, 'k-', lw=1.5)
                    ax5.text(0.05, 0., 'Lo.', fontsize=20, horizontalalignment='left',
                                  verticalalignment='center')
                    ax6 = fig.add_subplot(186)
                    P.axis('off')
                    im = ax6.contourf(xx, yy, self.coriolis[k, ...], cs, 
                                     cmap=cm, extend='both')
                    ax6.plot(xxout, yyout, 'k-', lw=1.5)
                    ax6.plot(xxin, yyin, 'k-', lw=1.5)
                    ax6.text(0.05, 0., 'Cor.', fontsize=20, horizontalalignment='left',
                                  verticalalignment='center')
                    ax7 = fig.add_subplot(187)
                    P.axis('off')
                    im = ax7.contourf(xx, yy, bal, cs, 
                                     cmap=cm, extend='both')
                    ax7.plot(xxout, yyout, 'k-', lw=1.5)
                    ax7.plot(xxin, yyin, 'k-', lw=1.5)
                    ax7.text(0.05, 0., 'sum', fontsize=20, horizontalalignment='left',
                                  verticalalignment='center')
                    ax8 = fig.add_subplot(188)
                    P.axis('off')
                    im = ax8.contourf(xx, yy, self.dtVp[k, ...], cs, 
                                     cmap=cm, extend='both')
                    ax8.plot(xxout, yyout, 'k-', lw=1.5)
                    ax8.plot(xxin, yyin, 'k-', lw=1.5)
                    ax8.text(0.05, 0., 'dvp/dt', fontsize=20, horizontalalignment='left',
                                  verticalalignment='center')

                    man = P.get_current_fig_manager()
                    man.canvas.draw()
                else:
                    P.cla()
                    print k
                    im = ax1.contourf(xx, yy, self.asVphi[k, ...], cs1, 
                                     cmap=cm, extend='both')
                    im = ax2.contourf(xx, yy, self.adv[k, ...], cs, 
                                     cmap=cm, extend='both')
                    im = ax3.contourf(xx, yy, self.rey[k, ...], cs, 
                                     cmap=cm, extend='both')
                    im = ax4.contourf(xx, yy, self.visc[k, ...], cs, 
                                     cmap=cm, extend='both')
                    im = ax5.contourf(xx, yy, self.lorentz[k, ...], cs, 
                                     cmap=cm, extend='both')
                    im = ax6.contourf(xx, yy, self.coriolis[k, ...], cs, 
                                     cmap=cm, extend='both')
                    im = ax7.contourf(xx, yy, bal, cs, 
                                     cmap=cm, extend='both')
                    im = ax8.contourf(xx, yy, self.dtVp[k, ...], cs, 
                                     cmap=cm, extend='both')

                    P.axis('off')
                    man.canvas.draw()

            #P.ioff()

            



if __name__ == '__main__':
    file ='TO_mov.test'
    TOMovie(file=file)
    P.show()
