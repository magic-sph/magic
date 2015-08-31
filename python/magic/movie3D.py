# -*- coding: utf-8 -*-
import glob
import re,os
import numpy as N
from magic.libmagic import symmetrize
from magic import ExtraPot
from .npfile import *
import sys
if sys.version_info.major == 3:
    from vtklib3 import *
elif  sys.version_info.major == 2:
    from vtklib2 import *


__author__  = "$Author$"
__date__   = "$Date$"
__version__ = "$Revision$"


class Movie3D:
    """
    This class allows to transform a B_3D.TAG file into a series of
    vts files
    """

    def __init__(self, file=None, step=1, lastvar=None, nvar='all', nrout=48,
                 ratio_out=2., extrapot=False):
        """
        :param nvar: the number of lines of the movie file we want to plot
                     starting from the last line
        :param lastvar: the rank of the last line to be read
        :param step: the stepping between two lines             
        """
        if file == None:
            dat = glob.glob('*_mov.*')
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
        n_fields = int(n_fields)
        n_surface = int(n_surface)
        if n_fields == 1:
            movtype = infile.fort_read('f')
            self.movtype = int(movtype)
        else:
            movtype = infile.fort_read('f')

        # RUN PARAMETERS
        runid = infile.fort_read('|S64')
        n_r_mov_tot, n_r_max, n_theta_max, n_phi_tot, minc, ra, \
             ek, pr, prmag, radratio, tScale = infile.fort_read('f')
        minc = int(minc)
        self.n_r_max = int(n_r_max)
        self.n_theta_max = int(n_theta_max)
        self.n_phi_tot = int(n_phi_tot)
        n_r_mov_tot = int(n_r_mov_tot)

        # GRID
        self.radius = infile.fort_read('f')
        self.radius = self.radius[:self.n_r_max] # remove inner core
        self.radius = self.radius[::-1]/(1.-radratio)
        self.theta = infile.fort_read('f')
        self.phi = infile.fort_read('f')

        shape = (n_r_mov_tot+2, self.n_theta_max, self.n_phi_tot)

        self.time = N.zeros(self.nvar, 'f')

        if not os.path.exists('vtsFiles'):
            startdir = os.getcwd()
            os.mkdir('vtsFiles')
            os.chdir('vtsFiles')
        for i in range(self.var2-self.nvar):
            n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                 movieDipLon, movieDipStrength, \
                 movieDipStrengthGeo = infile.fort_read('f')
            vecr = infile.fort_read('f', shape=shape)
            vect = infile.fort_read('f', shape=shape)
            vecp = infile.fort_read('f', shape=shape)
        for k in range(self.nvar):
            n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                 movieDipLon, movieDipStrength, \
                 movieDipStrengthGeo = infile.fort_read('f')
            self.time[k] = t_movieS
            if k % step == 0:
                #print(k+self.var2-self.nvar)
                vecr = infile.fort_read('f', shape=shape)
                vect = infile.fort_read('f', shape=shape)
                vecp = infile.fort_read('f', shape=shape)
                filename = 'B3D_%05d' % k
                vecr = vecr[:self.n_r_max, ...] # remove inner core 
                vecr = vecr[::-1, ...]
                vect = vect[:self.n_r_max, ...] # remove inner core 
                vect = vect[::-1, ...]
                vecp = vecp[:self.n_r_max, ...] # remove inner core 
                vecp = vecp[::-1, ...]
                br = vecr.T
                bt = vect.T
                bp = vecp.T
                brCMB = br[..., -1]
                rcmb = self.radius[-1]
                if extrapot:
                    pot = ExtraPot(rcmb, brCMB, minc, ratio_out=ratio_out, 
                                   nrout=nrout, cutCMB=True)

                    br = symmetrize(br, minc)
                    bt = symmetrize(bt, minc)
                    bp = symmetrize(bp, minc)
                    radii = N.concatenate((self.radius, pot.rout))
                    br = N.concatenate((br, pot.brout), axis=2)
                    bt = N.concatenate((bt, pot.btout), axis=2)
                    bp = N.concatenate((bp, pot.bpout), axis=2)
                else:
                    radii = self.radius
                vts(filename, radii, br, bt, bp)
                print('write %s.vts' % filename)
            else: # Otherwise we read
                vecr = infile.fort_read('f', shape=shape)
                vect = infile.fort_read('f', shape=shape)
                vecp = infile.fort_read('f', shape=shape)

        os.chdir(startdir)


if __name__ == '__main__':
    from magic import MagicGraph

    t1 = Movie3D(file='B_3D_mov.CJ2')

