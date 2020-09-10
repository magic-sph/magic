# -*- coding: utf-8 -*-
import glob
import re
import os
import numpy as np
from magic.libmagic import symmetrize
from magic.setup import buildSo
from magic import ExtraPot
from .npfile import *
import sys
if buildSo:
    from .vtklib import *


class Movie3D:
    """
    This class allows to read the 3D movie files :ref:`(B|V)_3D_.TAG<secMovieFile>`  and
    transform them into a series of VTS files ``./vtsFiles/B3D_#.TAG`` that can be further
    read using paraview.

    >>> Movie3D(file='B_3D.TAG')
    """

    def __init__(self, file=None, step=1, lastvar=None, nvar='all', nrout=48,
                 ratio_out=2., potExtra=False, precision=np.float32):
        """
        :param file: file name
        :type file: str
        :param nvar: the number of timesteps of the movie file we want to plot
                     starting from the last line
        :type nvar: int
        :param lastvar: the number of the last timestep to be read
        :type lastvar: int
        :param step: the stepping between two timesteps
        :type step: int
        :param precision: precision of the input file, np.float32 for single
                          precision, np.float64 for double precision
        :type precision: str
        :param potExtra: when set to True, potential extrapolation of the
                         magnetic field outside the fluid domain is also
                         computed
        :type potExtra: bool
        :param ratio_out: ratio of desired external radius to the CMB radius.
                          This is is only used when potExtra=True
        :type ratio_out: float
        :param nrout: number of additional radial grid points to compute the
                      potential extrapolation. This is only used when
                      potExtra=True
        :type nrout: int
        """
        if file is None:
            dat = glob.glob('*_mov.*')
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

        vecNames = np.r_[3]
        scalNames = np.r_[-1]

        # READ the movie file
        infile = npfile(filename, endian='B')
        # HEADER
        version = infile.fort_read('|S64')
        n_type, n_surface, const, n_fields = infile.fort_read(precision)
        n_fields = int(n_fields)
        n_surface = int(n_surface)
        if n_fields == 1:
            movtype = infile.fort_read(precision)
            self.movtype = int(movtype)
        else:
            movtype = infile.fort_read(precision)

        # RUN PARAMETERS
        runid = infile.fort_read('|S64')
        n_r_mov_tot, n_r_max, n_theta_max, n_phi_tot, minc, ra, \
            ek, pr, prmag, radratio, tScale = infile.fort_read(precision)
        minc = int(minc)
        self.n_r_max = int(n_r_max)
        self.n_theta_max = int(n_theta_max)
        self.n_phi_tot = int(n_phi_tot)
        n_r_mov_tot = int(n_r_mov_tot)

        # GRID
        if potExtra:
            self.radius = np.zeros((n_r_mov_tot+1+nrout), np.float32)
        else:
            self.radius = np.zeros((n_r_mov_tot+2), np.float32)
        tmp = infile.fort_read(precision)/(1.-radratio)
        rcmb = tmp[0]
        self.radius[2:n_r_mov_tot+2] = tmp[::-1]
        self.theta = infile.fort_read(precision)
        self.phi = infile.fort_read(precision)

        shape = (n_r_mov_tot+2, self.n_theta_max, self.n_phi_tot)

        self.time = np.zeros(self.nvar, precision)
        if potExtra:
            scals = np.zeros((n_r_mov_tot+1+nrout, self.n_theta_max,
                              self.n_phi_tot*minc+1, 1), np.float32)
        else:
            scals = np.zeros((n_r_mov_tot+2, self.n_theta_max,
                              self.n_phi_tot*minc+1, 1), np.float32)
        vecr = np.zeros_like(scals)
        vect = np.zeros_like(scals)
        vecp = np.zeros_like(scals)

        # Potential extrapolation
        radii = self.radius
        scals = scals.T
        scals[0, ...] = radii

        startdir = os.getcwd()
        if not os.path.exists('vtsFiles'):
            os.mkdir('vtsFiles')
            os.chdir('vtsFiles')
        else:
            os.chdir('vtsFiles')
        for i in range(self.var2-self.nvar):
            n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                 movieDipLon, movieDipStrength, \
                 movieDipStrengthGeo = infile.fort_read(precision)
            tmp = infile.fort_read(precision, shape=shape)
            vecr[0:n_r_mov_tot+2, :, :, 0] = symmetrize(tmp, minc,
                                                        reversed=True)
            tmp = infile.fort_read(precision, shape=shape)
            vect[0:n_r_mov_tot+2, :, :, 0] = symmetrize(tmp, minc,
                                                        reversed=True)
            tmp = infile.fort_read(precision, shape=shape)
            vecp[0:n_r_mov_tot+2, :, :, 0] = symmetrize(tmp, minc,
                                                        reversed=True)
        for k in range(self.nvar):
            n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                 movieDipLon, movieDipStrength, \
                 movieDipStrengthGeo = infile.fort_read(precision)
            self.time[k] = t_movieS
            if k % step == 0:
                # print(k+self.var2-self.nvar)
                tmp = infile.fort_read(precision, shape=shape)
                brCMB = np.zeros((self.n_theta_max, self.n_phi_tot), np.float32)
                brCMB = tmp[0, :, :]
                brCMB = brCMB.T
                if potExtra:
                    vecr[nrout-1:nrout+n_r_mov_tot+2, :, :, 0] = \
                        symmetrize(tmp, minc, reversed=True)
                    tmp = infile.fort_read(precision, shape=shape)
                    vect[nrout-1:nrout+n_r_mov_tot+2, :, :, 0] = \
                        symmetrize(tmp, minc, reversed=True)
                    tmp = infile.fort_read(precision, shape=shape)
                    vecp[nrout-1:nrout+n_r_mov_tot+2, :, :, 0] = \
                        symmetrize(tmp, minc, reversed=True)
                else:
                    vecr[0:n_r_mov_tot+2, :, :, 0] = \
                        symmetrize(tmp, minc, reversed=True)
                    tmp = infile.fort_read(precision, shape=shape)
                    vect[0:n_r_mov_tot+2, :, :, 0] = \
                        symmetrize(tmp, minc, reversed=True)
                    tmp = infile.fort_read(precision, shape=shape)
                    vecp[0:n_r_mov_tot+2, :, :, 0] = \
                        symmetrize(tmp, minc, reversed=True)
                filename = 'B3D_{:05d}'.format(k)
                vecr = vecr[::-1, ...]
                vect = vect[::-1, ...]
                vecp = vecp[::-1, ...]
                br = vecr.T
                bt = vect.T
                bp = vecp.T
                if potExtra:
                    pot = ExtraPot(rcmb, brCMB, minc, ratio_out=ratio_out,
                                   nrout=nrout, cutCMB=True, deminc=True)
                    br[0, ..., n_r_mov_tot+2:] = pot.brout
                    bt[0, ..., n_r_mov_tot+2:] = pot.btout
                    bp[0, ..., n_r_mov_tot+2:] = pot.bpout
                    radii[n_r_mov_tot+2:] = pot.rout
                else:
                    radii = self.radius
                vts(filename, radii, br, bt, bp, scals, scalNames, vecNames, 1)
                print('write {}.vts'.format(filename))
            else:  # Otherwise we read
                vecr = infile.fort_read(precision, shape=shape)
                vect = infile.fort_read(precision, shape=shape)
                vecp = infile.fort_read(precision, shape=shape)

        os.chdir(startdir)


if __name__ == '__main__':
    t1 = Movie3D(file='B_3D_mov.CJ2')
