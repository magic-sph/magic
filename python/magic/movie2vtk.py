# -*- coding: utf-8 -*-
import glob
import os
import re
import numpy as np
from .npfile import *
from magic.movie import getNlines
from magic.libmagic import symmetrize

try:
    try: # Version 2 changed naming convention of functions
        #import evtk
        from evtk.hl import structuredToVTK
        #gridToVTK = evtk.hl.structuredToVTK
        gridToVTK = structuredToVTK
    except:
        import evtk
        gridToVTK = evtk.hl.gridToVTK
except:
    print("movie2vtk requires the use of evtk library!")
    print("You can get it from https://github.com/paulo-herrera/PyEVTK")

class Movie2Vtk:
    """
    This class allows to transform an input Movie file
    :ref:`movie files <secMovieFile>` into a series of vts files that
    paraview will be able to render.

    >>> m = Movie2Vtk()
    >>> # This returns a list of the available movie files in the
    >>> # current working directory
    """

    def __init__(self, file=None, step=1, lastvar=None, nvar='all',
                 fluct=False, normRad=False, precision=np.float32,
                 deminc=True, ifield=0, dir='movie2vts', store_idx=0,
                 rmin=-1., rmax=-1., closePoles=False):
        """
        :param file: the name of the movie file one wants to load
        :type file: str
        :param nvar: the number of timesteps of the movie file we want to plot
                     starting from the last line
        :type nvar: int
        :param lastvar: the number of the last timesteps to be read
        :type lastvar: int
        :param step: the stepping between two timesteps
        :type step: int
        :param fluct: if fluct=True, substract the axisymmetric part
        :type fluct: bool
        :param normRad: if normRad=True, then we normalise for each radial
                        level
        :type normRad: bool
        :param precision: precision of the input file, np.float32 for single
                          precision, np.float64 for double precision
        :type precision: str
        :param deminc: a logical to indicate if one wants do get rid of the
                       possible azimuthal symmetry
        :type deminc: bool
        :param ifield: in case of a multiple-field movie file, you can change
                       the default field displayed using the parameter ifield
        :type ifield: int
        :param store_idx: starting index for storage
        :type store_idx: int
        :param rmin: minimum radial level
        :type rmin: float
        :param rmax: maximum radial level
        :type rmax: float
        :param closePoles: when set to True, the colatitudes are replaced by
                           a linspace between 0 and pi, this allows to close
                           the little holes at the poles in paraview
        :type closePoles: bool
        """

        self.rmin = rmin
        self.rmax = rmax

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

        fieldName = filename.split('_')[0]

        ff = filename.split('_')
        if len(ff) > 2:
            pPattern = re.compile(r'P=([\+\-]?[0-9]+)_([0-9]+)')
            rPattern = re.compile(r'R=([\+\-]?[0-9]+)_([0-9]+)')
            tPattern = re.compile(r'T=([\+\-]?[0-9]+)_([0-9]+)')
            if pPattern.match(ff[1]+'_'+ff[2]):
                a = pPattern.search(ff[1]+'_'+ff[2]).groups()[0]
                b = pPattern.search(ff[1]+'_'+ff[2]).groups()[1]
                self.phiCut = float(a+'.'+b)
            elif rPattern.match(ff[1]+'_'+ff[2]):
                a = rPattern.search(ff[1]+'_'+ff[2]).groups()[0]
                b = rPattern.search(ff[1]+'_'+ff[2]).groups()[1]
                self.rCut = float(a+'.'+b)
            elif tPattern.match(ff[1]+'_'+ff[2]):
                a = tPattern.search(ff[1]+'_'+ff[2]).groups()[0]
                b = tPattern.search(ff[1]+'_'+ff[2]).groups()[1]
                self.thetaCut = float(a+'.'+b)

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
        self.radius *= rout
        self.radius_ic *= rout
        # rin = self.radratio/(1.-self.radratio)
        # self.radius = chebgrid(self.n_r_max-1, rout, rin)
        self.theta = infile.fort_read(precision)
        if closePoles:
            self.theta = np.linspace(0., np.pi, self.n_theta_max)
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
        self.var2 = int(self.var2)

        if n_surface == 0:
            self.surftype = '3d volume'
            if self.movtype in [1, 2, 3]:
                shape = (n_r_mov_tot+2, self.n_theta_max, self.n_phi_tot)
            else:
                shape = (self.n_r_max, self.n_theta_max, self.n_phi_tot)
        elif n_surface == 1:
            self.surftype = 'r_constant'
            shape = (self.n_theta_max, self.n_phi_tot)
        elif n_surface == 2:
            self.surftype = 'theta_constant'
            if self.movtype in [1, 2, 3, 14]:  # read inner core
                shape = (n_r_mov_tot+2, self.n_phi_tot)
            else:
                shape = (self.n_r_max, self.n_phi_tot)
        elif n_surface == 3:
            self.surftype = 'phi_constant'
            if self.movtype in [1, 2, 3, 14]:  # read inner core
                shape = (n_r_mov_tot+2, self.n_theta_max)
                self.n_theta_plot = 2*self.n_theta_max
            elif self.movtype in [8, 9]:
                shape = (n_r_mov_tot+2, self.n_theta_max)
                self.n_theta_plot = self.n_theta_max
            elif self.movtype in [4, 5, 6, 7, 16, 17, 18, 47, 54]:
                shape = (self.n_r_max, self.n_theta_max)
                self.n_theta_plot = 2*self.n_theta_max
            elif self.movtype in [10, 11, 12, 19, 92, 94, 95]:
                shape = (self.n_r_max, self.n_theta_max)
                self.n_theta_plot = self.n_theta_max

        if not os.path.exists(dir):
            os.mkdir(dir)

        # Read the data and store it into series of vts files

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
            for ll in range(self.n_fields):
                dat = infile.fort_read(precision)
                if n_surface == 0:
                    dat = dat.reshape(shape)
                    fname = '{}{}{}_3D_{:05d}'.format(dir, os.sep, fieldName, k+1+store_idx)
                    if self.movtype in [1, 2, 3]:
                        # datic = dat[self.n_r_max:, ...].T
                        dat = dat[:self.n_r_max, ...].T
                        self.scal3D2vtk(fname, dat, fieldName)
                    else:
                        self.scal3D2vtk(fname, dat.T, fieldName)
                elif n_surface == 2:
                    fname = '{}{}{}_eq_{:05d}'.format(dir, os.sep, fieldName, k+1+store_idx)
                    dat = dat.reshape(shape)
                    if self.movtype in [1, 2, 3, 14]:
                        # datic = dat[self.n_r_max:, :].T
                        dat = dat[:self.n_r_max, :].T
                        self.equat2vtk(fname, dat, fieldName)
                    else:
                        self.equat2vtk(fname, dat.T, fieldName)
                elif n_surface == 3:
                    if self.movtype in [1, 2, 3, 14]:
                        len1 = (self.n_r_max*self.n_theta_max*2)
                        datoc = dat[:len1]
                        # datic = dat[len1:]

                        datoc0 = datoc[:len(datoc)//2].reshape(self.n_r_max,
                                                               self.n_theta_max)
                        datoc1 = datoc[len(datoc)//2:].reshape(self.n_r_max,
                                                               self.n_theta_max)
                        fname = '{}{}{}_pcut{}_{:05d}'.format(dir, os.sep, fieldName,
                                                              str(self.phiCut), k+1+store_idx)
                        self.mer2vtk(fname, datoc0.T, self.phiCut, fieldName)
                        name = str(self.phiCut+np.pi)
                        if len(name) > 8:
                            name = name[:8]
                        fname = '{}{}{}_pcut{}_{:05d}'.format(dir, os.sep, fieldName,
                                                              name, k+1, k+1+store_idx)
                        self.mer2vtk(fname, datoc1.T, self.phiCut+np.pi,
                                     fieldName)
                        # datic0 = datic[:len(datic)/2].reshape(self.n_r_ic_max+2,
                        #                                       self.n_theta_max)
                        # datic1 = datic[len(datic)/2:].reshape(self.n_r_ic_max+2,
                        #                                       self.n_theta_max)
                    elif self.movtype in [4, 5, 6, 7, 16, 17, 18, 47, 54, 91]:
                        dat0 = dat[:len(dat)//2].reshape(shape)
                        dat1 = dat[len(dat)//2:].reshape(shape)
                        fname = '{}{}{}_pcut{}_{:05d}'.format(dir, os.sep, fieldName,
                                                              str(self.phiCut), k+1+store_idx)
                        self.mer2vtk(fname, dat0.T, self.phiCut, fieldName)
                        name = str(self.phiCut+np.pi)
                        if len(name) > 8:
                            name = name[:8]
                        fname = '{}{}{}_pcut{}_{:05d}'.format(dir, os.sep, fieldName,
                                                              name, k+1+store_idx)
                        self.mer2vtk(fname, dat1.T, self.phiCut+np.pi,
                                     fieldName)
                else:
                    fname = '{}{}{}_rcut{}_{:05d}'.format(dir, os.sep, fieldName,
                                                          str(self.rCut), k+1+store_idx)
                    dat = dat.reshape(shape)
                    self.rcut2vtk(fname, dat.T, self.rCut, fieldName)

        infile.close()

    def scal3D2vtk(self, fname, data, name):
        """
        This routine transforms a 3-D scalar field of dimension
        (n_phi,n_theta,n_r) into a vts file.

        :param fname: file name
        :type fname: str
        :param data: input data
        :type data: numpy.ndarray
        :param r: radius of the selected cut
        :type r: float
        :param name: name of the physical field stored in the vts file
        :type name: str
        """
        if self.rmin > 0 and self.rmax > 0:
            mask = np.where(abs(self.radius-self.rmin)==abs(self.radius-self.rmin).min(), 1, 0)
            idx2 = np.nonzero(mask)[0][0]
            mask = np.where(abs(self.radius-self.rmax)==abs(self.radius-self.rmax).min(), 1, 0)
            idx1 = np.nonzero(mask)[0][0]
            nr = idx2-idx1+1
        else:
            nr = data.shape[1]
            idx1 = 0
            idx2 = data.shape[1]

        data = symmetrize(data, self.minc)
        phi = np.linspace(0., 2.*np.pi, data.shape[0])
        X = np.zeros(data.shape, dtype=np.float32)
        Y = np.zeros_like(X)
        Z = np.zeros_like(X)
        ttheta, pphi, rr = np.meshgrid(self.theta, phi, self.radius[idx1:idx2+1])
        X = rr*np.sin(ttheta)*np.cos(pphi)
        Y = rr*np.sin(ttheta)*np.sin(pphi)
        Z = rr*np.cos(ttheta)*np.ones_like(pphi)

        point_data = {}
        point_data[name] = data[..., idx1:idx2]

        gridToVTK(fname, X, Y, Z, pointData=point_data)
        print('Store {}.vts'.format(fname))

    def rcut2vtk(self, fname, data, r, name):
        """
        This routine transforms a radial cut of dimension (n_phi,n_theta)
        into a vts file.

        :param fname: file name
        :type fname: str
        :param data: input data
        :type data: numpy.ndarray
        :param r: radius of the selected cut
        :type r: float
        :param name: name of the physical field stored in the vts file
        :type name: str
        """
        data = symmetrize(data, self.minc)
        phi = np.linspace(0., 2.*np.pi, data.shape[0])
        X = np.zeros((data.shape[0], data.shape[1], 1), dtype=np.float32)
        Y = np.zeros_like(X)
        Z = np.zeros_like(X)
        ttheta, pphi = np.meshgrid(self.theta, phi)
        X[:, :, 0] = r*np.sin(ttheta)*np.cos(pphi)
        Y[:, :, 0] = r*np.sin(ttheta)*np.sin(pphi)
        Z[:, :, 0] = r*np.cos(ttheta)

        dat = np.zeros_like(X)
        point_data = {}
        dat[..., 0] = data
        point_data[name] = dat

        gridToVTK(fname, X, Y, Z, pointData=point_data)
        print('Store {}.vts'.format(fname))

    def mer2vtk(self, fname, data, phi0, name):
        """
        This routine transforms a meridional cut of dimension (n_theta,n_r)
        into a vts file.

        :param fname: file name
        :type fname: str
        :param data: input data
        :type data: numpy.ndarray
        :param phi0: longitude of the selected cut
        :type phi0: float
        :param name: name of the physical field stored in the vts file
        :type name: str
        """
        if self.rmin > 0 and self.rmax > 0:
            mask = np.where(abs(self.radius-self.rmin)==abs(self.radius-self.rmin).min(), 1, 0)
            idx2 = np.nonzero(mask)[0][0]
            mask = np.where(abs(self.radius-self.rmax)==abs(self.radius-self.rmax).min(), 1, 0)
            idx1 = np.nonzero(mask)[0][0]
            nr = idx2-idx1+1
        else:
            nr = data.shape[1]
            idx1 = 0
            idx2 = data.shape[1]

        X = np.zeros((1, data.shape[0], nr), dtype=np.float32)
        Y = np.zeros_like(X)
        Z = np.zeros_like(X)
        rr, ttheta = np.meshgrid(self.radius[idx1:idx2+1], self.theta)

        X[0, :, :] = rr*np.sin(ttheta)*np.cos(phi0)
        Y[0, :, :] = rr*np.sin(ttheta)*np.sin(phi0)
        Z[0, :, :] = rr*np.cos(ttheta)

        dat = np.zeros_like(X)
        dat[0, ...] = data[:, idx1:idx2+1]

        point_data = {}
        point_data[name] = dat

        gridToVTK(fname, X, Y, Z, pointData=point_data)
        print('Store {}.vts'.format(fname))

    def equat2vtk(self, fname, data, name):
        """
        This routine transforms an equatorial cut of dimension (n_phi,n_r)
        into a vts file.

        :param fname: file name
        :type fname: str
        :param data: input data
        :type data: numpy.ndarray
        :param name: name of the physical field stored in the vts file
        :type name: str
        """
        if self.rmin > 0 and self.rmax > 0:
            mask = np.where(abs(self.radius-self.rmin)==abs(self.radius-self.rmin).min(), 1, 0)
            idx2 = np.nonzero(mask)[0][0]
            mask = np.where(abs(self.radius-self.rmax)==abs(self.radius-self.rmax).min(), 1, 0)
            idx1 = np.nonzero(mask)[0][0]
            nr = idx2-idx1+1
        else:
            nr = data.shape[1]
            idx1 = 0
            idx2 = data.shape[1]
        data = symmetrize(data, self.minc)
        phi = np.linspace(0., 2.*np.pi, data.shape[0])
        X = np.zeros((1, data.shape[0], nr), dtype=np.float32)
        Y = np.zeros_like(X)
        Z = np.zeros_like(X)
        rr, pphi = np.meshgrid(self.radius[idx1:idx2+1], phi)
        X[0, :, :] = rr*np.cos(pphi)
        Y[0, :, :] = rr*np.sin(pphi)

        dat = np.zeros_like(X)
        dat[0, ...] = data[:, idx1:idx2+1]

        point_data = {}
        point_data[name] = dat

        gridToVTK(fname, X, Y, Z, pointData=point_data)
        print('Store {}.vts'.format(fname))
