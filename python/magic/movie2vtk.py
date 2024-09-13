# -*- coding: utf-8 -*-
import glob
import os
import re
import numpy as np
from .npfile import npfile
from magic.movie import getNlines, Movie
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

class Movie2Vtk(Movie):
    """
    This class allows to transform an input Movie file
    :ref:`movie files <secMovieFile>` into a series of vts files that
    paraview will be able to render.

    >>> m = Movie2Vtk()
    >>> # This returns a list of the available movie files in the
    >>> # current working directory
    """

    def __init__(self, file=None, step=1, lastvar=None, nvar='all', fluct=False,
                 normRad=False, precision=np.float32, deminc=True, ifield=0,
                 dir='movie2vts', store_idx=0, datadir='.', rmin=-1., rmax=-1.,
                 closePoles=False, mean_field=False):
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
        :param mean_field: if mean_field=True, only retain the axisymmetric part
                           (this necessitates the storage of the relevant axisymmetric
                           movie file for the phi-slices)
        :type mean_field: bool
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
        :param datadir: working directory
        :type datadir : str
        """

        if mean_field: # Cannot be both at the same time!
            fluct = False
        if fluct:
            mean_field = False

        self.rmin = rmin
        self.rmax = rmax

        if file is None:
            dat = glob.glob(os.path.join(datadir, '*[Mm]ov.*'))
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
        end = mot.findall(os.path.join(datadir, filename))[0]

        fieldName = filename.split('_')[0]

        ff = filename.split('_')
        field = ff[0]
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
        infile = npfile(os.path.join(datadir, filename), endian='B')

        # Header
        self._read_header(infile, ifield, precision)
        if closePoles:
            self.theta = np.linspace(0., np.pi, self.n_theta_max)

        # Get the number of lines
        self._get_nlines(datadir, filename, nvar, lastar)

        # Determine the shape of the data
        self._get_data_shape(precision, allocate=False)

        if self.n_surface == 3:
            if fluct or mean_field:
                if 'B' in field:
                    prefix = 'AB_mov'
                elif 'V' in field:
                    prefix = 'AV_mov'
                elif 'T' in field:
                    prefix = 'ATmov'
                elif 'C' in field:
                    prefix = 'ACmov'
                tag = filename.split('.')[-1]
                av_mov = prefix + '.' + tag
                mov_mean = Movie(file=av_mov, datadir=datadir, iplot=False)

        if not os.path.exists(os.path.join(datadir, dir)):
            os.mkdir(os.path.join(datadir, dir))

        # Read the data and store it into series of vts files

        # If one skip the beginning, nevertheless read but do not store
        for i in range(self.var2-self.nvar):
            if self.version == 2:
                n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                    movieDipLon, movieDipStrength, movieDipStrengthGeo = \
                infile.fort_read(precision)
            else:
                t_movieS = infile.fort_read(precision)
            for ll in range(self.n_fields):
                dat = infile.fort_read(precision, shape=self.shape, order='F')
        # then read the remaining requested nvar lines
        for k in range(self.nvar):
            if self.version == 2:
                n_frame, t_movieS, omega_ic, omega_ma, movieDipColat, \
                    movieDipLon, movieDipStrength, movieDipStrengthGeo = \
                    infile.fort_read(precision)
            else:
                t_movieS = infile.fort_read(precision)
            for ll in range(self.n_fields):
                dat = infile.fort_read(precision, shape=self.shape, order='F')
                if self.n_surface == 0:
                    fname = '{}{}{}_3D_{:05d}'.format(dir, os.sep, fieldName,
                                                      k+1+store_idx)
                    dat = dat[:, :, :self.n_r_max]
                    if fluct:
                        dat -= dat.mean(axis=0)
                    elif mean_field:
                        tmp = dat.mean(axis=0)
                        for ip in range(self.n_phi_tot):
                            dat[ip, :, :] = tmp[:, :]
                    fname = os.path.join(datadir, fname)
                    self.scal3D2vtk(fname, dat, fieldName)
                elif self.n_surface == 2:
                    fname = '{}{}{}_eq_{:05d}'.format(dir, os.sep, fieldName,
                                                      k+1+store_idx)
                    dat = dat[:, :self.n_r_max]
                    if fluct:
                        dat -= dat.mean(axis=0)
                    elif mean_field:
                        tmp = np.zeros_like(dat)
                        for ip in range(self.n_phi_tot):
                            tmp[ip, :] = dat.mean(axis=0)
                        dat = tmp
                    fname = os.path.join(datadir, fname)
                    self.equat2vtk(fname, dat, fieldName)
                elif self.n_surface == 3:
                    if fluct or mean_field:
                        field_m = mov_mean.data[0, k, ...]
                    datoc0 = dat[:, :self.n_r_max]
                    datoc1 = dat[:, self.n_r_max:2*self.n_r_max]
                    if fluct:
                        datoc0 -= field_m
                        datoc1 -= field_m
                    elif mean_field:
                        datoc0 = field_m
                        datoc1 = field_m

                    fname = '{}{}{}_pcut{}_{:05d}'.format(dir, os.sep,
                                                          fieldName,
                                                          str(self.phiCut),
                                                          k+1+store_idx)
                    fname = os.path.join(datadir, fname)
                    self.mer2vtk(fname, datoc0, self.phiCut, fieldName)
                    name = str(self.phiCut+np.pi)
                    if len(name) > 8:
                        name = name[:8]
                    fname = '{}{}{}_pcut{}_{:05d}'.format(dir, os.sep,
                                                          fieldName,
                                                          name, k+1,
                                                          k+1+store_idx)
                    fname = os.path.join(datadir, fname)
                    self.mer2vtk(fname, datoc1, self.phiCut+np.pi,
                                 fieldName)
                else:
                    fname = '{}{}{}_rcut{}_{:05d}'.format(dir, os.sep,
                                                          fieldName,
                                                          str(self.rCut),
                                                          k+1+store_idx)
                    fname = os.path.join(datadir, fname)
                    if fluct:
                        dat -= dat.mean(axis=0)
                    elif mean_field:
                        tmp = dat.mean(axis=0)
                        for ip in range(self.n_phi_tot):
                            dat[ip, :] = tmp
                    self.rcut2vtk(fname, dat, self.rCut, fieldName)

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
        if self.rmin >= 0 and self.rmax > 0:
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
        if self.rmin >= 0 and self.rmax > 0:
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
        if self.rmin >= 0 and self.rmax > 0:
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
