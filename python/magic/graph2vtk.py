# -*- coding: utf-8 -*-
from scipy.ndimage import map_coordinates
from scipy.interpolate import interp1d
from magic import MagicGraph, BLayers, ExtraPot
from .setup import labTex
from .libmagic import symmetrize, thetaderavg, rderavg, phideravg
import matplotlib.pyplot as plt
import numpy as np
from .vtklib import *

def sph2cart_scal(scals, radius, nx=96, ny=96, nz=96, minc=1):
    """
    This function interpolates a series of scalar fields from the spherical
    coordinates to the cartesian coordinates.

    :param scals: an array that contains the different scalar quantities
    :type scals: numpy.ndarray[nscals,nphi,ntheta,nr]
    :param radius: the input radius
    :type radius: numpy.ndarray
    :param nx: number of grid points in the x direction
    :type nx: int
    :param ny: number of grid points in the x direction
    :type ny: int
    :param nz: number of grid points in the x direction
    :type nz: int
    :param minc: azimuthal symmetry
    :type minc: int
    :returns: a tuple that contains the scalars, the max of the grid
              and the grid spacing
    :rtype: (numpy.ndarray[nscals,nz,ny,nx],float,float)
    """
    nscals, nphi, nt, nr = scals.shape
    phi = np.linspace(-np.pi/minc, np.pi/minc, nphi)
    theta = np.linspace(0., np.pi, nt)
    # Cube: take care of the sqrt(3.) !!!
    gridMax = radius.max()
    spacing = 2.*gridMax/(nx-1)
    Z,Y,X = np.mgrid[-1:1:nz*1j,-1:1:ny*1j,-1:1:nx*1j]*gridMax
    new_r = np.sqrt(X**2+Y**2+Z**2)#.ravel()
    new_phi = np.arctan2(Y, X)#.ravel()
    new_theta = np.arctan2(np.sqrt(X**2+Y**2), Z)#.ravel()
    del X,Y,Z

    ir = interp1d(radius, np.arange(len(radius)), bounds_error=False)
    it = interp1d(theta, np.arange(len(theta)), bounds_error=False)
    ip = interp1d(phi, np.arange(len(phi)), bounds_error=False)

    new_ir = ir(new_r)
    new_it = it(new_theta)
    new_ip = ip(new_phi)

    new_ir[new_r < radius.min()] = 0.
    new_it[new_r < radius.min()] = 0.
    new_it[new_r < radius.min()] = 0.

    coords = np.array([new_ip, new_it, new_ir])

    scals_cart = np.zeros((nscals, nz, ny, nx), 'f')
    for iscal in range(nscals):
        if iscal == 0: # radius has been already calculated
            scals_cart[iscal, ...] = new_r
        else:
            scals_cart[iscal, ...] = map_coordinates(scals[iscal, ...], coords)
            scals_cart[iscal, new_r < radius.min()] = 0.
            scals_cart[iscal, new_r > radius.max()] = 0.
    del coords

    del new_theta, new_phi

    return scals_cart,gridMax,spacing

def sph2cart_vec(vecr, vect, vecp, radius, nx=96, ny=96, nz=96, minc=1):
    """
    This function interpolates a series of vector fields from the spherical
    coordinates to the cartesian coordinates.

    :param vecr: the radial components of the different vector fields
    :type vecr: numpy.ndarray[nvecs,nphi,ntheta,nr]
    :param vect: the latitudinal components of the different vector fields
    :type vect: numpy.ndarray[nvecs,nphi,ntheta,nr]
    :param vecp: the azimuthal components of the different vector fields
    :type vecp: numpy.ndarray[nvecs,nphi,ntheta,nr]
    :param radius: the input radius
    :type radius: numpy.ndarray
    :param nx: number of grid points in the x direction
    :type nx: int
    :param ny: number of grid points in the x direction
    :type ny: int
    :param nz: number of grid points in the x direction
    :type nz: int
    :param minc: azimuthal symmetry
    :type minc: int
    :returns: a tuple that contains the three vectors components
    :rtype: (numpy.ndarray[nvecs,nz,ny,nx],...)
    """
    nvecs, nphi, nt, nr = vecr.shape
    phi = np.linspace(-np.pi/minc, np.pi/minc, nphi)
    theta = np.linspace(0., np.pi, nt)
    # Cube: take care of the sqrt(3.) !!!
    gridMax = radius.max()
    spacing = 2.*gridMax/(nx-1)
    Z,Y,X = np.mgrid[-1:1:nz*1j,-1:1:ny*1j,-1:1:nx*1j]*gridMax
    new_r = np.sqrt(X**2+Y**2+Z**2)#.ravel()
    new_phi = np.arctan2(Y, X)#.ravel()
    new_theta = np.arctan2(np.sqrt(X**2+Y**2), Z)#.ravel()
    del X,Y,Z

    ir = interp1d(radius, np.arange(len(radius)), bounds_error=False)
    it = interp1d(theta, np.arange(len(theta)), bounds_error=False)
    ip = interp1d(phi, np.arange(len(phi)), bounds_error=False)

    new_ir = ir(new_r)
    new_it = it(new_theta)
    new_ip = ip(new_phi)

    new_ir[new_r < radius.min()] = 0.
    new_it[new_r < radius.min()] = 0.
    new_it[new_r < radius.min()] = 0.

    coords = np.array([new_ip, new_it, new_ir])

    vecx = np.zeros((nvecs, nz, ny, nx), 'f')
    vecy = np.zeros_like(vecx)
    vecz = np.zeros_like(vecx)
    for ivec in range(nvecs):
        br_cart = map_coordinates(vecr[ivec, ...], coords)
        bt_cart = map_coordinates(vect[ivec, ...], coords)
        bp_cart = map_coordinates(vecp[ivec, ...], coords)

        br_cart[new_r < radius.min()] = 0.
        bt_cart[new_r < radius.min()] = 0.
        bp_cart[new_r < radius.min()] = 0.

        br_cart[new_r > radius.max()] = 0.
        bt_cart[new_r > radius.max()] = 0.
        bp_cart[new_r > radius.max()] = 0.

        bs = br_cart*np.sin(new_theta)+bt_cart*np.cos(new_theta)
        vecx[ivec, ...] = bs*np.cos(new_phi) - bp_cart*np.sin(new_phi)
        vecy[ivec, ...] = bs*np.sin(new_phi) + bp_cart*np.cos(new_phi)
        vecz[ivec, ...] = br_cart*np.cos(new_theta) - bt_cart*np.sin(new_theta)

    del coords, br_cart, bp_cart

    del new_theta, new_phi

    return vecx,vecy,vecz

class Graph2Vtk:
    """
    This class allows to transform an input graphic file to a file format
    readable by paraview/visit or mayavi. It also allows to compute a possible
    potential extrapolation of the field lines in an arbitrary outer
    spherical shell domain

    >>> # Load a graphic file
    >>> gr = MagicGraph(ivar=1)
    >>> # store myOut.vts
    >>> Graph2Vtk(gr, 'myOut', outType='vts')
    >>> # store u' and B for the vector fields and vortz and T for the scalars
    >>> Graph2Vtk(gr, scals=['temp', 'vortz'], vecs=['ufluct', 'B'])
    >>> # store only T'
    >>> Graph2Vtk(gr, scals=['tempfluct'], vecs=[])
    >>> # store only B with its potential extrapolation up to 3*r_cmb
    >>> Graph2Vtk(gr, scals=[], vecs=['B'], potExtra=True, ratio_out=3)
    >>> # Extrapolate on a cartesian grid of size 128^3
    >>> Graph2Vtk(gr, outType='vti', nx=128, ny=128, nz=128)
    """

    def __init__(self, gr, filename='out', scals=['vr', 'emag', 'tfluct'],
                 vecs=['u', 'B'], potExtra=False, ratio_out=2, nrout=32,
                 deminc=True, outType='vts', nFiles=1, nx=96, ny=96, nz=96,
                 labFrame=False):
        """
        :param filename: the file name of the output (without extension)
        :type filename: str
        :param gr: the input graphic file one wants to transform to vts/vti
        :type gr: magic.MagicGraph
        :param scals: a list that contains the possible input scalars: 'entropy',
                      'vr', 'vp', 'tfluct', 'vortz', 'vortzfluct', 'ekin',
                      'emag', 'vortr', 'colat'
        :type scals: list(str)
        :param vecs: a list that contains the possible input vectors: 'u',
                     'b', 'ufluct', 'bfluct'
        :type vecs: list(str)
        :param potExtra: when set to True, calculates the potential extrapolation
                         of the magnetic field up to ratio_out*r_cmb
        :type potExtra: bool
        :param ratio_out: in case of potential extrapolation, this is the ratio
                          of the external outer radius to r_cmb (rout/rcmb)
        :type ratio_out: float
        :param nrout: in case of potential extrapolation, this input allows
                      to specify thenumber of radial grid points in the
                      outer spherical envelope
        :type nrout: integer
        :param deminc: a logical to indicate if one wants do get rid of the
                       possible azimuthal symmetry
        :type deminc: bool
        :param outType: nature of the VTK file produced. This can be either
                        'vts' for the spherical grid or 'vti' for an extrapolation
                        on a cartesian grid
        :type outType: str
        :param nFiles: number of output chunks in case of parallel vts file format
                       (pvts)
        :type nFiles: int
        :param nx: number of grid points in the x direction
        :type nx: int
        :param ny: number of grid points in the x direction
        :type ny: int
        :param nz: number of grid points in the x direction
        :type nz: int
        :param labFrame: when set to True, transform the velocity to the lab frame
        :type labFrame: bool
        """


        if deminc:
            self.minc = 1
        else:
            self.minc = gr.minc

        if labFrame:
            th3D = np.zeros_like(gr.vphi)
            rr3D = np.zeros_like(th3D)
            for i in range(gr.ntheta):
                th3D[:, i, :] = gr.colatitude[i]
            for i in range(gr.nr):
                rr3D[:, :, i] = gr.radius[i]
            s3D = rr3D * np.sin(th3D)
            gr.vphi = gr.vphi+s3D/gr.ek

        keyScal = {}
        keyScal['radius'] = -1
        keyScal['entropy'] = 1
        keyScal['temperature'] = 1
        keyScal['temp'] = 1
        keyScal['s'] = 1
        keyScal['t'] = 1
        keyScal['T'] = 1
        keyScal['S'] = 1
        keyScal['emag'] = 2
        keyScal['vortz'] = 3
        keyScal['vr'] = 4
        keyScal['ur'] = 4
        keyScal['vp'] = 5
        keyScal['vphi'] = 5
        keyScal['uphi'] = 5
        keyScal['vphi'] = 5
        keyScal['tfluct'] = 6
        keyScal['entropyfluct'] = 6
        keyScal['tempfluct'] = 6
        keyScal['vortzfluct'] = 7
        keyScal['vortr'] = 8
        keyScal['flucttemp'] = 9
        keyScal['ekin'] = 10
        keyScal['br'] = 11
        keyScal['Br'] = 11
        keyScal['vs'] = 12
        keyScal['Vs'] = 12
        keyScal['us'] = 12
        keyScal['colat'] = 13
        keyScal['theta'] = 13
        keyScal['xi'] = 14
        keyScal['Xi'] = 14
        keyScal['xifluct'] = 15
        keyScal['Xifluct'] = 15

        # Change default scalars and vectors in non-magnetic cases
        if hasattr(gr, 'mode'):
            if gr.mode == 1 or gr.mode == 7 or gr.mode == 10:
                if 'emag' in keyScal:
                    keyScal.__delitem__('emag')
                if 'br' in keyScal:
                    keyScal.__delitem__('br')
                if 'Br' in keyScal:
                    keyScal.__delitem__('Br')

        self.scalNames = np.zeros(len(scals), 'i')
        for k, scal in enumerate(scals):
            if scal in keyScal:
                self.scalNames[k] = keyScal[scal]
            else:
                self.scalNames[k] = 0
        # Include the radius as the first scalar quantity
        self.scalNames = np.insert(self.scalNames, 0, -1)
        # Remove the possible 0
        self.scalNames = self.scalNames[self.scalNames!=0]
        self.nscals = len(self.scalNames)

        keyVec = {}
        keyVec['u'] = 1
        keyVec['U'] = 1
        keyVec['v'] = 1
        keyVec['vel'] = 1
        keyVec['V'] = 1
        keyVec['ufluct'] = 2
        keyVec['vfluct'] = 2
        keyVec['Ufluct'] = 2
        keyVec['Vfluct'] = 2
        keyVec['b'] = 3
        keyVec['B'] = 3
        keyVec['mag'] = 3
        keyVec['bfluct'] = 4
        keyVec['Bfluct'] = 4

        if hasattr(gr, 'mode'):
            if gr.mode == 1 or gr.mode == 7 or gr.mode == 10:
                if 'bfluc' in keyVec:
                    keyVec.__delitem__('bfluct')
                if 'Bfluct' in keyVec:
                    keyVec.__delitem__('Bfluct')
                if 'b' in keyVec:
                    keyVec.__delitem__('b')
                if 'B' in keyVec:
                    keyVec.__delitem__('B')

        self.vecNames = np.zeros(len(vecs), 'i')
        for k, vec in enumerate(vecs):
            if vec in keyVec:
                self.vecNames[k] = keyVec[vec]
            else:
                self.vecNames[k] = 0
        # Remove the possible 0
        self.vecNames = self.vecNames[self.vecNames!=0]
        self.nvecs = len(self.vecNames)

        self.radius = gr.radius[::-1]

        # Potential extrapolation
        if potExtra and nrout != 0:
            rcmb = gr.radius[0]
            brCMB = gr.Br[..., 0]
            pot = ExtraPot(rcmb, brCMB, gr.minc, ratio_out=ratio_out,
                           nrout=nrout, cutCMB=True, deminc=deminc)
            self.radius = np.concatenate((gr.radius[::-1], pot.rout))
        #----------------------------------------
        # Scalars
        #----------------------------------------
        if deminc:
            if potExtra:
                self.scals = np.zeros((self.nscals, gr.nphi, gr.ntheta, gr.nr+nrout-1), 'f')
            else:
                self.scals = np.zeros((self.nscals, gr.nphi, gr.ntheta, gr.nr), 'f')
        else:
            if potExtra:
                self.scals = np.zeros((self.nscals, gr.npI, gr.ntheta, gr.nr+nrout-1), 'f')
            else:
                self.scals = np.zeros((self.nscals, gr.npI, gr.ntheta, gr.nr), 'f')

        for k, index in enumerate(self.scalNames):
            if index == -1: # radius
                self.scals[k, ...] = self.radius
            if index == 1: # entropy
                if deminc:
                    self.scals[k, :, :, 0:gr.nr] = symmetrize(gr.entropy[..., ::-1],
                                                              gr.minc)
                else:
                    self.scals[k, :, :, 0:gr.nr] = gr.entropy[..., ::-1]
            elif index == 2: # magnetic energy
                if deminc:
                    self.scals[k, :, :, 0:gr.nr] = symmetrize(0.5*(gr.Br[..., ::-1]**2+\
                                                        gr.Btheta[..., ::-1]**2+\
                                                        gr.Bphi[..., ::-1]**2), gr.minc)
                else:
                    self.scals[k, :, :, 0:gr.nr] = 0.5*(gr.Br[..., ::-1]**2+\
                                                        gr.Btheta[..., ::-1]**2+\
                                                        gr.Bphi[..., ::-1]**2)
            elif index == 3 or index == 7: # z-vorticity
                th3D = np.zeros_like(gr.vphi)
                rr3D = np.zeros_like(th3D)
                for i in range(gr.ntheta):
                    th3D[:, i, :] = gr.colatitude[i]
                for i in range(gr.nr):
                    rr3D[:, :, i] = gr.radius[i]
                s3D = rr3D * np.sin(th3D)
                dtheta = thetaderavg(gr.vphi*s3D)
                dr = rderavg(gr.vphi*s3D, eta=gr.radratio, spectral=True,
                             exclude=False)
                vs = gr.vr * np.sin(th3D) + gr.vtheta * np.cos(th3D) # 'vs'
                ds = np.sin(th3D)*dr + np.cos(th3D)/rr3D*dtheta
                vortz = -1./s3D*phideravg(vs, minc=gr.minc)+ds/s3D
                if deminc:
                    self.scals[k, :, :, 0:gr.nr] = symmetrize(vortz[..., ::-1],
                                                              gr.minc)
                else:
                    self.scals[k, :, :, 0:gr.nr] = vortz[..., ::-1]
                if index == 7: # Fluctuation
                    self.scals[k, ...] = self.scals[k, ...]-\
                                         self.scals[k, ...].mean(axis=0)
            elif index == 4: # radial velocity
                if deminc:
                    self.scals[k, :, :, 0:gr.nr] = symmetrize(gr.vr[..., ::-1],
                                                              gr.minc)
                else:
                    self.scals[k, :, :, 0:gr.nr] = gr.vr[..., ::-1]
            elif index == 5: # zonal velocity
                if deminc:
                    self.scals[k, :, :, 0:gr.nr] = symmetrize(gr.vphi[..., ::-1],
                                                              gr.minc)
                else:
                    self.scals[k, :, :, 0:gr.nr] = gr.vphi[..., ::-1]
            elif index == 6: # fluctuation of entropy/temperature
                if deminc:
                    self.scals[k, :, :, 0:gr.nr] = symmetrize(gr.entropy[..., ::-1],
                                                              gr.minc)
                else:
                    self.scals[k, :, :, 0:gr.nr] = gr.entropy[..., ::-1]
                self.scals[k, ...] = self.scals[k, ...]-\
                                     self.scals[k, ...].mean(axis=0)
            elif index == 8: # r-vorticity
                th3D = np.zeros_like(gr.vphi)
                rr3D = np.zeros_like(th3D)
                for i in range(gr.ntheta):
                    th3D[:, i, :] = gr.colatitude[i]
                for i in range(gr.nr):
                    rr3D[:, :, i] = gr.radius[i]
                s3D = rr3D * np.sin(th3D)
                vortr = 1./s3D*(thetaderavg(np.sin(th3D)*gr.vphi)-\
                        phideravg(gr.vtheta, gr.minc))
                if deminc:
                    self.scals[k, :, :, 0:gr.nr] = symmetrize(vortr[..., ::-1],
                                                              gr.minc)
                else:
                    self.scals[k, :, :, 0:gr.nr] = vortr[..., ::-1]
            elif index == 9: # fluctuation of entropy/temperature (substract 1-D time-average)
                bl = BLayers(iplot=False)
                if deminc:
                    self.scals[k, :, :, 0:gr.nr] = symmetrize(gr.entropy[..., ::-1],
                                                              gr.minc)
                else:
                    self.scals[k, :, :, 0:gr.nr] = gr.entropy[..., ::-1]
                self.scals[k, ...] = self.scals[k, ...]-bl.ss[::-1]
            elif index == 10: # kinetic energy
                if deminc:
                    self.scals[k, :, :, 0:gr.nr] = symmetrize(0.5*(gr.vr[..., ::-1]**2+\
                                                        gr.vtheta[..., ::-1]**2+\
                                                        gr.vphi[..., ::-1]**2), gr.minc)
                else:
                    self.scals[k, :, :, 0:gr.nr] = 0.5*(gr.vr[..., ::-1]**2+\
                                                        gr.vtheta[..., ::-1]**2+\
                                                        gr.vphi[..., ::-1]**2)

            elif index == 11: # radial mag field
                if deminc:
                    self.scals[k, :, :, 0:gr.nr] = symmetrize(gr.Br[..., ::-1],
                                                              gr.minc)
                else:
                    self.scals[k, :, :, 0:gr.nr] = gr.Br[..., ::-1]
            elif index == 12: #cylindrical velocity
                th3D = np.zeros_like(gr.vphi)
                rr3D = np.zeros_like(th3D)
                for i in range(gr.ntheta):
                    th3D[:, i, :] = gr.colatitude[i]
                vs = gr.vr * np.sin(th3D) + gr.vtheta * np.cos(th3D)

                if deminc:
                    self.scals[k, :, :, 0:gr.nr] = symmetrize(vs[..., ::-1], gr.minc)
                else:
                    self.scals[k, :, :, 0:gr.nr] = vs[..., ::-1]

            elif index == 13: # Colatitude
                th3D = np.zeros_like(gr.vphi)
                for i in range(gr.ntheta):
                    th3D[:, i, :] = gr.colatitude[i]

                if deminc:
                    self.scals[k, :, :, 0:gr.nr] = symmetrize(th3D[..., ::-1], gr.minc)
                else:
                    self.scals[k, :, :, 0:gr.nr] = th3D[..., ::-1]

            if index == 14: # chemical composition
                if deminc:
                    self.scals[k, :, :, 0:gr.nr] = symmetrize(gr.xi[..., ::-1],
                                                              gr.minc)
                else:
                    self.scals[k, :, :, 0:gr.nr] = gr.xi[..., ::-1]

            elif index == 15: # fluctuation of chemical composition
                if deminc:
                    self.scals[k, :, :, 0:gr.nr] = symmetrize(gr.xi[..., ::-1],
                                                              gr.minc)
                else:
                    self.scals[k, :, :, 0:gr.nr] = gr.xi[..., ::-1]
                self.scals[k, ...] = self.scals[k, ...]-\
                                     self.scals[k, ...].mean(axis=0)

            if potExtra:
                if index == 2:
                    self.scals[k, :, :, gr.nr:] = 0.5*(pot.brout**2+pot.btout**2+\
                                                       pot.bpout**2)
                if index == 11:
                    self.scals[k, :, :, gr.nr:] = pot.brout


        #----------------------------------------
        # Vectors
        #----------------------------------------
        if deminc:
            if potExtra:
                self.vecr = np.zeros((self.nvecs, gr.nphi, gr.ntheta, gr.nr+nrout-1), 'f')
                self.vect = np.zeros_like(self.vecr)
                self.vecp = np.zeros_like(self.vecr)
            else:
                self.vecr = np.zeros((self.nvecs, gr.nphi, gr.ntheta, gr.nr), 'f')
                self.vect = np.zeros_like(self.vecr)
                self.vecp = np.zeros_like(self.vecr)
        else:
            if potExtra:
                self.vecr = np.zeros((self.nvecs, gr.npI, gr.ntheta, gr.nr+nrout-1), 'f')
                self.vect = np.zeros_like(self.vecr)
                self.vecp = np.zeros_like(self.vecr)
            else:
                self.vecr = np.zeros((self.nvecs, gr.npI, gr.ntheta, gr.nr), 'f')
                self.vect = np.zeros_like(self.vecr)
                self.vecp = np.zeros_like(self.vecr)


        for k, index in enumerate(self.vecNames):
            if index == 1:
                if deminc:
                    self.vecr[k, :, :, 0:gr.nr] = symmetrize(gr.vr[..., ::-1],
                                                             gr.minc)
                    self.vect[k, :, :, 0:gr.nr] = symmetrize(gr.vtheta[..., ::-1],
                                                             gr.minc)
                    self.vecp[k, :, :, 0:gr.nr] = symmetrize(gr.vphi[..., ::-1],
                                                             gr.minc)
                else:
                    self.vecr[k, :, :, 0:gr.nr] = gr.vr[..., ::-1]
                    self.vect[k, :, :, 0:gr.nr] = gr.vtheta[..., ::-1]
                    self.vecp[k, :, :, 0:gr.nr] = gr.vphi[..., ::-1]
            elif index == 2:
                if deminc:
                    self.vecr[k, :, :, 0:gr.nr] = symmetrize(gr.vr[..., ::-1],
                                                             gr.minc)
                    self.vect[k, :, :, 0:gr.nr] = symmetrize(gr.vtheta[..., ::-1],
                                                             gr.minc)
                    self.vecp[k, :, :, 0:gr.nr] = symmetrize(gr.vphi[..., ::-1],
                                                             gr.minc)
                else:
                    self.vecr[k, :, :, 0:gr.nr] = gr.vr[..., ::-1]
                    self.vect[k, :, :, 0:gr.nr] = gr.vtheta[..., ::-1]
                    self.vecp[k, :, :, 0:gr.nr] = gr.vphi[..., ::-1]
                self.vecr[k, ...] = self.vecr[k, ...]-self.vecr[k, ...].mean(axis=0)
                self.vect[k, ...] = self.vect[k, ...]-self.vect[k, ...].mean(axis=0)
                self.vecp[k, ...] = self.vecp[k, ...]-self.vecp[k, ...].mean(axis=0)
            elif index == 3:
                if deminc:
                    self.vecr[k, :, :, 0:gr.nr] = symmetrize(gr.Br[..., ::-1],
                                                             gr.minc)
                    self.vect[k, :, :, 0:gr.nr] = symmetrize(gr.Btheta[..., ::-1],
                                                             gr.minc)
                    self.vecp[k, :, :, 0:gr.nr] = symmetrize(gr.Bphi[..., ::-1],
                                                             gr.minc)
                else:
                    self.vecr[k, :, :, 0:gr.nr] = gr.Br[..., ::-1]
                    self.vect[k, :, :, 0:gr.nr] = gr.Btheta[..., ::-1]
                    self.vecp[k, :, :, 0:gr.nr] = gr.Bphi[..., ::-1]
            elif index == 4:
                if deminc:
                    self.vecr[k, :, :, 0:gr.nr] = symmetrize(gr.Br[..., ::-1],
                                                             gr.minc)
                    self.vect[k, :, :, 0:gr.nr] = symmetrize(gr.Btheta[..., ::-1],
                                                             gr.minc)
                    self.vecp[k, :, :, 0:gr.nr] = symmetrize(gr.Bphi[..., ::-1],
                                                             gr.minc)
                else:
                    self.vecr[k, :, :, 0:gr.nr] = gr.Br[..., ::-1]
                    self.vect[k, :, :, 0:gr.nr] = gr.Btheta[..., ::-1]
                    self.vecp[k, :, :, 0:gr.nr] = gr.Bphi[..., ::-1]
                self.vecr[k, ...] = self.vecr[k, ...]-self.vecr[k, ...].mean(axis=0)
                self.vect[k, ...] = self.vect[k, ...]-self.vect[k, ...].mean(axis=0)
                self.vecp[k, ...] = self.vecp[k, ...]-self.vecp[k, ...].mean(axis=0)

            if potExtra:
                if index == 3:
                    self.vecr[k, :, :, gr.nr:] = pot.brout
                    self.vect[k, :, :, gr.nr:] = pot.btout
                    self.vecp[k, :, :, gr.nr:] = pot.bpout
                elif index == 4:
                    pot.brout = pot.brout-pot.brout.mean(axis=0)
                    pot.btout = pot.btout-pot.btout.mean(axis=0)
                    pot.bpout = pot.bpout-pot.bpout.mean(axis=0)
                    self.vecr[k, :, :, gr.nr:] = pot.brout
                    self.vect[k, :, :, gr.nr:] = pot.btout
                    self.vecp[k, :, :, gr.nr:] = pot.bpout

        if outType == 'vts':
            self.writeVTS(filename, nFiles)

        elif outType == 'vti':
            self.writeVTI(filename, nx, ny, nz)

    def writeVTS(self, filename, nFiles):
        """
        This function stores the output on a structured-grid vts file.

        :param filename: the file name of the output (without extension)
        :type filename: str
        :param nFiles: number of outpute files (in case of pvts)
        :type nFiles: int
        """
        if self.nvecs == 0:
            if nFiles == 1 and np.size(self.scals)+3*np.size(self.vecr) < 512**3:
                vts_scal(filename, self.radius, self.scals, self.scalNames,
                         self.minc)
            elif nFiles > 1:
                pvts_scal(filename, self.radius, self.scals, self.scalNames,
                          nFiles, self.minc)

            else:
                pvts_scal(filename, self.radius, self.scals, self.scalNames,
                          8, self.minc)
        else:
            if nFiles == 1 and np.size(self.scals)+3*np.size(self.vecr) < 512**3:
                vts(filename, self.radius, self.vecr, self.vect, self.vecp,
                    self.scals, self.scalNames, self.vecNames, self.minc)
            elif nFiles > 1:
                pvts(filename, self.radius, self.vecr, self.vect, self.vecp,
                     self.scals, self.scalNames, self.vecNames, nFiles, self.minc)
            else:
                pvts(filename, self.radius, self.vecr, self.vect, self.vecp,
                     self.scals, self.scalNames, self.vecNames, 8, self.minc)

    def writeVTI(self, filename, nx=96, ny=96, nz=96):
        """
        In this case, the output is extrapolated on a cartesian grid
        and then written in a vti file.

        :param filename: the file name of the output (without extension)
        :type filename: str
        :param nx: number of grid points in the x direction
        :type nx: int
        :param ny: number of grid points in the x direction
        :type ny: int
        :param nz: number of grid points in the x direction
        :type nz: int
        """
        scals_cart, gridMax, spacing = sph2cart_scal(self.scals, \
                                                     self.radius, nx, ny, nz,
                                                     self.minc)
        if self.nvecs == 0:
            vti_scal(filename, scals_cart, self.scalNames, self.minc, \
                     gridMax, spacing)
        else:
            vecx, vecy, vecz = sph2cart_vec(self.vecr, self.vect, self.vecp, \
                                        self.radius, nx, ny, nz, self.minc)
            vti(filename, vecx, vecy, vecz, scals_cart, self.scalNames,
                self.vecNames, self.minc, gridMax, spacing)


if __name__ == '__main__':
    ivar = 1
    g = MagicGraph(ivar=ivar)
    Graph2Vtk(g, scals=['tfluct', 'emag'], vecs=['b'], potExtra=True)
