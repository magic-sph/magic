# -*- coding: utf-8 -*-
from scipy.ndimage import map_coordinates
from scipy.interpolate import interp1d
from magic import MagicGraph, BLayers
from magic.setup import labTex
from magic.libmagic import symmetrize, thetaderavg, rderavg, phideravg
import matplotlib.pyplot as P
import numpy as N
import sys
if sys.version_info.major == 3:
    from potential3 import *
    from vtklib3 import *
elif  sys.version_info.major == 2:
    from potential2 import *
    from vtklib2 import *

__author__  = "$Author$"
__date__   = "$Date$"
__version__ = "$Revision$"


def sph2cart(br, bt, bp, radius, nx=96, ny=96, nz=96, hydro=False):
    """
    A routine that transforms from spherical to cartesian
    coordinates.

    :param br,bt,bp: the three components of a vector field
                     in spherical coodrinates
    :param radius: the input radius
    :param nx,ny,nz: the resolution of the cartesian grid
    """
    np, nt, nr = br.shape
    phi = N.linspace(-N.pi, N.pi, np)
    theta = N.linspace(0., N.pi, nt)
    # Cube: take care of the sqrt(3.) !!!
    if hydro:
        gridMax = radius.max()/N.sqrt(3.)
    else:
        gridMax = radius.max()/N.sqrt(3.)
    spacing = 2.*gridMax/(nx-1)
    Z,Y,X = N.mgrid[-1:1:nz*1j,-1:1:ny*1j,-1:1:nx*1j]*gridMax
    new_r = N.sqrt(X**2+Y**2+Z**2).ravel()
    new_phi = N.arctan2(Y, X).ravel()
    new_theta = N.arctan2(N.sqrt(X**2+Y**2), Z).ravel()
    del X,Y,Z

    ir = interp1d(radius, N.arange(len(radius)), bounds_error=False)
    it = interp1d(theta, N.arange(len(theta)), bounds_error=False)
    ip = interp1d(phi, N.arange(len(phi)), bounds_error=False)

    new_ir = ir(new_r)
    new_it = it(new_theta)
    new_ip = ip(new_phi)

    new_ir[new_r < radius.min()] = 0.
    new_it[new_r < radius.min()] = 0.
    new_it[new_r < radius.min()] = 0.

    coords = N.array([new_ip, new_it, new_ir])

    br_cart = map_coordinates(br, coords)
    bt_cart = map_coordinates(bt, coords)
    bp_cart = map_coordinates(bp, coords)

    br_cart[new_r < radius.min()] = 0.
    bt_cart[new_r < radius.min()] = 0.
    bp_cart[new_r < radius.min()] = 0.

    if hydro:
        br_cart[new_r > radius.max()] = 0.
        bt_cart[new_r > radius.max()] = 0.
        bp_cart[new_r > radius.max()] = 0.

    del coords

    bs = br_cart*N.sin(new_theta)+bt_cart*N.cos(new_theta)
    bx = bs*N.cos(new_phi) - bp_cart*N.sin(new_phi)
    by = bs*N.sin(new_phi) + bp_cart*N.cos(new_phi)
    bz = br_cart*N.cos(new_theta) - bt_cart*N.sin(new_theta)
    del new_theta, new_phi

    return bx,by,bz,br_cart,bp_cart,new_r,gridMax,spacing


def proj(theta, phi):
    x = 2.*N.sqrt(2.)*N.sin(theta)*N.sin(phi/2.)/\
        N.sqrt(1.+N.sin(theta)*N.cos(phi/2.))
    y = N.sqrt(2.)*N.cos(theta)/N.sqrt(1.+N.sin(theta)*N.cos(phi/2.))
    return x, y

class ExtraPot:
    """
    This class is used to compute the potential field extrapolation,
    starting from the surface radial field. Thanks to two different methods,
    the output can be written in VTK file compatible with paraview or
    mayavi.
    """

    def __init__(self, rcmb, brcmb, minc, ratio_out=2., nrout=32, cutCMB=False):
        """
        :param bcmb: the surface radial field, array of dimension (np, nt)
        :param rcmb: the value of the radius at the surface
        :param minc: symmetry
        :param ratio_out: the ratio of maximum radius to surface radius
        :param nrout: the number of radial point (linearly spaced) of the
                      extrapolated field
        :param cutCMB: a boolean if one wants to remove the first grid point
                       (useful if one wants to merge the graphic file with
                       the extrapolation)
        """
        self.rcmb = rcmb
        self.brcmb = brcmb
        self.minc = minc
        self.nrout = nrout
        self.np, self.nt = self.brcmb.shape
        brcmbfour = N.fft.fft(self.brcmb, axis=0)/(4.*N.pi*self.np)

        self.rout = N.linspace(self.rcmb, ratio_out*rcmb, self.nrout)
        if cutCMB:
            self.rout = self.rout[1:]
            self.nrout = self.nrout -1

        self.brout = N.zeros((self.np, self.nt, self.nrout), 'f')
        self.btout = N.zeros_like(self.brout)
        self.bpout = N.zeros_like(self.brout)

        for k, radius  in enumerate(self.rout):
            radratio = self.rcmb/radius
            brm, btm, bpm =  extrapolate(brcmbfour, radratio, self.minc)
            brsurf = N.fft.ifft(brm, axis=0)*self.np
            self.brout[..., k] = brsurf.real
            btsurf = N.fft.ifft(btm, axis=0)*self.np
            self.btout[..., k] = btsurf.real
            bpsurf = N.fft.ifft(bpm, axis=0)*self.np
            self.bpout[..., k] = bpsurf.real

        self.brout = symmetrize(self.brout, self.minc)
        self.btout = symmetrize(self.btout, self.minc)
        self.bpout = symmetrize(self.bpout, self.minc)

    def avg(self, field='br', levels=12, cm='RdYlBu_r', normed=True,
            vmax=None, vmin=None):
        """
        A small routine to plot azimuthal averages of the extrapolated
        field.

        :param field: the field you want to plot
        :param levels: the number of levels of the colormap
        :param cm: the name of the colormap
        """

        if field == 'br':
            if labTex:
                label = r'$B_r$'
            else:
                label = 'Br'
            data = self.brout
        elif field == 'bp':
            if labTex:
                label = r'$B_\phi$'
            else:
                label = 'Bphi'
            data = self.bpout

        phiavg = data.mean(axis=0)
        th = N.linspace(0, N.pi, self.nt)
        rr, tth = N.meshgrid(self.rout, th)
        xx = rr * N.sin(tth)
        yy = rr * N.cos(tth)

        fig = P.figure(figsize=(4,8))
        fig.subplots_adjust(top=0.99, bottom=0.01, right=0.99, left=0.01)
        ax = fig.add_subplot(111, frameon=False)
        cmap = P.get_cmap(cm)
        if vmax is not None and vmin is not None:
            normed = False
            cs = N.linspace(vmin, vmax, levels)
            im = ax.contourf(xx, yy, phiavg, cs, cmap=cmap, extend='both')
        else:
            cs = levels
            im = ax.contourf(xx, yy, phiavg, cs, cmap=cmap)
        ax.plot(self.rcmb*N.sin(th), self.rcmb*N.cos(th), 'k-')
        ax.axis('off')

    def writeVTS(self, filename):
        """
        In this case, the output is directly written on the spherical
        grid, i.e. a vts file.

        :param filename: the file name of the output (without extension)
        """
        vts(filename, self.rout,self.brout,self.btout,self.bpout)

    def writeVTI(self, filename, nx=96, ny=96, nz=96, nscals=3, nvecs=1):
        """
        In this case, the output is extrapolated on a cartesian grid
        and then written in a vti file.

        nx,ny,nz are used to specify the resolution of the cartesian grid.
        """
        bx, by, bz, br_cart, bp_cart, new_r, gridMax, spacing = \
          sph2cart(self.brout, self.btout, self.bpout, self.rout, nx, ny, nz)
        field = N.vstack((bx, by, bz))
        scals = N.zeros((nscals, by.shape[0]), 'f')
        vecs = N.zeros((nvecs, 3, by.shape[0]), 'f')
        scals[0, :] = new_r
        scals[1, :] = br_cart
        scals[2, :] = bx**2+by**2+bz**2
        vecs[0, ...] = field
        writevtr(filename, scals, vecs, gridMax, spacing, nx, ny, nz)


class TotalField:

    def __init__(self, g, ratio_out=2, nrout=0, hydro=False, vortOut=False,
                 ttOut=False, fluct=False, deminc=True):
        """
        For hydro run, specify  hydro=True

        :param ratio_out: for potential extrapolation (rout/rcmb)
        :param nrout: number of radial grid points in the outer spherical shell
                      (potential extrapolation)
        :param hydro: a logical to indicate if this is a non-magnetic setup
        :param vortOut: a logical to compute z-vorticity instead of vr
        :param ttOut: a logical to compute entropy instead of vr
        :param fluct: a logical to substract the axisymmetric part
        :param deminc: a logical to indicate if one wants do "de-minc"
        """
        self.hydro = hydro
        self.vortOut = vortOut
        if not self.hydro:
            rcmb = g.radius[0]
            brCMB = g.Br[..., 0]
            if nrout != 0:
                pot = ExtraPot(rcmb, brCMB, g.minc, ratio_out=ratio_out, 
                               nrout=nrout, cutCMB=True)

            if deminc:
                self.br = symmetrize(g.Br[..., ::-1], g.minc)
                self.bt = symmetrize(g.Btheta[..., ::-1], g.minc)
                self.bp = symmetrize(g.Bphi[..., ::-1], g.minc)
            else:
                self.br = g.Br[..., ::-1]
                self.bt = g.Btheta[..., ::-1]
                self.bp = g.Bphi[..., ::-1]

            if fluct:
                self.br = self.br-self.br.mean(axis=0)
                self.bt = self.bt-self.bt.mean(axis=0)
                self.bp = self.bp-self.bp.mean(axis=0)

            if nrout != 0:
                self.radius = N.concatenate((g.radius[::-1], pot.rout))
                self.br = N.concatenate((self.br, pot.brout), axis=2)
                self.bt = N.concatenate((self.bt, pot.btout), axis=2)
                self.bp = N.concatenate((self.bp, pot.bpout), axis=2)
            else:
                self.radius = g.radius[::-1]
        else:
            self.radius = g.radius[::-1]


            if deminc:
                self.br = symmetrize(g.vr[..., ::-1], g.minc)
                self.bt = symmetrize(g.vtheta[..., ::-1], g.minc)
                self.bp = symmetrize(g.vphi[..., ::-1], g.minc)
            else:
                self.br = g.vr[..., ::-1]
                self.bt = g.vtheta[..., ::-1]
                self.bp = g.vphi[..., ::-1]

            if fluct:
                self.br = self.br-self.br.mean(axis=0)
                self.bt = self.bt-self.bt.mean(axis=0)
                self.bp = self.bp-self.bp.mean(axis=0)

            if ttOut:
                if deminc:
                    self.entropy = symmetrize(g.entropy[..., ::-1], g.minc)
                else:
                    self.entropy = g.entropy[..., ::-1]
                if fluct:
                    self.entropy = self.entropy-self.entropy.mean(axis=0)
                #bl = BLayers(iplot=False)
                #self.entropy = self.entropy-bl.ss[::-1]

            if self.vortOut:
                th3D = N.zeros_like(g.vphi)
                rr3D = N.zeros_like(th3D)
                for i in range(g.ntheta):
                    th3D[:, i, :] = g.colatitude[i]
                for i in range(g.nr):
                    rr3D[:, :, i] = g.radius[i]
                s3D = rr3D * N.sin(th3D)
                if fluct:
                    dtheta = thetaderavg((g.vphi-g.vphi.mean(axis=0))*s3D)
                    dr = rderavg((g.vphi-g.vphi.mean(axis=0))*s3D,
                                 eta=g.radratio, spectral=True, exclude=False)
                    vs = (g.vr-g.vr.mean(axis=0))*N.sin(th3D) + \
                         (g.vtheta-g.vtheta.mean(axis=0))*N.cos(th3D) # 'vs'
                else:
                    dtheta = thetaderavg(g.vphi*s3D)
                    dr = rderavg(g.vphi*s3D, eta=g.radratio, spectral=True,
                                 exclude=False)
                    vs = g.vr * N.sin(th3D) + g.vtheta * N.cos(th3D) # 'vs'

                ds = N.sin(th3D)*dr + N.cos(th3D)/rr3D*dtheta
                self.vortz = -1./s3D*phideravg(vs)+ds/s3D

                self.vortr = 1./s3D*(thetaderavg(N.sin(th3D)*g.vphi)-phideravg(g.vtheta))

                del dr, dtheta, ds, rr3D, th3D, s3D

                if deminc:
                    self.vortz = symmetrize(self.vortz[..., ::-1], g.minc)
                    self.vortr = symmetrize(self.vortr[..., ::-1], g.minc)
                else:
                    self.vortz = self.vortz[..., ::-1]
                    self.vortr = self.vortr[..., ::-1]


    def writeVTI(self, filename, nx=96, ny=96, nz=96, nscals=3, nvecs=1):
        """
        In this case, the output is extrapolated on a cartesian grid
        and then written in a vti file.

        nx,ny,nz are used to specify the resolution of the cartesian grid.
        """
        filename += '_tot'
        if not self.vortOut:
            bx, by, bz, br_cart, bp_cart, new_r, gridMax, spacing = \
                    sph2cart(self.br, self.bt, self.bp, self.radius, 
                             nx, ny, nz, self.hydro)
        else:
            bx, by, bz, br_cart, bp_cart, new_r, gridMax, spacing = \
                    sph2cart(self.vortz, self.bt, self.bp, self.radius, 
                             nx, ny, nz, self.hydro)
        if not self.hydro:
            field = N.vstack((bx, by, bz))
            scals = N.zeros((nscals, by.shape[0]), 'f')
            vecs = N.zeros((nvecs, 3, by.shape[0]), 'f')
            scals[0, :] = new_r
            scals[1, :] = br_cart
            scals[2, :] = bx**2+by**2+bz**2
            vecs[0, ...] = field

        else:
            field = N.vstack((bx, by, bz))
            scals = N.zeros((nscals, by.shape[0]), 'f')
            vecs = N.zeros((1, 3, by.shape[0]), 'f')
            scals[0, :] = new_r
            scals[1, :] = br_cart
            scals[2, :] = bp_cart
            vecs[0, ...] = field
        writevtr(filename, scals, vecs, gridMax, spacing, nx, ny, nz)

    def writeVTS(self, filename, minc=1):
        """
        In this case, the output is directly written on the spherical
        grid, i.e. a vts file.

        :param filename: the file name of the output (without extension)
        :param minc: azymuthal symmetry (default 1)
        """
        filename += '_tot'
        if self.vortOut:
            #vts(filename, self.radius, self.vortz, self.bt, self.bp, minc)
            vts(filename, self.radius, self.br, self.bt, self.bp, abs(self.br),
                self.entropy, minc)
            #pvts(filename, self.radius, self.vortz, 16, minc)
        else:
            #vts(filename, self.radius, self.br, self.bt, self.bp, minc)
            pvts(filename, self.radius, self.entropy, 8, minc)


if __name__ == '__main__':
    ivar = 1
    g = MagicGraph(ivar=ivar)
    filename = 'G_%i_%s' % (ivar, g.tag)

    t = TotalField(g, nrout=64, ratio_out=2.)
    t.writeVTI(filename)
    t.writeVTS(filename)

    """
    rcmb = g.radius[0]
    brCMB = g.Br[..., 0]
    pot = ExtraPot(rcmb, brCMB, g.minc,ratio_out=2, nrout=64)

    pot.writeVTS(filename)
    pot.writeVTI(filename)
    """
