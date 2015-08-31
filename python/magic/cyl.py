# -*- coding: utf-8 -*-
import numpy as N
import matplotlib.pyplot as P
from .libmagic import anelprof, cylSder, cylZder, phideravg, symmetrize, cut
from magic import MagicGraph, MagicSetup
from magic.setup import labTex
from scipy.ndimage import map_coordinates
from scipy.interpolate import interp1d
import os, pickle

__author__  = "$Author$"
__date__   = "$Date$"
__version__ = "$Revision$"



def sph2cyl_plane(data, rad, ns, nz):
    """
    Subroutine that extrapolates a phi-slice of a spherical shell on 
    a cylindrical grid

    :param data: a list of 2-D arrays [(ntheta, nr), (ntheta, nr), ...]
    :param rad: radius
    :param ns: number of grid points in s direction
    :param nz: number of grid points in z direction
    """
    ntheta, nr = data[0].shape

    radius = rad[::-1]

    theta = N.linspace(0., N.pi, ntheta)

    Z, S = N.mgrid[-radius.max():radius.max():nz*1j,0:radius.max():ns*1j]

    new_r = N.sqrt(S**2+Z**2).ravel()
    new_theta = N.arctan2(S, Z).ravel()
    ir = interp1d(radius, N.arange(len(radius)), bounds_error=False)
    it = interp1d(theta, N.arange(len(theta)), bounds_error=False)

    new_ir = ir(new_r)
    new_it = it(new_theta)
    new_ir[new_r > radius.max()] = len(radius)-1.
    new_ir[new_r < radius.min()] = 0.

    coords = N.array([new_it, new_ir])

    output = []
    for dat in data:
        dat_cyl = map_coordinates(dat[:, ::-1], coords, order=3)
        dat_cyl[new_r > radius.max()] = 0.
        dat_cyl[new_r < radius.min()] = 0.
        dat_cyl = dat_cyl.reshape((nz, ns))
        output.append(dat_cyl)

    return Z, S, output


def zavg(input, radius, ns, minc, save=True, filename='vp.pickle', normed=True):
    """
    A routine that computes a z-integration of a list of input arrays (on the spherical
    grid). List works well for 2-D (phi-slice) arrays, only one single in 3-D (too
    demanding otherwise).

    :param input: a list of 2-D or 3-D arrays 
    :param radius: spherical radius
    :param ns: radial resolution of the cylindrical grid (nz=2*ns)
    :param minc: symmetry
    :param save: a boolean to specify if one wants to save the outputs into 
                 a pickle (default is True)
    :param filename: name of the output pickle
    :param normed: a boolean to specify if ones wants to simply integrate over
                   z or compute an average (default is True: average)
    """
    nz = 2*ns
    ro = radius[0]
    ri = radius[-1]
    z = N.linspace(-ro, ro, nz)
    cylRad = N.linspace(0., ro, ns)
    cylRad = cylRad[1:-1]
    
    height = N.zeros_like(cylRad)
    height[cylRad>=ri] = 2.*N.sqrt(ro**2-cylRad[cylRad>=ri]**2)
    height[cylRad<ri] = 2.*(N.sqrt(ro**2-cylRad[cylRad<ri]**2)\
                                 -N.sqrt(ri**2-cylRad[cylRad<ri]**2))

    if len(input[0].shape) == 3:
        nphi = input[0].shape[0]
        phi = N.linspace(0., 2.*N.pi/minc, nphi)
        output = N.zeros((nphi, ns-2), 'f')
        for iphi in range(nphi):
            print(iphi)
            Z, S, out2D = sph2cyl_plane([input[0][iphi, ...]], radius, ns, nz)
            S = S[:, 1:-1]
            Z = Z[:, 1:-1]
            output[iphi, :] = N.trapz(out2D[0][:, 1:-1], z, axis=0)
            if normed:
                output[iphi, :] /= height

        if save:
            nphi, ntheta, nr = input[0].shape
            file = open(filename, 'wb')
            pickle.dump([cylRad, phi, output], file) # cylindrical average
            pickle.dump([radius, phi, input[0][:, ntheta/2, :]], file) # equatorial cut
            file.close()
        return height, cylRad, phi, output
    elif len(input[0].shape) == 2:
        Z, S, out2D = sph2cyl_plane(input, radius, ns, nz)
        S = S[:, 1:-1]
        Z = Z[:, 1:-1]
        output = []
        outIntZ = N.zeros((ns-2), 'f')
        for k,out in enumerate(out2D):
            outIntZ = N.trapz(out[:, 1:-1], z, axis=0)
            if normed:
                outIntZ /= height
            output.append(outIntZ)

        if save:
            file = open(filename, 'wb')
            pickle.dump([radius,  cylRad, height], file) # cylindrical average
            for k,out in enumerate(output):
                pickle.dump(out, file) # cylindrical average
                ntheta, nr = input[k].shape
                pickle.dump(input[k][ntheta/2, :], file) # equatorial cut
            file.close()

        return height, cylRad, output


def sph2cyl(g, ns=None, nz=None):
    if ns is None or nz is None:
        ns = g.nr ; nz = 2*ns

    theta = N.linspace(0., N.pi, g.ntheta)
    radius = g.radius[::-1]

    Z, S = N.mgrid[-radius.max():radius.max():nz*1j,0:radius.max():ns*1j]

    new_r = N.sqrt(S**2+Z**2).ravel()
    new_theta = N.arctan2(S, Z).ravel()
    ir = interp1d(radius, N.arange(len(radius)), bounds_error=False)
    it = interp1d(theta, N.arange(len(theta)), bounds_error=False)

    new_ir = ir(new_r)
    new_it = it(new_theta)
    new_ir[new_r > radius.max()] = len(radius)-1.
    new_ir[new_r < radius.min()] = 0.

    coords = N.array([new_it, new_ir])

    vr_cyl = N.zeros((g.npI, nz, ns), dtype='f')
    vp_cyl = N.zeros_like(vr_cyl)
    vt_cyl = N.zeros_like(vr_cyl)
    for k in range(g.npI):
        print(k)
        dat = map_coordinates(g.vphi[k, :, ::-1], coords, order=3)
        dat[new_r > radius.max()] = 0.
        dat[new_r < radius.min()] = 0.
        vp_cyl[k, ...] = dat.reshape((nz, ns))
        dat = map_coordinates(g.vtheta[k, :, ::-1], coords, order=3)
        dat[new_r > radius.max()] = 0.
        dat[new_r < radius.min()] = 0.
        vt_cyl[k, ...] = dat.reshape((nz, ns))
        dat = map_coordinates(g.vr[k, :, ::-1], coords, order=3)
        dat[new_r > radius.max()] = 0.
        dat[new_r < radius.min()] = 0.
        vr_cyl[k, ...] = dat.reshape((nz, ns))

    th3D = N.zeros((g.npI, nz, ns), dtype='f')
    for i in range(g.npI):
        th3D[i, ...] = N.arctan2(S, Z)
    vs = vr_cyl * N.sin(th3D) + vt_cyl * N.cos(th3D)
    vz = vr_cyl * N.cos(th3D) - vt_cyl * N.sin(th3D)

    return S, Z, vs, vp_cyl, vz


class Cyl(MagicSetup):

    def __init__(self, lastvar=1, avg=False, datadir='.', ns=None):
        MagicSetup.__init__(self, datadir)
         
        self.datadir = datadir

        filename = '%sG_%i.%s' % ('cyl', lastvar, self.tag)
        if not os.path.exists(filename):
            print("sph2cyl...")
            gr = MagicGraph(ivar=lastvar, datadir=self.datadir)
            if ns is None:
                self.ns = gr.nr
                self.nz = 2*self.ns
            else:
                self.ns = ns
                self.nz = 2*ns
            self.nphi = gr.nphi
            self.npI = gr.npI
            self.minc = gr.minc
            self.ro = gr.radius[0]
            self.ri = gr.radius[-1]
            self.S, self.Z, self.vs, self.vphi, self.vz = sph2cyl(gr, 
                                    self.ns, self.nz)
            file = open(filename, 'wb')
            pickle.dump([self.ns, self.nz, self.nphi, self.npI, self.minc], file)
            pickle.dump([self.ro, self.ri], file)
            pickle.dump([self.S, self.Z, self.vs, self.vphi, self.vz], 
                        file)
            file.close()
        else:
            print("read cyl file")
            file = open(filename, 'r')
            self.ns, self.nz, self.nphi, self.npI, self.minc = pickle.load(file)
            self.ro, self.ri = pickle.load(file)
            self.S, self.Z, self.vs, self.vphi, self.vz = \
                         pickle.load(file)
            file.close()
        self.radius = N.linspace(0., self.ro, self.ns)
        temp0, rho0, beta0 = anelprof(N.linspace(self.ro, self.ri, self.ns), 
                                     self.strat, self.polind)
        rho = N.zeros((self.nphi/2, self.ns), dtype='f')
        beta = N.zeros_like(rho)
        for i in range(self.nphi/2):
            rho[i, :] = rho0
            beta[i, :] = beta0
        Z, S, [rho, beta] = sph2cyl_plane([rho,beta], 
                                 N.linspace(self.ro, self.ri, self.ns), 
                                 self.ns, self.nz)
        self.rho = N.zeros_like(self.vs)
        self.beta = N.zeros_like(self.vs)
        for i in range(self.npI):
            self.rho[i, ...] = rho
            self.beta[i, ...] = beta
        self.z = N.linspace(-self.ro, self.ro, self.nz)

    def surf(self, field='Bphi', r=0.85, vmin=None, vmax=None, 
             levels=16, cm='RdYlBu_r', normed=True, figsize=None):
        r /= (1-self.ri/self.ro) # as we give a normalised radius
        ind = N.nonzero(N.where(abs(self.radius-r) \
                        == min(abs(self.radius-r)), 1, 0))
        indPlot = ind[0][0]

        if field in ('Vr', 'vr', 'Ur', 'ur'):
            data = self.vp
            label = 'Radial velocity'
        elif field in ('Vphi', 'vphi', 'Uphi', 'uphi', 'up', 'Up', 'Vp', 'vp'):
            data = self.vphi
            if labTex:
                label = r'$V_{\phi}$'
            else:
                label = 'vphi'
        elif field in ('Vs', 'vs'):
            data = self.vs
            label = 'Vs'
        elif field in ('Vz', 'vz'):
            data = self.vz
            label = 'Vz'

        phi = N.linspace(0., 2.*N.pi, self.nphi)

        data[..., indPlot] = cut(data[..., indPlot], vmax, vmin)
        data = symmetrize(data, self.minc)

        cmap = P.get_cmap(cm)

        fig = P.figure()
        ax = fig.add_subplot(111)
        im = ax.contourf(phi, self.z, data[..., indPlot].T, levels, cmap=cmap, 
                        aa=True)
        rad = self.radius[indPlot] * (1. - self.ri/self.ro)
        if labTex:
            ax.set_xlabel(r'$\phi$', fontsize=18)
            ax.set_ylabel(r'$z$', fontsize=18)
            ax.set_title('%s: $r/r_o$ = %.3f' % (label, rad), fontsize=24)
        else:
            ax.set_xlabel('phi', fontsize=18)
            ax.set_ylabel('z', fontsize=18)
            ax.set_title('%s: r/ro = %.3f' % (label, rad), fontsize=24)
        cbar = P.colorbar(im)

        if field not in ['entropy', 's', 'S'] and normed is True:
            im.set_clim(-max(abs(data[..., indPlot].max()), 
                             abs(data[..., indPlot].min())), 
                         max(abs(data[..., indPlot].max()),
                             abs(data[..., indPlot].min())))

    def equat(self, field='vs', levels=16, cm='RdYlBu_r', normed=True, vmax=None,
              vmin=None):
        """
        Plot the equatorial plane of a given field
        """
        if field in ('Vr', 'vr', 'Ur', 'ur'):
            data = self.vp
            label = 'Radial velocity'
        elif field in ('beta'):
            data = self.beta
            if labTex:
                label = r'$\beta$'
            else:
                label = 'beta'
        elif field in ('Vphi', 'vphi', 'Uphi', 'uphi', 'up', 'Up', 'Vp', 'vp'):
            data = self.vphi
            if labTex:
                label = r'$v_{\phi}$'
            else:
                label = 'vphi'
        elif field in ('Vs', 'vs'):
            data = self.vs
            label = r'$v_s$'
        elif field in ('Vz', 'vz'):
            data = self.vz
            if labTex:
                label = r'$v_z$'
            else:
                label = 'vz'
        elif field in ('dvz'):
            data =  cylZder(self.z, self.vz)
            if labTex:
                label = r'$\partial v_z/\partial z$'
            else:
                label = 'dvzdz'
        elif field in ('anel'):
            betas = cylSder(self.radius, N.log(self.rho))
            betaz = cylZder(self.z, N.log(self.rho))
            data = self.vs * betas + self.vz * betaz
            if labTex:
                label = r'$\beta v_r$'
            else:
                label = 'beta vr'
        elif field in ('Cr', 'cr'):
            vp = self.vphi.copy()-self.vphi.mean(axis=0) # convective vp
            data =  self.rho * self.vs * vp 
            if labTex:
                label = r'$\langle \rho v_s v_\phi\rangle$'
            else:
                label = 'rho vs vphi'

        equator = data[:, self.nz/2,:]
        equator = cut(equator, vmax, vmin)
        equator = symmetrize(equator, self.minc)

        phi = N.linspace(0., 2.*N.pi, self.nphi)
        rr, pphi = N.meshgrid(self.radius, phi)
        xx = rr * N.cos(pphi)
        yy = rr * N.sin(pphi)

        fig = P.figure(figsize=(8.25, 6))
        ax = fig.add_subplot(111, frameon=False)
        cmap = P.get_cmap(cm)
        im = ax.contourf(xx, yy, equator, levels, cmap=cmap)
        ax.plot(self.ri * N.cos(phi), self.ri*N.sin(phi), 'k-')
        ax.plot(self.ro * N.cos(phi), self.ro*N.sin(phi), 'k-')
        ax.set_title(label, fontsize=24)
        ax.axis('off')
        fig.colorbar(im)

        if field not in ['entropy', 's', 'S'] and normed is True:
            im.set_clim(-max(abs(equator.max()), abs(equator.min())), 
                         max(abs(equator.max()), abs(equator.min())))

    def avg(self, field='Bphi', levels=16, cm='RdYlBu_r', normed=True,
            vmax=None, vmin=None):
        """
        Plot the azimutal average of a given field.
        """
        if field in ('Vr', 'vr', 'Ur', 'ur'):
            data = self.vp
            label = 'Radial velocity'
        elif field in ('Vphi', 'vphi', 'Uphi', 'uphi', 'up', 'Up', 'Vp', 'vp'):
            data = self.vphi
            if labTex:
                label = r'$V_{\phi}$'
            else:
                label = 'vphi'
        elif field in ('Vs', 'vs'):
            data = self.vs
            label = 'Vs'
        elif field in ('Vz', 'vz'):
            data = self.vz
            label = 'Vz'
        elif field in ('rho'):
            data = self.rho
            if labTex:
                label = r'$\rho$'
            else:
                label = 'rho'
        elif field in ('Cr', 'cr'):
            vp = self.vphi.copy()-self.vphi.mean(axis=0) # convective vp
            data =  self.vs * vp
            denom = N.sqrt(N.mean(self.vs**2, axis=0)* N.mean(vp**2, axis=0))
            if labTex:
                label = r'$\langle v_s v_\phi\rangle$'
            else:
                label = 'vs vphi'

        th = N.linspace(0., N.pi, 128)

        if field not in ('Cr', 'cr'):
            phiavg = data.mean(axis=0)
        else:
            mask = N.where(denom == 0, 1, 0)
            phiavg = data.mean(axis=0)/(denom+mask)
            m1 = N.sqrt(self.S**2+self.Z**2) >= self.ri
            m2 = N.sqrt(self.S**2+self.Z**2) <= self.ro
            m3 = self.S <= self.ri
            m4 = self.S >= self.ri
            print('Correlation', phiavg[m1*m2].mean())
            print('Correlation out TC', phiavg[m1*m2*m4].mean())
            print('Correlation in TC', phiavg[m1*m2*m3].mean())

        phiavg = cut(phiavg, vmax, vmin)

        fig = P.figure(figsize=(5.5, 8))
        ax = fig.add_subplot(111, frameon=False)
        cmap = P.get_cmap(cm)
        im = ax.contourf(self.S, self.Z, phiavg, levels, cmap=cmap)
        ax.plot(self.ri*N.sin(th), self.ri*N.cos(th), 'k-')
        ax.plot(self.ro*N.sin(th), self.ro*N.cos(th), 'k-')
        ax.plot([0., 0], [self.ri, self.ro], 'k-')
        ax.plot([0., 0], [-self.ri, -self.ro], 'k-')
        ax.set_title(label, fontsize=24)
        ax.axis('off')
        fig.colorbar(im)

        if field not in ['entropy', 's', 'S'] and normed is True:
            im.set_clim(-max(abs(phiavg.max()), abs(phiavg.min())), 
                         max(abs(phiavg.max()), abs(phiavg.min())))

    def avgz(self, field='vs', levels=16, cm='RdYlBu_r', normed=True, vmin=None,
             vmax=None, avg=False):
        """
        Plot the z-average of a given field
        """
        phi = N.linspace(0., 2.*N.pi, self.nphi)
        rr, pphi = N.meshgrid(self.radius, phi)
        xx = rr * N.cos(pphi)
        yy = rr * N.sin(pphi)
        if field in ('Vr', 'vr', 'Ur', 'ur'):
            data = self.vphi
            label = 'Radial velocity'
        elif field in ('betaz'):
            betaz = cylZder(self.z, N.log(self.rho))
            data =  self.vz * betaz
            data *= self.vs
            if labTex:
                label = r'$\beta_z u_z$'
            else:
                label = 'betaz uz'
        elif field in ('betas'):
            betas = cylSder(self.radius, N.log(self.rho))
            data =  self.vs * betas
            data *= self.vs
            if labTex:
                label = r'$\beta_s u_s$'
            else:
                label = 'betas us'
        elif field in ('rho'):
            data = self.rho
            if labTex:
                label = r'$\rho$'
            else:
                label = 'rho'
        elif field in ('anel'):
            vp = self.vphi.copy()-self.vphi.mean(axis=0) # convective vp
            betas = cylSder(self.radius, N.log(self.rho))
            betaz = cylZder(self.z, N.log(self.rho))
            data = self.vs * betas + self.vz * betaz
            data1 = cylSder(self.radius, self.vphi*self.S)-phideravg(self.vs, self.minc)
            mask = N.where(self.S == 0, 1, 0)
            data1 = data1/(self.S+mask)
            data *= data1
            if labTex:
                label = r'$\beta u_r$'
            else:
                label = 'beta ur'
        elif field in ('vortz'):
            data = cylSder(self.radius, self.vphi*self.S)-phideravg(self.vs, self.minc)
            mask = N.where(self.S == 0, 1, 0)
            data = data/(self.S+mask)
            if labTex:
                label = r'$\omega_z$'
            else:
                label = 'omegaz'
        elif field in ('vopot'):
            data = cylSder(self.radius, self.vphi*self.S)-phideravg(self.vs, self.minc)
            mask = N.where(self.S == 0, 1, 0)
            data = data/(self.S+mask)
            data = data-2./self.ek*N.log(self.rho)
            label = 'vopot'
        elif field in ('Vphi', 'vphi', 'Uphi', 'uphi', 'up', 'Up', 'Vp', 'vp'):
            data = self.vphi
            if labTex:
                label = r'$V_{\phi}$'
            else:
                label = 'vphi'
        elif field in ('Vs', 'vs'):
            data = self.vs
            if labTex:
                label = r'$v_s$'
            else:
                label = 'vs'
        elif field in ('Vz', 'vz'):
            data = self.vz
            if labTex:
                label = r'$v_z$'
            else:
                label = 'vz'
        elif field in ('vpc'):
            data = self.vphi.copy()-self.vphi.mean(axis=0) # convective vp
            if labTex:
                label = r'$v_p$ conv'
            else:
                label = 'vphi conv'
        elif field in ('Cr', 'cr'):
            vp = self.vphi.copy()-self.vphi.mean(axis=0) # convective vp
            data =  self.rho * self.vs * vp
            denom = N.zeros((self.npI, self.ns), dtype='f')
            if labTex:
                label = r'$\langle \rho v_s v_\phi\rangle$'
            else:
                label = 'rho vs vphi'
        elif field in ('reynolds'):
            vp = self.vphi.copy()-self.vphi.mean(axis=0) # convective vp
            phi = N.linspace(0., 2.*N.pi, self.npI)
            data =  self.rho * self.vs * vp
            if labTex:
                label = r'$\rho v_s v_\phi$'
            else:
                label = 'rho vs vphi'
        elif field in ('vsvp'):
            vp = self.vphi.copy()-self.vphi.mean(axis=0) # convective vp
            phi = N.linspace(0., 2.*N.pi, self.npI)
            data =  self.vs * vp
            if labTex:
                label = r'$v_s v_\phi$'
            else:
                label = 'vs vphi'
        elif field in ('vrvs'):
            th2D = N.arctan2(self.S, self.Z)
            vr = self.vs * N.sin(th2D) + self.vz * N.cos(th2D)
            phi = N.linspace(0., 2.*N.pi, self.npI)
            data =  self.vs * vr
            denom = N.zeros((self.npI, self.ns), dtype='f')
            if labTex:
                label = r'$\rho v_s v_r$'
            else:
                label = 'rho vs vr'
        elif field in ('dvz'):
            data =  cylZder(self.z, self.vz)
            data1 = cylSder(self.radius, self.vphi*self.S)-phideravg(self.vs, self.minc)
            mask = N.where(self.S == 0, 1, 0)
            data1 = data1/(self.S+mask)
            data *= data1
            if labTex:
                label = r'$\partial v_z/\partial z$'
            else:
                label = 'dvz/dz'
        elif field in ('balance'):
            if labTex:
                label = r'$\partial v_z/\partial z+\beta v_r$'
            else:
                label = 'dvz/dz + beta*vr'
            data =  cylZder(self.z, self.vz)
            betas = cylSder(self.radius, N.log(self.rho))
            betaz = cylZder(self.z, N.log(self.rho))
            data1 = self.vs * betas + self.vz * betaz
            data += data1
            data2 = cylSder(self.radius, self.vphi*self.S)-phideravg(self.vs, self.minc)
            mask = N.where(self.S == 0, 1, 0)
            data2 = data2/(self.S+mask)
            data *= data2

        equator = N.zeros((self.npI, self.ns), dtype='f')
        for i, rad in enumerate(self.radius):
            if rad <= self.ri:
                zo = N.sqrt(self.ro**2-rad**2) 
                zi = N.sqrt(self.ri**2-rad**2) 
                m1 = abs(self.z) <= zo
                m2 = abs(self.z) >= zi
                equator[:, i] = data[:, m1*m2, i].mean(axis=1)
                if field  in ('Cr', 'cr'):
                    denom[:, i] = N.sqrt( \
                   N.mean(self.rho[:, m1*m2, i]*self.vs[:, m1*m2, i]**2, axis=1)\
                 * N.mean(self.rho[:, m1*m2, i]*vp[:, m1*m2, i]**2, axis=1))
                elif field in ('vrvs'):
                    denom[:, i] = N.sqrt( \
                   N.mean(vr[:, m1*m2, i]**2, axis=1)\
                 * N.mean(self.vs[:, m1*m2, i]**2, axis=1))
            elif rad > self.ri and rad < self.ro:
                zo = N.sqrt(self.ro**2-rad**2) 
                m1 = self.z >= -zo
                m2 = self.z <= zo
                equator[:, i] = data[:, m1*m2, i].mean(axis=1)
                if field  in ('Cr', 'cr'):
                    denom[:, i] = N.sqrt( \
                   N.mean(self.rho[:, m1*m2, i]*self.vs[:, m1*m2, i]**2, axis=1)\
                 * N.mean(self.rho[:, m1*m2, i]*vp[:, m1*m2, i]**2, axis=1))
                elif field in ('vrvs'):
                    denom[:, i] = N.sqrt( \
                   N.mean(vr[:, m1*m2, i]**2, axis=1)\
                 * N.mean(self.vs[:, m1*m2, i]**2, axis=1))
        if field  in ('Cr', 'cr', 'vrvs'):
            mask = N.where(denom == 0, 1, 0)
            equator /= (denom+mask)

        equator = cut(equator, vmax, vmin)
        equator = symmetrize(equator, self.minc)


        fig = P.figure(figsize=(8.25, 6))
        ax = fig.add_subplot(111, frameon=False)
        cmap = P.get_cmap(cm)
        im = ax.contourf(xx, yy, equator, levels, cmap=cmap)
        ax.plot(self.ri * N.cos(phi), self.ri*N.sin(phi), 'k-')
        ax.plot(self.ro * N.cos(phi), self.ro*N.sin(phi), 'k-')
        ax.set_title(label, fontsize=24)
        ax.axis('off')
        fig.colorbar(im)
        if avg:
            fig = P.figure()
            ax = fig.add_subplot(111)
            if field in ('vphi'):
                dat  = N.mean(equator, axis=0)
                dat = dat[:-1]
                ax.plot(self.radius[:-1], dat)
            else:
                ax.plot(self.radius, N.mean(equator, axis=0))
            ax.set_xlabel('Radius', fontsize=18)
            ax.set_xlim(0, self.radius.max())

        if field not in ['entropy', 's', 'S'] and normed is True:
            im.set_clim(-max(abs(equator.max()), abs(equator.min())), 
                         max(abs(equator.max()), abs(equator.min())))

    def slice(self, field='Bphi', lon_0=0., levels=16, cm='RdYlBu_r', 
              normed=True):
        """
        Plot an azimuthal slice of a given field.
        """
        if field in ('Vr', 'vr', 'Ur', 'ur'):
            data = self.vp
            label = 'Radial velocity'
        elif field in ('Vphi', 'vphi', 'Uphi', 'uphi', 'up', 'Up', 'Vp', 'vp'):
            data = self.vphi
            if labTex:
                label = r'$v_{phi}$'
            else:
                label = 'vphi'
        elif field in ('Vs', 'vs'):
            data = self.vs
            if labTex:
                label = r'$v_s$'
            else:
                label = 'vs'
        elif field in ('Vz', 'vz'):
            data = self.vz
            if labTex:
                label = r'$v_z$'
            else:
                label = 'vz'

        data = symmetrize(data, self.minc)

        th = N.linspace(-N.pi/2, N.pi/2, 128)
        phi = N.linspace(0., 360, self.nphi)

        lon_0 = N.asarray(lon_0)
        
        cmap = P.get_cmap(cm)

        if len(lon_0) > 1:
            fig = P.figure(figsize=(3.5*len(lon_0), 5.1))
            for k, lon in enumerate(lon_0):
                ind = N.nonzero(N.where(abs(phi-lon) \
                                == min(abs(phi-lon)), 1, 0))
                indPlot = ind[0][0]
                phislice = data[indPlot, ...]
                ax = fig.add_subplot(1,len(lon_0),k+1, frameon=False)

                im = ax.contourf(self.S, self.Z, phislice, levels, cmap=cmap)
                ax.plot(self.ro*N.cos(th), self.ro*N.sin(th), 'k-')
                ax.plot(self.ri*N.cos(th), self.ri*N.sin(th), 'k-')
                ax.plot([0., 0], [self.ri, self.ro], 'k-')
                ax.plot([0., 0], [-self.ri, -self.ro], 'k-')
                ax.axis('off')
                ax.set_title(label+r' $%i^\circ$' % lon)
                #fig.colorbar(im, orientation='horizontal')

        else:
            ind = N.nonzero(N.where(abs(phi-lon_0[0]) \
                            == min(abs(phi-lon_0[0])), 1, 0))
            indPlot = ind[0][0]
            phislice = data[indPlot, ...]

            fig = P.figure(figsize=(5.5, 8))
            ax = fig.add_subplot(111, frameon=False)
            im = ax.contourf(self.S, self.Z, phislice, levels, cmap=cmap)
            ax.plot(self.ro*N.cos(th), self.ro*N.sin(th), 'k-')
            ax.plot(self.ri*N.cos(th), self.ri*N.sin(th), 'k-')
            ax.plot([0., 0], [self.ri, self.ro], 'k-')
            ax.plot([0., 0], [-self.ri, -self.ro], 'k-')
            ax.set_title(label, fontsize=24)
            ax.axis('off')
            fig.colorbar(im)

        if field not in ['entropy', 's', 'S'] and normed is True:
            im.set_clim(-max(abs(phislice.max()), abs(phislice.min())), 
                         max(abs(phislice.max()), abs(phislice.min())))




if __name__ == '__main__':
    c = Cyl(lastvar=1)
    c.equat(field='vs', normed=False)
    P.show()
