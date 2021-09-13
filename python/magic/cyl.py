# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from .libmagic import anelprof, cylSder, cylZder, phideravg, symmetrize, \
                      progressbar
from .plotlib import cut
from magic import MagicGraph, MagicSetup
from magic.setup import labTex, buildSo
from scipy.ndimage import map_coordinates
from scipy.interpolate import interp1d
import os
import pickle

if buildSo:
    try:
        from magic.cylavg import *
        zavgMode = 'f2py'
    except ImportError:
        zavgMode = 'python'
else:
    zavgMode = 'python'

def sph2cyl_plane(data, rad, ns, nz):
    """
    This function extrapolates a phi-slice of a spherical shell on
    a cylindrical grid

    >>> # Read G_1.test
    >>> gr = MagicGraph(ivar=1, tag='test')
    >>> # phi-average v_\phi and s
    >>> vpm = gr.vphi.mean(axis=0)
    >>> sm = gr.entropy.mean(axis=0)
    >>> # Interpolate on a cylindrical grid
    >>> Z, S, outputs = sph2cyl_plane([vpm, sm], gr.radius, 512, 1024)
    >>> vpm_cyl, sm_cyl = outputs

    :param data: a list of 2-D arrays [(ntheta, nr), (ntheta, nr), ...]
    :type data: list(numpy.ndarray)
    :param rad: radius
    :type rad: numpy.ndarray
    :param ns: number of grid points in s direction
    :type ns: int
    :param nz: number of grid points in z direction
    :type nz: int
    :returns: a python tuple that contains two numpy.ndarray and a list (S,Z,output).
              S[nz,ns] is a meshgrid that contains the radial coordinate.
              Z[nz,ns] is a meshgrid that contains the vertical coordinate.
              output=[arr1[nz,ns], ..., arrN[nz,ns]] is a list of the interpolated
              array on the cylindrical grid.
    :rtype: tuple
    """
    ntheta, nr = data[0].shape

    radius = rad[::-1]

    theta = np.linspace(0., np.pi, ntheta)

    Z, S = np.mgrid[-radius.max():radius.max():nz*1j,0:radius.max():ns*1j]

    new_r = np.sqrt(S**2+Z**2).ravel()
    new_theta = np.arctan2(S, Z).ravel()
    ir = interp1d(radius, np.arange(len(radius)), bounds_error=False)
    it = interp1d(theta, np.arange(len(theta)), bounds_error=False)

    new_ir = ir(new_r)
    new_it = it(new_theta)
    new_ir[new_r > radius.max()] = len(radius)-1.
    new_ir[new_r < radius.min()] = 0.

    coords = np.array([new_it, new_ir])

    output = []
    for dat in data:
        dat_cyl = map_coordinates(dat[:, ::-1], coords, order=3)
        dat_cyl[new_r > radius.max()] = 0.
        dat_cyl[new_r < radius.min()] = 0.
        dat_cyl = dat_cyl.reshape((nz, ns))
        output.append(dat_cyl)

    return Z, S, output

if zavgMode == 'f2py':

    def zavg(input, radius, ns, minc, save=True, filename='vp.pickle',
             normed=True):
        """
        This function computes a z-integration of a list of input arrays 
        (on the spherical grid). This works well for 2-D (phi-slice) arrays.
        In case of 3-D arrays, only one element is allowed (too demanding
        otherwise).

        :param input: a list of 2-D or 3-D arrays
        :type input: list(numpy.ndarray)
        :param radius: spherical radius
        :type radius: numpy.ndarray
        :param ns: radial resolution of the cylindrical grid (nz=2*ns)
        :type ns: int
        :param minc: azimuthal symmetry
        :type minc: int
        :param save: a boolean to specify if one wants to save the outputs into
                     a pickle (default is True)
        :type save: bool
        :param filename: name of the output pickle when save=True
        :type filename: str
        :param normed: a boolean to specify if ones wants to simply integrate 
                       over z or compute a z-average (default is True: average)
        :type normed: bool
        :returns: a python tuple that contains two numpy.ndarray and a list
                  (height,cylRad,output) height[ns] is the height of the
                  spherical shell for all radii. cylRad[ns] is the cylindrical
                  radius. output=[arr1[ns], ..., arrN[ns]] contains
                  the z-integrated output arrays.
        :rtype: tuple
        """
        ro = radius[0]
        ri = radius[-1]
        cylRad = np.linspace(ro, 0., ns)
        if len(input[0].shape) == 3:
            ntheta = input[0].shape[1]
        elif len(input[0].shape) == 2:
            ntheta = input[0].shape[0]
        theta = np.linspace(0., np.pi, ntheta)

        height = np.zeros_like(cylRad)
        height[cylRad >= ri] = 2.*np.sqrt(ro**2-cylRad[cylRad >= ri]**2)
        height[cylRad < ri] = 2.*(np.sqrt(ro**2-cylRad[cylRad < ri]**2)
                                  -np.sqrt(ri**2-cylRad[cylRad < ri]**2))

        if len(input[0].shape) == 3:
            nphi = input[0].shape[0]
            phi = np.linspace(0., 2.*np.pi/minc, nphi)
            output = np.zeros((nphi, ns), dtype=input[0].dtype)
            for iphi in progressbar(range(nphi)):
                output[iphi, :] = cylmean(input[0][iphi, ...], radius, cylRad,
                                          theta)
                if not normed:
                    output[iphi, :] *= height

            if save:
                nphi, ntheta, nr = input[0].shape
                file = open(filename, 'wb')
                pickle.dump([cylRad, phi, output], file) # cylindrical average
                pickle.dump([radius, phi, input[0][:, ntheta/2, :]], file) # equatorial cut
                file.close()
            return height, cylRad, phi, output

        elif len(input[0].shape) == 2:
            output = []
            for dat in input:
                outIntZ = np.zeros(ns, dtype=dat.dtype)
                outIntZ = cylmean(dat, radius, cylRad, theta)
                if not normed:
                    outIntZ *= height
                output.append(outIntZ)

            if save:
                file = open(filename, 'wb')
                pickle.dump([radius,  cylRad, height], file) # cylindrical average
                for k,out in enumerate(output):
                    pickle.dump(out, file) # cylindrical average
                    ntheta, nr = input[k].shape
                    pickle.dump(input[k][ntheta//2, :], file) # equatorial cut
                file.close()

            return height, cylRad, output

else:

    def zavg(input, radius, ns, minc, save=True, filename='vp.pickle',
             normed=True):
        """
        This function computes a z-integration of a list of input arrays 
        (on the spherical grid). This works well for 2-D (phi-slice) 
        arrays. In case of 3-D arrays, only one element is allowed
        (too demanding otherwise).

        :param input: a list of 2-D or 3-D arrays
        :type input: list(numpy.ndarray)
        :param radius: spherical radius
        :type radius: numpy.ndarray
        :param ns: radial resolution of the cylindrical grid (nz=2*ns)
        :type ns: int
        :param minc: azimuthal symmetry
        :type minc: int
        :param save: a boolean to specify if one wants to save the outputs into
                     a pickle (default is True)
        :type save: bool
        :param filename: name of the output pickle when save=True
        :type filename: str
        :param normed: a boolean to specify if ones wants to simply integrate
                       over z or compute a z-average (default is True: average)
        :type normed: bool
        :returns: a python tuple that contains two numpy.ndarray and a
                  list (height,cylRad,output) height[ns] is the height of the
                  spherical shell for all radii. cylRad[ns] is the cylindrical
                  radius. output=[arr1[ns], ..., arrN[ns]] contains
                  the z-integrated output arrays.
        :rtype: tuple
        """
        nz = 2*ns
        ro = radius[0]
        ri = radius[-1]
        z = np.linspace(-ro, ro, nz)
        cylRad = np.linspace(0., ro, ns)
        cylRad = cylRad[1:-1]

        height = np.zeros_like(cylRad)
        height[cylRad >= ri] = 2.*np.sqrt(ro**2-cylRad[cylRad >= ri]**2)
        height[cylRad < ri] = 2.*(np.sqrt(ro**2-cylRad[cylRad < ri]**2)
                                  -np.sqrt(ri**2-cylRad[cylRad < ri]**2))

        if len(input[0].shape) == 3:
            nphi = input[0].shape[0]
            phi = np.linspace(0., 2.*np.pi/minc, nphi)
            output = np.zeros((nphi, ns-2), dtype=input[0].dtype)
            for iphi in progressbar(range(nphi)):
                Z, S, out2D = sph2cyl_plane([input[0][iphi, ...]], radius, ns,
                                            nz)
                S = S[:, 1:-1]
                Z = Z[:, 1:-1]
                output[iphi, :] = np.trapz(out2D[0][:, 1:-1], z, axis=0)
                if normed:
                    output[iphi, :] /= height

            if save:
                nphi, ntheta, nr = input[0].shape
                file = open(filename, 'wb')
                pickle.dump([cylRad, phi, output], file) # cylindrical average
                pickle.dump([radius, phi, input[0][:, ntheta//2, :]], file) # equatorial cut
                file.close()
            return height, cylRad, phi, output
        elif len(input[0].shape) == 2:
            Z, S, out2D = sph2cyl_plane(input, radius, ns, nz)
            S = S[:, 1:-1]
            Z = Z[:, 1:-1]
            output = []
            outIntZ = np.zeros((ns-2), dtype=input[0].dtype)
            for k,out in enumerate(out2D):
                outIntZ = np.trapz(out[:, 1:-1], z, axis=0)
                if normed:
                    outIntZ /= height
                output.append(outIntZ)

            if save:
                file = open(filename, 'wb')
                pickle.dump([radius,  cylRad, height], file) # cylindrical average
                for k,out in enumerate(output):
                    pickle.dump(out, file) # cylindrical average
                    ntheta, nr = input[k].shape
                    pickle.dump(input[k][ntheta//2, :], file) # equatorial cut
                file.close()

            return height, cylRad, output


def sph2cyl(g, ns=None, nz=None):
    """
    This function interpolates the three flow (or magnetic field)
    component of a :ref:`G_#.TAG <secGraphFile>` file
    on a cylindrical grid of size (ns, nz).

    .. warning:: This might be really slow!

    :param g: input graphic output file
    :type g: :py:class:`magic.MagicGraph`
    :param ns: number of grid points in the radial direction
    :type ns: int
    :param nz: number of grid points in the vertical direction
    :type nz: int
    :returns: a python tuple of five numpy.ndarray (S,Z,vs,vp_cyl,vz).
              S[nz,ns] is a meshgrid that contains the radial coordinate.
              Z[nz,ns] is a meshgrid that contains the vertical coordinate.
              vs[nz,ns] is the radial component of the velocity (or magnetic
              field), vp_cyl[nz,ns] the azimuthal component and vz[nz,ns] the
              vertical component.
    :rtype: tuple
    """
    if ns is None or nz is None:
        ns = g.nr ; nz = 2*ns

    theta = np.linspace(0., np.pi, g.ntheta)
    radius = g.radius[::-1]

    Z, S = np.mgrid[-radius.max():radius.max():nz*1j,0:radius.max():ns*1j]

    new_r = np.sqrt(S**2+Z**2).ravel()
    new_theta = np.arctan2(S, Z).ravel()
    ir = interp1d(radius, np.arange(len(radius)), bounds_error=False)
    it = interp1d(theta, np.arange(len(theta)), bounds_error=False)

    new_ir = ir(new_r)
    new_it = it(new_theta)
    new_ir[new_r > radius.max()] = len(radius)-1.
    new_ir[new_r < radius.min()] = 0.

    coords = np.array([new_it, new_ir])

    vr_cyl = np.zeros((g.npI, nz, ns), dtype=g.vr.dtype)
    vp_cyl = np.zeros_like(vr_cyl)
    vt_cyl = np.zeros_like(vr_cyl)
    for k in progressbar(range(g.npI)):
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

    th3D = np.zeros((g.npI, nz, ns), dtype=g.vr.dtype)
    for i in range(g.npI):
        th3D[i, ...] = np.arctan2(S, Z)
    vs = vr_cyl * np.sin(th3D) + vt_cyl * np.cos(th3D)
    vz = vr_cyl * np.cos(th3D) - vt_cyl * np.sin(th3D)

    return S, Z, vs, vp_cyl, vz


class Cyl(MagicSetup):
    """
    This class allows to extrapolate a given :ref:`graphic file <secGraphFile>`
    on a cylindrical grid. Once done, the extrapolated file is stored in
    a python.pickle file. It is then possible to display 2-D cuts of the extrapolated
    arrays (radial cuts, phi-averages, equatorial cuts, z-averages and phi-slices)

    .. warning:: This process is actually **very demanding** and it might take a lot
                 of time to extrapolate the G_#.TAG file. Be careful when choosing the
                 input value of ns!

    >>> # Extrapolate the G file to the cylindrical grid (ns=128, nz=2*ns)
    >>> c = Cyl(ivar=1, ns=128)
    >>> # Radial cut of v_r
    >>> c.surf(field='vr', r=0.8)
    >>> # Vertical average of B_\phi
    >>> c.avgz(field='Bphi', cm='seismic', levels=33)
    >>> # Azimuthal average of v_\phi
    >>> c.avg(field='Bphi')
    >>> # Equatorial cut of of v_theta
    >>> c.equat(field='vtheta')
    """

    def __init__(self, ivar=1, datadir='.', ns=None):
        """
        :param ivar: the number of the Graphic file
        :type ivar: int
        :param datadir: working directory
        :type datadir: str
        :param ns: number of grid points in the radial direction
        :type ns: int
        """
        MagicSetup.__init__(self, datadir)

        self.datadir = datadir

        filename = '{}G_{}.{}'.format('cyl', ivar, self.tag)
        if not os.path.exists(filename):
            print("sph2cyl...")
            gr = MagicGraph(ivar=ivar, datadir=self.datadir)
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
        self.radius = np.linspace(0., self.ro, self.ns)
        temp0, rho0, beta0 = anelprof(np.linspace(self.ro, self.ri, self.ns),
                                     self.strat, self.polind)
        rho = np.zeros((self.nphi/2, self.ns), dtype=self.vr.dtype)
        beta = np.zeros_like(rho)
        for i in range(self.nphi/2):
            rho[i, :] = rho0
            beta[i, :] = beta0
        Z, S, [rho, beta] = sph2cyl_plane([rho,beta],
                                 np.linspace(self.ro, self.ri, self.ns),
                                 self.ns, self.nz)
        self.rho = np.zeros_like(self.vs)
        self.beta = np.zeros_like(self.vs)
        for i in range(self.npI):
            self.rho[i, ...] = rho
            self.beta[i, ...] = beta
        self.z = np.linspace(-self.ro, self.ro, self.nz)

    def surf(self, field='Bphi', r=0.85, vmin=None, vmax=None,
             levels=16, cm='RdYlBu_r', normed=True, figsize=None):
        """
        Plot the surface distribution of an input field at a given
        input radius (normalised by the outer boundary radius).

            >>> c = Cyl(ns=65)
            >>> # Surface plot of B_\phi from -10 to 10
            >>> c.surf(field='Bphi', r=0.6, vmin=-10, vmax=10, levels=65)

        :param field:  name of the input field
        :type field: str
        :param r: radial level (normalised to the outer boundary radius)
        :type r: float
        :param levels: number of contour levels
        :type levels: int
        :param cm: name of the color map
        :type cm: str
        :param normed: when set to True, the contours are normalised fro
                       -max(field), max(field)
        :type normed: bool
        :param vmin: truncate the contour levels to values > vmin
        :type vmin: float
        :param vmax: truncate the contour levels to values < vmax
        :type vmax: float
        """
        r /= (1-self.ri/self.ro) # as we give a normalised radius
        ind = np.nonzero(np.where(abs(self.radius-r) \
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

        phi = np.linspace(0., 2.*np.pi, self.nphi)

        data[..., indPlot] = cut(data[..., indPlot], vmax, vmin)
        data = symmetrize(data, self.minc)

        cmap = plt.get_cmap(cm)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        im = ax.contourf(phi, self.z, data[..., indPlot].T, levels, cmap=cmap,
                        aa=True)
        rad = self.radius[indPlot] * (1. - self.ri/self.ro)
        if labTex:
            ax.set_xlabel(r'$\phi$', fontsize=18)
            ax.set_ylabel(r'$z$', fontsize=18)
            ax.set_title('{}: $r/r_o$ = {:.3f}'.format(label, rad), fontsize=24)
        else:
            ax.set_xlabel('phi', fontsize=18)
            ax.set_ylabel('z', fontsize=18)
            ax.set_title('{}: r/ro = {:.3f}'.format(label, rad), fontsize=24)
        cbar = plt.colorbar(im)

        if field not in ['entropy', 's', 'S'] and normed is True:
            im.set_clim(-max(abs(data[..., indPlot].max()),
                             abs(data[..., indPlot].min())),
                         max(abs(data[..., indPlot].max()),
                             abs(data[..., indPlot].min())))

    def equat(self, field='vs', levels=16, cm='RdYlBu_r', normed=True, vmax=None,
              vmin=None):
        """
        Plot an input field in the equatorial plane.

            >>> c = Cyl(ns=65)
            >>> # Equatorial cut of v_\phi
            >>> c.equat(field='vphi', cm='seismic', levels=33)

        :param field:  name of the input field
        :type field: str
        :param levels: number of contour levels
        :type levels: int
        :param cm: name of the color map
        :type cm: str
        :param normed: when set to True, the contours are normalised fro
                       -max(field), max(field)
        :type normed: bool
        :param vmin: truncate the contour levels to values > vmin
        :type vmin: float
        :param vmax: truncate the contour levels to values < vmax
        :type vmax: float
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
            betas = cylSder(self.radius, np.log(self.rho))
            betaz = cylZder(self.z, np.log(self.rho))
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

        phi = np.linspace(0., 2.*np.pi, self.nphi)
        rr, pphi = np.meshgrid(self.radius, phi)
        xx = rr * np.cos(pphi)
        yy = rr * np.sin(pphi)

        fig = plt.figure(figsize=(8.25, 6))
        ax = fig.add_subplot(111, frameon=False)
        cmap = plt.get_cmap(cm)
        im = ax.contourf(xx, yy, equator, levels, cmap=cmap)
        ax.plot(self.ri * np.cos(phi), self.ri*np.sin(phi), 'k-')
        ax.plot(self.ro * np.cos(phi), self.ro*np.sin(phi), 'k-')
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

            >>> c = Cyl(ns=65)
            >>> # Azimuthal average of B_r
            >>> c.avg(field='Br', cm='seismic', levels=33)

        :param field: name of the input field
        :type field: str
        :param levels: number of contour levels
        :type levels: int
        :param cm: name of the color map
        :type cm: str
        :param normed: when set to True, the contours are normalised fro
                       -max(field), max(field)
        :type normed: bool
        :param vmin: truncate the contour levels to values > vmin
        :type vmin: float
        :param vmax: truncate the contour levels to values < vmax
        :type vmax: float
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
            denom = np.sqrt(np.mean(self.vs**2, axis=0)* np.mean(vp**2, axis=0))
            if labTex:
                label = r'$\langle v_s v_\phi\rangle$'
            else:
                label = 'vs vphi'

        th = np.linspace(0., np.pi, 128)

        if field not in ('Cr', 'cr'):
            phiavg = data.mean(axis=0)
        else:
            mask = np.where(denom == 0, 1, 0)
            phiavg = data.mean(axis=0)/(denom+mask)
            m1 = np.sqrt(self.S**2+self.Z**2) >= self.ri
            m2 = np.sqrt(self.S**2+self.Z**2) <= self.ro
            m3 = self.S <= self.ri
            m4 = self.S >= self.ri
            print('Correlation', phiavg[m1*m2].mean())
            print('Correlation out TC', phiavg[m1*m2*m4].mean())
            print('Correlation in TC', phiavg[m1*m2*m3].mean())

        phiavg = cut(phiavg, vmax, vmin)

        fig = plt.figure(figsize=(5.5, 8))
        ax = fig.add_subplot(111, frameon=False)
        cmap = plt.get_cmap(cm)
        im = ax.contourf(self.S, self.Z, phiavg, levels, cmap=cmap)
        ax.plot(self.ri*np.sin(th), self.ri*np.cos(th), 'k-')
        ax.plot(self.ro*np.sin(th), self.ro*np.cos(th), 'k-')
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
        Plot the vertical average of a given field.

            >>> c = Cyl(ns=65)
            >>> # Vertical average of v_s
            >>> c.avg(field='vs', cm='seismic', levels=33)

        :param field: name of the input field
        :type field: str
        :param levels: number of contour levels
        :type levels: int
        :param cm: name of the color map
        :type cm: str
        :param normed: when set to True, the contours are normalised fro
                       -max(field), max(field)
        :type normed: bool
        :param vmin: truncate the contour levels to values > vmin
        :type vmin: float
        :param vmax: truncate the contour levels to values < vmax
        :type vmax: float
        :param avg: when set to True, an additional figure with the phi-average
                    profile is also displayed
        :type avg: bool
        """
        phi = np.linspace(0., 2.*np.pi, self.nphi)
        rr, pphi = np.meshgrid(self.radius, phi)
        xx = rr * np.cos(pphi)
        yy = rr * np.sin(pphi)
        if field in ('Vr', 'vr', 'Ur', 'ur'):
            data = self.vphi
            label = 'Radial velocity'
        elif field in ('betaz'):
            betaz = cylZder(self.z, np.log(self.rho))
            data =  self.vz * betaz
            data *= self.vs
            if labTex:
                label = r'$\beta_z u_z$'
            else:
                label = 'betaz uz'
        elif field in ('betas'):
            betas = cylSder(self.radius, np.log(self.rho))
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
            betas = cylSder(self.radius, np.log(self.rho))
            betaz = cylZder(self.z, np.log(self.rho))
            data = self.vs * betas + self.vz * betaz
            data1 = cylSder(self.radius, self.vphi*self.S)-phideravg(self.vs, self.minc)
            mask = np.where(self.S == 0, 1, 0)
            data1 = data1/(self.S+mask)
            data *= data1
            if labTex:
                label = r'$\beta u_r$'
            else:
                label = 'beta ur'
        elif field in ('vortz'):
            data = cylSder(self.radius, self.vphi*self.S)-phideravg(self.vs, self.minc)
            mask = np.where(self.S == 0, 1, 0)
            data = data/(self.S+mask)
            if labTex:
                label = r'$\omega_z$'
            else:
                label = 'omegaz'
        elif field in ('vopot'):
            data = cylSder(self.radius, self.vphi*self.S)-phideravg(self.vs, self.minc)
            mask = np.where(self.S == 0, 1, 0)
            data = data/(self.S+mask)
            data = data-2./self.ek*np.log(self.rho)
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
            denom = np.zeros((self.npI, self.ns), dtype=vp.dtype)
            if labTex:
                label = r'$\langle \rho v_s v_\phi\rangle$'
            else:
                label = 'rho vs vphi'
        elif field in ('reynolds'):
            vp = self.vphi.copy()-self.vphi.mean(axis=0) # convective vp
            phi = np.linspace(0., 2.*np.pi, self.npI)
            data =  self.rho * self.vs * vp
            if labTex:
                label = r'$\rho v_s v_\phi$'
            else:
                label = 'rho vs vphi'
        elif field in ('vsvp'):
            vp = self.vphi.copy()-self.vphi.mean(axis=0) # convective vp
            phi = np.linspace(0., 2.*np.pi, self.npI)
            data =  self.vs * vp
            if labTex:
                label = r'$v_s v_\phi$'
            else:
                label = 'vs vphi'
        elif field in ('vrvs'):
            th2D = np.arctan2(self.S, self.Z)
            vr = self.vs * np.sin(th2D) + self.vz * np.cos(th2D)
            phi = np.linspace(0., 2.*np.pi, self.npI)
            data =  self.vs * vr
            denom = np.zeros((self.npI, self.ns), dtype=vr.dtype)
            if labTex:
                label = r'$\rho v_s v_r$'
            else:
                label = 'rho vs vr'
        elif field in ('dvz'):
            data =  cylZder(self.z, self.vz)
            data1 = cylSder(self.radius, self.vphi*self.S)-phideravg(self.vs, self.minc)
            mask = np.where(self.S == 0, 1, 0)
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
            betas = cylSder(self.radius, np.log(self.rho))
            betaz = cylZder(self.z, np.log(self.rho))
            data1 = self.vs * betas + self.vz * betaz
            data += data1
            data2 = cylSder(self.radius, self.vphi*self.S)-phideravg(self.vs, self.minc)
            mask = np.where(self.S == 0, 1, 0)
            data2 = data2/(self.S+mask)
            data *= data2

        equator = np.zeros((self.npI, self.ns), dtype=self.vs.dtype)
        for i, rad in enumerate(self.radius):
            if rad <= self.ri:
                zo = np.sqrt(self.ro**2-rad**2)
                zi = np.sqrt(self.ri**2-rad**2)
                m1 = abs(self.z) <= zo
                m2 = abs(self.z) >= zi
                equator[:, i] = data[:, m1*m2, i].mean(axis=1)
                if field  in ('Cr', 'cr'):
                    denom[:, i] = np.sqrt( \
                   np.mean(self.rho[:, m1*m2, i]*self.vs[:, m1*m2, i]**2, axis=1)\
                 * np.mean(self.rho[:, m1*m2, i]*vp[:, m1*m2, i]**2, axis=1))
                elif field in ('vrvs'):
                    denom[:, i] = np.sqrt( \
                   np.mean(vr[:, m1*m2, i]**2, axis=1)\
                 * np.mean(self.vs[:, m1*m2, i]**2, axis=1))
            elif rad > self.ri and rad < self.ro:
                zo = np.sqrt(self.ro**2-rad**2)
                m1 = self.z >= -zo
                m2 = self.z <= zo
                equator[:, i] = data[:, m1*m2, i].mean(axis=1)
                if field  in ('Cr', 'cr'):
                    denom[:, i] = np.sqrt( \
                   np.mean(self.rho[:, m1*m2, i]*self.vs[:, m1*m2, i]**2, axis=1)\
                 * np.mean(self.rho[:, m1*m2, i]*vp[:, m1*m2, i]**2, axis=1))
                elif field in ('vrvs'):
                    denom[:, i] = np.sqrt( \
                   np.mean(vr[:, m1*m2, i]**2, axis=1)\
                 * np.mean(self.vs[:, m1*m2, i]**2, axis=1))
        if field  in ('Cr', 'cr', 'vrvs'):
            mask = np.where(denom == 0, 1, 0)
            equator /= (denom+mask)

        equator = cut(equator, vmax, vmin)
        equator = symmetrize(equator, self.minc)


        fig = plt.figure(figsize=(8.25, 6))
        ax = fig.add_subplot(111, frameon=False)
        cmap = plt.get_cmap(cm)
        im = ax.contourf(xx, yy, equator, levels, cmap=cmap)
        ax.plot(self.ri * np.cos(phi), self.ri*np.sin(phi), 'k-')
        ax.plot(self.ro * np.cos(phi), self.ro*np.sin(phi), 'k-')
        ax.set_title(label, fontsize=24)
        ax.axis('off')
        fig.colorbar(im)
        if avg:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            if field in ('vphi'):
                dat  = np.mean(equator, axis=0)
                dat = dat[:-1]
                ax.plot(self.radius[:-1], dat)
            else:
                ax.plot(self.radius, np.mean(equator, axis=0))
            ax.set_xlabel('Radius', fontsize=18)
            ax.set_xlim(0, self.radius.max())

        if field not in ['entropy', 's', 'S'] and normed is True:
            im.set_clim(-max(abs(equator.max()), abs(equator.min())),
                         max(abs(equator.max()), abs(equator.min())))

    def slice(self, field='Bphi', lon_0=0., levels=16, cm='RdYlBu_r',
              normed=True):
        """
        Plot an azimuthal slice of a given field.

            >>> c = Cyl(ns=65)
            >>> # Slices of v_r at 30 and 60 degrees
            >>> c.slice(field='vr', lon_0=[30, 60])

        :param field: name of the input field
        :type field: str
        :param lon_0: the longitude of the slice in degrees, or a list of longitudes
        :type lon_0: float or list
        :param levels: number of contour levels
        :type levels: int
        :param cm: name of the color map
        :type cm: str
        :param normed: when set to True, the contours are normalised fro
                       -max(field), max(field)
        :type normed: bool
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

        th = np.linspace(-np.pi/2, np.pi/2, 128)
        phi = np.linspace(0., 360, self.nphi)

        lon_0 = np.asarray(lon_0)

        cmap = plt.get_cmap(cm)

        if len(lon_0) > 1:
            fig = plt.figure(figsize=(3.5*len(lon_0), 5.1))
            for k, lon in enumerate(lon_0):
                ind = np.nonzero(np.where(abs(phi-lon) \
                                == min(abs(phi-lon)), 1, 0))
                indPlot = ind[0][0]
                phislice = data[indPlot, ...]
                ax = fig.add_subplot(1,len(lon_0),k+1, frameon=False)

                im = ax.contourf(self.S, self.Z, phislice, levels, cmap=cmap)
                ax.plot(self.ro*np.cos(th), self.ro*np.sin(th), 'k-')
                ax.plot(self.ri*np.cos(th), self.ri*np.sin(th), 'k-')
                ax.plot([0., 0], [self.ri, self.ro], 'k-')
                ax.plot([0., 0], [-self.ri, -self.ro], 'k-')
                ax.axis('off')
                ax.set_title(label+r' ${}^\circ$'.format(lon))
                #fig.colorbar(im, orientation='horizontal')

        else:
            ind = np.nonzero(np.where(abs(phi-lon_0[0]) \
                            == min(abs(phi-lon_0[0])), 1, 0))
            indPlot = ind[0][0]
            phislice = data[indPlot, ...]

            fig = plt.figure(figsize=(5.5, 8))
            ax = fig.add_subplot(111, frameon=False)
            im = ax.contourf(self.S, self.Z, phislice, levels, cmap=cmap)
            ax.plot(self.ro*np.cos(th), self.ro*np.sin(th), 'k-')
            ax.plot(self.ri*np.cos(th), self.ri*np.sin(th), 'k-')
            ax.plot([0., 0], [self.ri, self.ro], 'k-')
            ax.plot([0., 0], [-self.ri, -self.ro], 'k-')
            ax.set_title(label, fontsize=24)
            ax.axis('off')
            fig.colorbar(im)

        if field not in ['entropy', 's', 'S'] and normed is True:
            im.set_clim(-max(abs(phislice.max()), abs(phislice.min())),
                         max(abs(phislice.max()), abs(phislice.min())))




if __name__ == '__main__':
    c = Cyl(ivar=1)
    c.equat(field='vs', normed=False)
    plt.show()
