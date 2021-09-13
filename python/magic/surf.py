# -*- coding: utf-8 -*-
from magic import MagicGraph, MagicSetup, MagicRadial
from magic.setup import labTex, defaultCm, defaultLevels
from .libmagic import *
from .plotlib import equatContour, merContour, radialContour
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.integrate import trapz


class Surf:
    """
    This class allows to display the content of a graphic file
    (:ref:`G_#.TAG <secGraphFile>` or G_ave.TAG). It allows to plot
    radial, azimuthal and equatorial cuts as well as phi-averages.

    >>> # To read G_1.test
    >>> s = Surf(ivar=1, ave=False, tag='test')
    >>> # To read the latest G file in the working directory (double precision)
    >>> s = Surf(precision=np.float64)

    >>> # Possible plots
    >>> s.equat(field='vr')
    >>> s.avg(field='vp')
    >>> s.surf(field='entropy', r=0.8)
    >>> s.slice(field='Br', lon_0=[0, 30])
    """

    def __init__(self, ivar=None, datadir='.', vort=False, ave=False, tag=None,
                 precision=np.float32):
        """
        :param ivar: index of the graphic file
        :type ivar: int
        :param ave: when set to True, it tries to read a time-averaged graphic file
        :type ave: bool
        :param tag: TAG suffix extension of the graphic file
        :type tag: str
        :param vort: a boolean to specify whether one wants to compute the 3-D
                     vorticiy components (take care of the memory imprint)
        :type vort: bool
        :param datadir: the working directory
        :type datadir: str
        :param precision: the storage precision of the graphic file (single or
                          double precision). Default is np.float32 (single)
        :type precision: str
        """
        self.precision = precision
        self.datadir = datadir
        self.gr = MagicGraph(ivar=ivar, datadir=self.datadir, ave=ave, tag=tag,
                             precision=self.precision)

        if vort:
            thlin = self.gr.colatitude
            th3D = np.zeros_like(self.gr.vphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]

            s3D = rr3D * np.sin(th3D)
            dtheta = thetaderavg(self.gr.vphi*s3D)
            dr = rderavg(self.gr.vphi*s3D, eta=self.gr.radratio, spectral=True,
                         exclude=False)
            ds = np.sin(th3D)*dr + np.cos(th3D)/rr3D*dtheta
            vs = self.gr.vr * np.sin(th3D) + self.gr.vtheta * np.cos(th3D) # 'vs'
            self.vortz = -1./s3D*phideravg(vs, self.gr.minc)+ds/s3D
            del dr, dtheta, ds, rr3D, th3D, s3D

    def surf(self, field='Bphi', proj='hammer', lon_0=0., r=0.85, vmax=None,
             vmin=None, lat_0=30., levels=defaultLevels, cm=defaultCm, ic=False,
             lon_shift=0, normed=True, cbar=True, tit=True, lines=False):
        """
        Plot the surface distribution of an input field at a given
        input radius (normalised by the outer boundary radius).

           >>> s = Surf()
           >>> # Radial flow component at ``r=0.95 r_o``, 65 contour levels
           >>> s.surf(field='vr', r=0.95, levels=65, cm='seismic')

           >>> # Minimal plot (no cbar, not title)
           >>> s.surf(field='entropyfluct', r=0.6, tit=False, cbar=False)

           >>> # Control the limit of the colormap from -1e3 to 1e3
           >>> s.surf(field='vp', r=1., vmin=-1e3, vmax=1e3, levels=33)

           >>> # If basemap is installed, additional projections are available
           >>> s.surf(field='Br', r=0.95, proj='ortho', lat_0=45, lon_0=45)

        :param field: the name of the field you want to display
        :type field: str
        :param proj: the type of projection. Default is Hammer, in case
                     you want to use 'ortho' or 'moll', then Basemap is
                     required.
        :type proj: str
        :param r: the radius at which you want to display the input
                  data (in normalised units with the radius of the outer boundary)
        :type r: float
        :param levels: the number of levels in the contour
        :type levels: int
        :param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
        :type cm: str
        :param lon_shift: translate map in azimuth (in degrees)
        :type lon_shift: int
        :param lon_0: central azimuth (only used with Basemap)
        :type lon_0: float
        :param lat_0: central latitude (only used with Basemap)
        :type lat_0: float
        :param tit: display the title of the figure when set to True
        :param tit: display the title of the figure when set to True
        :type tit: bool
        :param cbar: display the colorbar when set to True
        :type cbar: bool
        :param lines: when set to True, over-plot solid lines to highlight
                      the limits between two adjacent contour levels
        :type lines: bool
        :param vmax: maximum value of the contour levels
        :type vmax: float
        :param vmin: minimum value of the contour levels
        :type vmin: float
        :param normed: when set to True, the colormap is centered around zero.
                       Default is True, except for entropy/temperature plots.
        :type normed: bool
        :param lines: when set to True, over-plot solid lines to highlight
                      the limits between two adjacent contour levels
        :type lines: bool
        """

        if proj != 'ortho':
            lon_0 = 0.

        if field in ('Vs', 'vs'):
            vr = self.gr.vr
            vt = self.gr.vtheta
            thlin = np.linspace(0., np.pi, self.gr.ntheta)
            th3D = np.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            data = vr * np.sin(th3D) + vt * np.cos(th3D)
            data_ic = None
            if labTex:
                label = r'$v_s$'
            else:
                label = r'vs'
        elif field in ('Vz', 'vz'):
            vr = self.gr.vr
            vt = self.gr.vtheta
            thlin = np.linspace(0., np.pi, self.gr.ntheta)
            th3D = np.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            data = vr * np.cos(th3D) - vt * np.sin(th3D)
            data_ic = None
            if labTex:
                label = r'$v_z$'
            else:
                label = r'vz'
        elif field in ('thu'):
            data = self.gr.vr*(self.gr.entropy-self.gr.entropy.mean(axis=0))
            data_ic = None
            label = 'thu'
        elif field in ('flux'):
            data = rderavg(self.gr.entropy, eta=self.gr.radratio)
            data_ic = None
            label = 'flux'
        elif field in ('mag_pres_force_r'):
            data = -rderavg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, eta=self.gr.radratio)/2.0
            data_ic = None
            label = 'Rad. mag. pres. force'
        elif field in ('mag_pres_force_t'):
            rr3D = np.zeros_like(self.gr.Bphi)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            data = -thetaderavg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, order=2)/rr3D/2.0
            data_ic = None
            label = 'Lati. mag. pres. force'
        elif field in ('mag_pres_force_p'):
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = -phideravg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, self.gr.minc)/(rr3D*np.sin(th3D))/2.0
            data_ic = None
            label = 'Longi. mag. pres. force'
        elif field in ('mag_tens_force_r'):
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = self.gr.Br * rderavg(self.gr.Br, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Br, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Br, self.gr.minc) / np.sin(th3D) / rr3D - \
                   (self.gr.Btheta**2 + self.gr.Bphi**2) / rr3D
            data_ic = None
            label = 'Rad. tens. force'
        elif field in ('mag_tens_force_t'):
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = self.gr.Br * rderavg(self.gr.Btheta, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Btheta, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Btheta, self.gr.minc) / np.sin(th3D) / rr3D + \
                   self.gr.Btheta * self.gr.Br / rr3D - \
                   self.gr.Bphi**2 * np.arctan(th3D) / rr3D
            data_ic = None
            label = 'Lati. tens. force'
        elif field in ('mag_tens_force_p'):
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = self.gr.Br * rderavg(self.gr.Bphi, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Bphi, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Bphi, self.gr.minc) / np.sin(th3D) / rr3D+ \
                   self.gr.Bphi * self.gr.Br / rr3D + \
                   self.gr.Bphi * self.gr.Btheta * np.arctan(th3D) / rr3D
            data_ic = None
            label = 'Longi. tens. force'
        elif field in ('Lorentz_r'):
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = -rderavg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, eta=self.gr.radratio)/2.0 + \
                   self.gr.Br * rderavg(self.gr.Br, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Br, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Br, self.gr.minc) / np.sin(th3D) / rr3D - \
                   (self.gr.Btheta**2 + self.gr.Bphi**2) / rr3D
            data_ic = None
            label = 'Radial Lorentz force'
        elif field in ('Lorentz_t'):
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = -thetaderavg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, order=2)/rr3D/2.0 + \
                   self.gr.Br * rderavg(self.gr.Btheta, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Btheta, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Btheta, self.gr.minc) / np.sin(th3D) / rr3D + \
                   self.gr.Btheta * self.gr.Br / rr3D - \
                   self.gr.Bphi**2 * np.arctan(th3D) / rr3D
            data_ic = None
            label = 'Lati. Lorentz force'
        elif field in ('Lorentz_p'):
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = -phideravg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, self.gr.minc)/(rr3D*np.sin(th3D))/2.0 + \
                   self.gr.Br * rderavg(self.gr.Bphi, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Bphi, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Bphi, self.gr.minc) / np.sin(th3D) / rr3D+ \
                   self.gr.Bphi * self.gr.Br / rr3D + \
                   self.gr.Bphi * self.gr.Btheta * np.arctan(th3D) / rr3D
            data_ic = None
            label = 'Longi. Lorentz force'
        elif field in ('ohm'):
            label = 'Ohmic dissipation/1e6'
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            dphi = 2.*np.pi/self.gr.nphi

            Op = (np.roll(self.gr.Btheta,-1,axis=2)-np.roll(self.gr.Btheta,1,axis=2))/\
                 (np.roll(rr3D,-1,axis=2)-np.roll(rr3D,1,axis=2)) + \
                 self.gr.Btheta/rr3D - \
                 (np.roll(self.gr.Br,-1,axis=1)-np.roll(self.gr.Br,1,axis=1))/\
                 (rr3D*(np.roll(th3D,-1,axis=1)-np.roll(th3D,1,axis=1)))
            Ot = (np.roll(self.gr.Br,-1,axis=0)-np.roll(self.gr.Br,1,axis=0))/\
                  (rr3D*np.sin(th3D)*dphi) - \
                 (np.roll(self.gr.Bphi,-1,axis=2)-np.roll(self.gr.Bphi,1,axis=2))/\
                 (np.roll(rr3D,-1,axis=2)-np.roll(rr3D,1,axis=2)) - \
                 self.gr.Bphi/rr3D
            Or = (np.roll(self.gr.Bphi,-1,axis=1)-np.roll(self.gr.Bphi,1,axis=1))/\
                 (rr3D*(np.roll(th3D,-1,axis=1)-np.roll(th3D,1,axis=1))) + \
                 np.cos(th3D)*self.gr.Bphi/(rr3D*np.sin(th3D)) - \
                 (np.roll(self.gr.Btheta,-1,axis=0)-np.roll(self.gr.Btheta,1,axis=0))/\
                 (rr3D*np.sin(th3D)*dphi)

            Or[:, 0, :] = Or[:, 1, :]
            Or[:, -1, :] = Or[:, -2, :]
            Ot[..., 0] = Ot[..., 1]
            Ot[..., -1] = Ot[..., -2]
            Op[..., 0] = Op[..., 1]
            Op[..., -1] = Op[..., -2]
            data = (Op**2+Ot**2+Or**2)/1e6
            data_ic = None
        elif field in ('vortzfluct'):
            th3D = np.zeros_like(self.gr.vphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            s3D = rr3D*np.sin(th3D)
            dth = thetaderavg((self.gr.vphi-self.gr.vphi.mean(axis=0))*rr3D*np.sin(th3D))
            dr = rderavg((self.gr.vphi-self.gr.vphi.mean(axis=0))*rr3D*np.sin(th3D), \
                         eta=self.gr.radratio, spectral=True, exclude=False)
            ds = np.sin(th3D)*dr + np.cos(th3D)/rr3D*dth
            data = -1./(rr3D*np.sin(th3D))*phideravg(self.gr.vr*np.sin(th3D)+self.gr.vtheta*np.cos(th3D), self.gr.minc)+ds/(rr3D*np.sin(th3D))

            del dr, dth, ds, rr3D, th3D

            data_ic = None

            if labTex:
                label = r"$\omega_z'$"
            else:
                label = 'vortzfluct'
        elif field in ('vortz'):
            th3D = np.zeros_like(self.gr.vphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            s3D = rr3D*np.sin(th3D)
            dth = thetaderavg(self.gr.vphi*rr3D*np.sin(th3D))
            dr = rderavg(self.gr.vphi*rr3D*np.sin(th3D), eta=self.gr.radratio,
                         spectral=True, exclude=False)
            ds = np.sin(th3D)*dr + np.cos(th3D)/rr3D*dth
            data = -1./(rr3D*np.sin(th3D))*phideravg(self.gr.vr*np.sin(th3D)+self.gr.vtheta*np.cos(th3D), self.gr.minc)+ds/(rr3D*np.sin(th3D))

            del dr, dth, ds, rr3D, th3D

            data_ic = None

            if labTex:
                label = r'$\omega_z$'
            else:
                label = 'vortz'
        else:
            data, data_ic, label = selectField(self.gr, field, labTex, ic=ic)

        if field in ['entropy', 's', 'S', 'u2', 'b2', 'nrj', 'temperature']:
            normed = False

        # Determine the radius 
        r /= (1-self.gr.radratio) # as we give a normalised radius
        ri = self.gr.radratio/(1.-self.gr.radratio)
        if r < ri and data_ic is not None:
            ind = np.nonzero(np.where(abs(self.gr.radius_ic-r) \
                             == min(abs(self.gr.radius_ic-r)), 1, 0))
            indPlot = ind[0][0]
            rad = self.gr.radius_ic[indPlot]*(1.-self.gr.radratio)
            rprof = data_ic[..., indPlot]
        else:
            ind = np.nonzero(np.where(abs(self.gr.radius-r) \
                             == min(abs(self.gr.radius-r)), 1, 0))
            indPlot = ind[0][0]
            rad = self.gr.radius[indPlot] * (1.-self.gr.radratio)
            rprof = data[..., indPlot]
        # Shifting the azimuth data by lon_shift
        lon_shift = int(lon_shift*self.gr.nphi/360)
        rprof = np.roll(rprof,lon_shift,axis=0)
        rprof = symmetrize(rprof, self.gr.minc)

        fig = radialContour(rprof, rad, label, proj, lon_0, vmax, vmin,
                            lat_0, levels, cm, normed, cbar, tit, lines)

    def equat(self, field='vr', levels=defaultLevels, cm=defaultCm,
              normed=True, vmax=None, vmin=None, cbar=True, tit=True,
              avg=False, normRad=False, ic=False):
        """
        Plot the equatorial cut of a given field

           >>> s = Surf()
           >>> # Equatorial cut of the z-vorticity, 65 contour levels
           >>> s.equat(field='vortz', levels=65, cm='seismic')

           >>> # Minimal plot (no cbar, not title)
           >>> s.equat(field='bphi', tit=False, cbar=False)

           >>> # Control the limit of the colormap from -1e3 to 1e3
           >>> s.equat(field='vr', vmin=-1e3, vmax=1e3, levels=33)

           >>> # Normalise the contour levels radius by radius
           >>> s.equat(field='jphi', normRad=True)

        :param field: the name of the input physical quantity you want to
                      display
        :type field: str
        :param avg: when set to True, an additional figure which shows
                    the radial profile of the input physical quantity
                    (azimuthal average) is also displayed
        :type avg: bool
        :param normRad: when set to True, the contour levels are normalised
                        radius by radius (default is False)
        :type normRad: bool
        :param levels: the number of levels in the contour
        :type levels: int
        :param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
        :type cm: str
        :param tit: display the title of the figure when set to True
        :type tit: bool
        :param cbar: display the colorbar when set to True
        :type cbar: bool
        :param vmax: maximum value of the contour levels
        :type vmax: float
        :param vmin: minimum value of the contour levels
        :type vmin: float
        :param normed: when set to True, the colormap is centered around zero.
                       Default is True, except for entropy/temperature plots.
        :type normed: bool
        :param ic: when set to True, also display the contour levels in
                   the inner core
        :type ic: bool
        """
        phi = np.linspace(0., 2.*np.pi, self.gr.nphi)
        rr, pphi = np.meshgrid(self.gr.radius, phi)
        xx = rr * np.cos(pphi)
        yy = rr * np.sin(pphi)

        if field in ('vortzfluct'):
            philoc = np.linspace(0., 2.*np.pi/self.gr.minc, self.gr.npI)
            rrloc, pphiloc = np.meshgrid(self.gr.radius, philoc)
            vpfluct = self.gr.vphi-self.gr.vphi.mean(axis=0)
            vrfluct = self.gr.vr-self.gr.vr.mean(axis=0)
            dr = rderavg(rrloc*vpfluct[:,self.gr.ntheta//2,:], spectral=False,
                         eta=self.gr.radratio, exclude=True)
            equator = 1./rrloc*(dr-phideravg(vrfluct[:,self.gr.ntheta//2,:], self.gr.minc))
            if labTex:
                label = r"$\omega_z'$"
            else:
                label = 'vortzfluct'
        elif field in ('vortz'):
            philoc = np.linspace(0., 2.*np.pi/self.gr.minc, self.gr.npI)
            rrloc, pphiloc = np.meshgrid(self.gr.radius, philoc)
            dr = rderavg(rrloc*self.gr.vphi[:,self.gr.ntheta//2,:], spectral=False,
                         eta=self.gr.radratio, exclude=True)
            equator = 1./rrloc*(dr - phideravg(self.gr.vr[:,self.gr.ntheta//2,:], self.gr.minc))
            if labTex:
                label = r'$\omega_z$'
            else:
                label = 'vortz'
        elif field in ('jz'):
            philoc = np.linspace(0., 2.*np.pi/self.gr.minc, self.gr.npI)
            rrloc, pphiloc = np.meshgrid(self.gr.radius, philoc)
            dr = rderavg(rrloc*self.gr.Bphi[:,self.gr.ntheta//2,:], spectral=False,
                         eta=self.gr.radratio, exclude=True)
            equator = 1./rrloc*(dr - phideravg(self.gr.Br[:,self.gr.ntheta//2,:], self.gr.minc))
            if labTex:
                label = r'$j_z$'
            else:
                label = 'jz'
        elif field in ('vopot'):
            philoc = np.linspace(0., 2.*np.pi/self.gr.minc, self.gr.npI)
            rrloc, pphiloc = np.meshgrid(self.gr.radius, philoc)
            dr = rderavg(rrloc*self.gr.vphi[:,self.gr.ntheta//2,:], spectral=False,
                         eta=self.gr.radratio, exclude=True)
            wz = 1./rrloc*(dr - phideravg(self.gr.vr[:,self.gr.ntheta//2,:], self.gr.minc))
            temp0, rho0, beta = anelprof(self.gr.radius, self.gr.strat,
                                         self.gr.polind, self.gr.g0, self.gr.g1,
                                         self.gr.g2)
            #equator = (wz + 2./(self.gr.ek))/(rho0)
            height = 2. * np.sqrt( self.gr.radius.max()**2-self.gr.radius**2 )
            equator = (wz + 2./(self.gr.ek))/(rho0*height)
            #equator = wz - 2./(self.gr.ek)*np.log(height)
            label = r'PV'
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            ax1.plot(self.gr.radius, equator.mean(axis=0))
            ax1.plot(self.gr.radius, 2./(self.gr.ek)/(rho0*height))
        elif field in ('rey'):
            temp0, rho0, beta = anelprof(self.gr.radius, self.gr.strat,
                                         self.gr.polind, self.gr.g0, self.gr.g1,
                                         self.gr.g2)
            vp = self.gr.vphi.copy()
            vp = self.gr.vphi- self.gr.vphi.mean(axis=0) # convective vp
            data =  rho0 * self.gr.vr * vp
            if labTex:
                label = r'$\rho v_s v_\phi$'
            else:
                label = r'rho vs vp'
        elif field in ('mr'):
            temp0, rho0, beta = anelprof(self.gr.radius, self.gr.strat,
                                         self.gr.polind, self.gr.g0, self.gr.g1,
                                         self.gr.g2)
            data =  rho0 * self.gr.vr
            if labTex:
                label = r'$\rho v_r$'
            else:
                label = r'rho vr'
        else:
            data, data_ic, label = selectField(self.gr, field, labTex, ic)

        if field not in ('vortz', 'vopot', 'jz', 'vortzfluct'):
            equator = data[:,int(self.gr.ntheta/2),:]
            if ic and data_ic is not None:
                equator_ic = data_ic[:,int(self.gr.ntheta/2),:]

        equator = symmetrize(equator, self.gr.minc)
        if ic and data_ic is not None:
            equator_ic = symmetrize(equator_ic, self.gr.minc)

        if field in ['entropy', 's', 'S', 'u2', 'b2', 'nrj', 'temperature']:
            normed = False

        fig, xx, yy = equatContour( equator, self.gr.radius, self.gr.minc, label,
                                    levels, cm, normed, vmax, vmin, cbar, tit,
                                    normRad)
        ax = fig.get_axes()[0]

        if ic and data_ic is not None:
            phi = np.linspace(0., 2.*np.pi, self.gr.nphi)
            rr, pphi = np.meshgrid(self.gr.radius_ic, phi)
            xx_ic = rr * np.cos(pphi)
            yy_ic = rr * np.sin(pphi)
            if vmax is not None and vmin is not None:
                cs = np.linspace(vmin, vmax, levels)
            else:
                cs = levels
            im_ic = ax.contourf(xx_ic, yy_ic, equator_ic, cs, cmap=cm,
                                extend='both')
            if normed:
                im_ic.set_clim(-max(abs(equator.max()), abs(equator.min())),
                                max(abs(equator.max()), abs(equator.min())))

        # Variable conductivity: add a dashed line
        if hasattr(self.gr, 'nVarCond'):
            if self.gr.nVarCond == 2:
                radi = self.gr.con_RadRatio * self.gr.radius[0]
                ax.plot(radi*np.cos(phi), radi*np.sin(phi), 'k--', lw=1.5)

        # If avg is requested, then display an additional figure
        # with azimutal average
        if avg:
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            ax1.plot(self.gr.radius, equator.mean(axis=0))
            ax1.set_xlabel('Radius')
            ax1.set_ylabel(label)
            ax1.set_xlim(self.gr.radius.min(), self.gr.radius.max())

    def avg(self, field='vphi', levels=defaultLevels, cm=defaultCm,
            normed=True, vmax=None, vmin=None, cbar=True, tit=True,
            pol=False, tor=False, mer=False, merLevels=16, polLevels=16,
            ic=False, lines=False):
        """
        Plot the azimutal average of a given field.

           >>> s = Surf()
           >>> # Axisymmetric zonal flows, 65 contour levels
           >>> s.avg(field='vp', levels=65, cm='seismic')

           >>> # Minimal plot (no cbar, not title)
           >>> s.avg(field='Br', tit=False, cbar=False)

           >>> # Axisymmetric Bphi + poloidal field lines
           >>> s.avg(field='Bp', pol=True, polLevels=8)

           >>> # Omega-effect, contours truncated from -1e3 to 1e3
           >>> s.avg(field='omeffect', vmax=1e3, vmin=-1e3)

        :param field: the field you want to display
        :type field: str
        :param levels: the number of levels in the contourf plot
        :type levels: int
        :param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
        :type cm: str
        :param tit: display the title of the figure when set to True
        :type tit: bool
        :param cbar: display the colorbar when set to True
        :type cbar: bool
        :param vmax: maximum value of the contour levels
        :type vmax: float
        :param vmin: minimum value of the contour levels
        :type vmin: float
        :param normed: when set to True, the colormap is centered around zero.
                       Default is True, except for entropy/temperature plots.
        :type normed: bool
        :param pol: diplay the poloidal field lines contours when set to True
        :type pol: bool
        :param tor: diplay the toroidal axisymmetric field contours when
                    set to True
        :type tor: bool
        :param mer: display the meridional circulation contours when
                    set to True
        :type mer: bool
        :param merLevels: number of contour levels to display meridional
                          circulation
        :type merLevels: int
        :param polLevels: number of contour levels to display poloidal
                          field lines
        :type polLevels: int
        :param ic: when set to True, also display the contour levels in
                   the inner core
        :type ic: bool
        :param lines: when set to True, over-plot solid lines to highlight
                      the limits between two adjacent contour levels
        :type lines: bool
        """
        if pol:
            if ic:
                rr2D = np.zeros((self.gr.ntheta, self.gr.nr+self.gr.n_r_ic_max-1),
                                dtype=self.precision)
                th2D = np.zeros_like(rr2D)
                data = np.zeros_like(rr2D)
                brm = self.gr.Br.mean(axis=0)
                brm_ic = self.gr.Br_ic.mean(axis=0)
                brm = np.concatenate((brm, brm_ic[:, 1:]), axis=-1)
                for i in range(self.gr.ntheta):
                    th2D[i, :] = self.gr.colatitude[i]+np.pi/2.
                for i in range(self.gr.nr):
                    rr2D[:, i] = self.gr.radius[i]
                for i in range(self.gr.n_r_ic_max-1):
                    rr2D[:, i+self.gr.nr] = self.gr.radius_ic[i+1]
                s2D = rr2D * np.abs(np.cos(th2D))
                data[0, :] = -0.5*s2D[0, :]*brm[0, :]*self.gr.colatitude[0]

                for i in range(1, self.gr.ntheta):
                    data[i, :] = data[i-1, :] \
                                - (s2D[i, :]*brm[i, :]+s2D[i-1, :]*brm[i-1, :]) *\
                                  (th2D[i, :]-th2D[i-1, :])
                dataerr = data[-1, :]-0.5*(s2D[-1, :]*brm[-1, :]) *\
                          (np.pi-self.gr.colatitude[-1])
                for i in range(self.gr.ntheta):
                    data[i, :] = data[i, :] - \
                                 dataerr*self.gr.colatitude[i]/np.pi
                poloLines = 0.5*data/np.cos(th2D)
            else:
                rr2D = np.zeros((self.gr.ntheta, self.gr.nr),
                                dtype=self.precision)
                th2D = np.zeros_like(rr2D)
                data = np.zeros_like(rr2D)
                brm = self.gr.Br.mean(axis=0)
                for i in range(self.gr.ntheta):
                    rr2D[i, :] = self.gr.radius
                    th2D[i, :] = self.gr.colatitude[i]+np.pi/2.
                s2D = rr2D * np.abs(np.cos(th2D))
                data[0, :] = -0.5*s2D[0, :]*brm[0, :]*self.gr.colatitude[0]

                for i in range(1, self.gr.ntheta):
                    data[i, :] = data[i-1, :] \
                                - (s2D[i, :]*brm[i, :]+s2D[i-1, :]*brm[i-1, :]) *\
                                  (th2D[i, :]-th2D[i-1, :])
                dataerr = data[-1, :]-0.5*(s2D[-1, :]*brm[-1, :]) *\
                          (np.pi-self.gr.colatitude[-1])
                for i in range(self.gr.ntheta):
                    data[i, :] = data[i, :] - \
                                 dataerr*self.gr.colatitude[i]/np.pi
                poloLines = 0.5*data/np.cos(th2D)

        if mer:
            rr2D = np.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
            th2D = np.zeros_like(rr2D)
            data = np.zeros_like(rr2D)
            if hasattr(self.gr, 'strat'):
                if (self.gr.strat != 0.):
                    temp, rho, beta = anelprof(self.gr.radius, self.gr.strat,
                                               self.gr.polind, g0=self.gr.g0,
                                               g1=self.gr.g1, g2=self.gr.g2)
                else:
                    rho = 1.
            else:
                rho = 1.
            vrm = self.gr.vr.mean(axis=0)*rho
            for i in range(self.gr.ntheta):
                rr2D[i, :] = self.gr.radius
                th2D[i, :] = self.gr.colatitude[i]+np.pi/2.
            s2D = rr2D * np.abs(np.cos(th2D))
            data[0, :] = -0.5*s2D[0, :]*vrm[0, :]*self.gr.colatitude[0]

            for i in range(1, self.gr.ntheta):
                data[i, :] = data[i-1, :] \
                            - (s2D[i, :]*vrm[i, :]+s2D[i-1, :]*vrm[i-1, :]) *\
                              (th2D[i, :]-th2D[i-1, :])
            dataerr = data[-1, :]-0.5*(s2D[-1, :]*vrm[-1, :]) *\
                      (np.pi-self.gr.colatitude[-1])
            for i in range(self.gr.ntheta):
                data[i, :] = data[i, :] - dataerr*self.gr.colatitude[i]/np.pi
            meriLines = 0.5*data/np.cos(th2D)

        if field in ('Vs', 'vs'):
            vr = self.gr.vr
            vt = self.gr.vtheta
            thlin = np.linspace(0., np.pi, self.gr.ntheta)
            th3D = np.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                            dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            data = vr * np.sin(th3D) + vt * np.cos(th3D)
            label = 'Vs'
        elif field in ('entropyreduced'):
            tt = self.gr.entropy.mean(axis=0).mean(axis=0)
            data = self.gr.entropy-tt
            label = 'tt'
        elif field in ('Vz', 'vz'):
            vr = self.gr.vr
            vt = self.gr.vtheta
            thlin = np.linspace(0., np.pi, self.gr.ntheta)
            th3D = np.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                            dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            data = vr * np.cos(th3D) - vt * np.sin(th3D)
            label = 'Vz'
        elif field in ('Omega'):
            if labTex:
                label = r'$\Omega$'
            else:
                label = 'omega'
            th2D = np.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
            rr2D = np.zeros_like(th2D)
            for i in range(self.gr.ntheta):
                th2D[i, :] = self.gr.colatitude[i]
                rr2D[i, :] = self.gr.radius
            s2D = rr2D * np.sin(th2D)
            data = self.gr.vphi/s2D + 1./self.gr.ek
        elif field in ('jphi'):
            if labTex:
                label = r'$j_\phi$'
            else:
                label = 'jphi'
            th2D = np.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
            rr2D = np.zeros_like(th2D)
            for i in range(self.gr.ntheta):
                th2D[i, :] = self.gr.colatitude[i]
                rr2D[i, :] = self.gr.radius
            Brm = self.gr.Br.mean(axis=0)
            Btm = self.gr.Btheta.mean(axis=0)
            data = 1./rr2D*(rderavg(rr2D*Btm, eta=self.gr.radratio) - \
                            thetaderavg(Brm))
        elif field in ('ohm'):
            if labTex:
                label = r'$\lambda\,j^2$'
            else:
                label = 'Ohmic dissipation'
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            dphi = 2.*np.pi/self.gr.nphi

            Op = (np.roll(self.gr.Btheta,-1,axis=2)-np.roll(self.gr.Btheta,1,axis=2))/\
                 (np.roll(rr3D,-1,axis=2)-np.roll(rr3D,1,axis=2)) + \
                 self.gr.Btheta/rr3D - \
                 (np.roll(self.gr.Br,-1,axis=1)-np.roll(self.gr.Br,1,axis=1))/\
                 (rr3D*(np.roll(th3D,-1,axis=1)-np.roll(th3D,1,axis=1)))
            Ot = (np.roll(self.gr.Br,-1,axis=0)-np.roll(self.gr.Br,1,axis=0))/\
                  (rr3D*np.sin(th3D)*dphi) - \
                 (np.roll(self.gr.Bphi,-1,axis=2)-np.roll(self.gr.Bphi,1,axis=2))/\
                 (np.roll(rr3D,-1,axis=2)-np.roll(rr3D,1,axis=2)) - \
                 self.gr.Bphi/rr3D
            Or = (np.roll(self.gr.Bphi,-1,axis=1)-np.roll(self.gr.Bphi,1,axis=1))/\
                 (rr3D*(np.roll(th3D,-1,axis=1)-np.roll(th3D,1,axis=1))) + \
                 np.cos(th3D)*self.gr.Bphi/(rr3D*np.sin(th3D)) - \
                 (np.roll(self.gr.Btheta,-1,axis=0)-np.roll(self.gr.Btheta,1,axis=0))/\
                 (rr3D*np.sin(th3D)*dphi)

            Or[:, 0, :] = Or[:, 1, :]
            Or[:, -1, :] = Or[:, -2, :]
            Ot[..., 0] = Ot[..., 1]
            Ot[..., -1] = Ot[..., -2]
            Op[..., 0] = Op[..., 1]
            Op[..., -1] = Op[..., -2]
            data = (Op**2+Ot**2+Or**2)
            rad = MagicRadial(field='varCond', iplot=False)
            data *= rad.lmbda[::-1] # it starts from ri in MagicRadial
        elif field in ('omeffect'):
            if labTex:
                label = r'$\Omega$-effect'
            else:
                label = r'omega-effect'
            rr2D = np.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
            th2D = np.zeros_like(rr2D)
            for i in range(self.gr.ntheta):
                th2D[i, :] = self.gr.colatitude[i]
                rr2D[i, :] = self.gr.radius
            brm = self.gr.Br.mean(axis=0)
            btm = self.gr.Btheta.mean(axis=0)
            bpm = self.gr.Bphi.mean(axis=0)
            vrm = self.gr.vr.mean(axis=0)
            vtm = self.gr.vtheta.mean(axis=0)
            vpm = self.gr.vphi.mean(axis=0)
            dvpdr = rderavg(vpm, eta=self.gr.radratio, spectral=True,
                            exclude=False)
            dvpdt = thetaderavg(vpm)
            # B. Brown
            # Phi component of <B> dot grad <u>
            #data = brm*dvpdr+btm/rr2D*dvpdt+vrm*bpm/rr2D+\
                   #vtm*bpm*np.cos(th2D)/(np.sin(th2D)*rr2D)
            # M. Schrinner and U. Christensen
            # Phi component of curl <Vphi> x <B>
            data = brm*dvpdr+btm/rr2D*dvpdt-vpm*brm/rr2D-\
                   vpm*btm*np.cos(th2D)/(np.sin(th2D)*rr2D)
        elif field in ('flux'):
            label = 'flux'
            temp0, rho0, beta0 = anelprof(self.gr.radius, self.gr.strat,
                                         self.gr.polind, self.gr.g0, self.gr.g1,
                                         self.gr.g2)
            ssm = self.gr.entropy.mean(axis=0)
            data = rderavg(ssm, self.gr.radratio, spectral=True, exclude=False)
        elif field in ('alphaeffect'):
            if labTex:
                label = r'$-\alpha \langle B_\phi\rangle$'
            else:
                label = 'alpha*Bphi'
            th3D = np.zeros_like(self.gr.vphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            dphi = 2.*np.pi/self.gr.nphi

            vt = self.gr.vtheta - self.gr.vtheta.mean(axis=0)
            vr = self.gr.vr-self.gr.vr.mean(axis=0)
            vp = self.gr.vphi-self.gr.vphi.mean(axis=0)
            wp = (np.roll(vt,-1,axis=2)-np.roll(vt,1,axis=2))/\
                 (np.roll(rr3D,-1,axis=2)-np.roll(rr3D,1,axis=2)) + \
                 vt/rr3D - \
                 (np.roll(vr,-1,axis=1)-np.roll(vr,1,axis=1))/\
                 (rr3D*(np.roll(th3D,-1,axis=1)-np.roll(th3D,1,axis=1)))
            wt = (np.roll(vr,-1,axis=0)-np.roll(vr,1,axis=0))/\
                  (rr3D*np.sin(th3D)*dphi) - \
                 (np.roll(vp,-1,axis=2)-np.roll(vp,1,axis=2))/\
                 (np.roll(rr3D,-1,axis=2)-np.roll(rr3D,1,axis=2)) - \
                 vp/rr3D
            wr = (np.roll(vp,-1,axis=1)-np.roll(vp,1,axis=1))/\
                 (rr3D*(np.roll(th3D,-1,axis=1)-np.roll(th3D,1,axis=1))) + \
                 np.cos(th3D)*vp/(rr3D*np.sin(th3D)) - \
                 (np.roll(vt,-1,axis=0)-np.roll(vt,1,axis=0))/\
                 (rr3D*np.sin(th3D)*dphi)

            wr[:, 0, :] = wr[:, 1, :]
            wr[:, -1, :] = wr[:, -2, :]
            wt[..., 0] = wt[..., 1]
            wt[..., -1] = wt[..., -2]
            wp[..., 0] = wp[..., 1]
            wp[..., -1] = wp[..., -2]
            data = -self.gr.Bphi.mean(axis=0)*(vr*wr+vp*wp+vt*wt)
        elif field in ('emf'):
            if labTex:
                label = r"$\langle u'\times B'\rangle_\phi$"
            else:
                label = 'emf'
            vrp = self.gr.vr-self.gr.vr.mean(axis=0)
            vtp = self.gr.vtheta-self.gr.vtheta.mean(axis=0)
            brp = self.gr.Br-self.gr.Br.mean(axis=0)
            btp = self.gr.Btheta-self.gr.Btheta.mean(axis=0)
            data = vrp*btp-vtp*brp
        elif field in ('hz'):
            if labTex:
                label = r'$H_z$'
            else:
                label = 'Hz'
            vr = self.gr.vr
            vt = self.gr.vtheta
            th3D = np.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            vz = vr * np.cos(th3D) - vt * np.sin(th3D)
            data = self.vortz * vz
            denom = np.sqrt(np.mean(vz**2, axis=0)* np.mean(self.vortz**2, axis=0))
        elif field in ('enstrophy'):
            label = 'Enstrophy'
            normed = False
            data = self.vortz**2
        elif field in ('helicity'):
            label = 'Helicity'
            th3D = np.zeros_like(self.gr.vphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            dphi = 2.*np.pi/self.gr.nphi

            wp = (np.roll(self.gr.vtheta,-1,axis=2)-np.roll(self.gr.vtheta,1,axis=2))/\
                 (np.roll(rr3D,-1,axis=2)-np.roll(rr3D,1,axis=2)) + \
                 self.gr.vtheta/rr3D - \
                 (np.roll(self.gr.vr,-1,axis=1)-np.roll(self.gr.vr,1,axis=1))/\
                 (rr3D*(np.roll(th3D,-1,axis=1)-np.roll(th3D,1,axis=1)))
            wt = (np.roll(self.gr.vr,-1,axis=0)-np.roll(self.gr.vr,1,axis=0))/\
                  (rr3D*np.sin(th3D)*dphi) - \
                 (np.roll(self.gr.vphi,-1,axis=2)-np.roll(self.gr.vphi,1,axis=2))/\
                 (np.roll(rr3D,-1,axis=2)-np.roll(rr3D,1,axis=2)) - \
                 self.gr.vphi/rr3D
            wr = (np.roll(self.gr.vphi,-1,axis=1)-np.roll(self.gr.vphi,1,axis=1))/\
                 (rr3D*(np.roll(th3D,-1,axis=1)-np.roll(th3D,1,axis=1))) + \
                 np.cos(th3D)*self.gr.vphi/(rr3D*np.sin(th3D)) - \
                 (np.roll(self.gr.vtheta,-1,axis=0)-np.roll(self.gr.vtheta,1,axis=0))/\
                 (rr3D*np.sin(th3D)*dphi)

            wr[:, 0, :] = wr[:, 1, :]
            wr[:, -1, :] = wr[:, -2, :]
            wt[..., 0] = wt[..., 1]
            wt[..., -1] = wt[..., -2]
            wp[..., 0] = wp[..., 1]
            wp[..., -1] = wp[..., -2]
            data = self.gr.vr*wr+self.gr.vphi*wp+self.gr.vtheta*wt
            self.hel = data.mean(axis=0)
        elif field in ('poloidal'):
            label = 'poloidal field lines'
            rr2D = np.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
            th2D = np.zeros_like(rr2D)
            data = np.zeros_like(rr2D)
            brm = self.gr.Br.mean(axis=0)
            for i in range(self.gr.ntheta):
                rr2D[i, :] = self.gr.radius
                th2D[i, :] = self.gr.colatitude[i]+np.pi/2.
            s2D = rr2D * np.abs(np.cos(th2D))
            data[0, :] = -0.5*s2D[0, :]*brm[0, :]*self.gr.colatitude[0]

            for i in range(1, self.gr.ntheta):
                data[i, :] = data[i-1, :] \
                            -(s2D[i, :]*brm[i, :]+s2D[i-1,:]*brm[i-1,:])* \
                             (th2D[i,:]-th2D[i-1,:])
            dataerr = data[-1, :]-0.5*(s2D[-1,:]*brm[-1,:])*\
                      (np.pi-self.gr.colatitude[-1])
            for i in range(self.gr.ntheta):
                data[i, :] = data[i, :] - dataerr*self.gr.colatitude[i]/np.pi
            data = 0.5*data/np.cos(th2D)
        elif field in ('meridional'):
            label = "meridional circulation"
            rr2D = np.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
            th2D = np.zeros_like(rr2D)
            data = np.zeros_like(rr2D)
            temp, rho, beta = anelprof(self.gr.radius, self.gr.strat,
                           self.gr.polind,
                           g0=self.gr.g0, g1=self.gr.g1, g2=self.gr.g2)
            vrm = self.gr.vr.mean(axis=0)*rho
            for i in range(self.gr.ntheta):
                rr2D[i, :] = self.gr.radius
                th2D[i, :] = self.gr.colatitude[i]+np.pi/2.
            s2D = rr2D * np.abs(np.cos(th2D))
            data[0, :] = -0.5*s2D[0, :]*vrm[0, :]*self.gr.colatitude[0]

            for i in range(1, self.gr.ntheta):
                data[i, :] = data[i-1, :] \
                            -(s2D[i, :]*vrm[i, :]+s2D[i-1,:]*vrm[i-1,:])* \
                             (th2D[i,:]-th2D[i-1,:])
            dataerr = data[-1, :]-0.5*(s2D[-1,:]*vrm[-1,:])*\
                      (np.pi-self.gr.colatitude[-1])
            for i in range(self.gr.ntheta):
                data[i, :] = data[i, :] - dataerr*self.gr.colatitude[i]/np.pi
            data = 0.5*data/np.cos(th2D)
        elif field in ('ra', 'ratio'):
            label = 'Ratio'
            data = self.gr.vphi**2#/(self.gr.vphi**2+\
                   #self.gr.vtheta**2+self.gr.vr**2)
            denom = np.mean(self.gr.vphi**2+ self.gr.vtheta**2+self.gr.vr**2,
                           axis=0)
            #denom = 1.
        elif field in ('beta'):
            if labTex:
                label = r'$\beta$'
            else:
                label = r'd ln rho/dr'
            temp0, rho0, beta = anelprof(self.gr.radius, self.gr.strat,
                                         self.gr.polind, self.gr.g0, self.gr.g1,
                                         self.gr.g2)
            data = beta * np.ones_like(self.gr.vr)#* self.gr.vr
        elif field in ('angular'):
            label = 'Angular momentum'
            th2D = np.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
            rr2D = np.zeros_like(th2D)
            rho2D = np.zeros_like(th2D)
            if hasattr(self.gr, 'strat'):
                temp0, rho0, beta = anelprof(self.gr.radius, self.gr.strat,
                                         self.gr.polind, self.gr.g0, self.gr.g1,
                                         self.gr.g2)
            else:
                rho0 = 1.
            for i in range(self.gr.ntheta):
                rho2D[i, :] = rho0
                rr2D[i, :] = self.gr.radius
                th2D[i, :] = self.gr.colatitude[i]
            s2D = rr2D * np.sin(th2D)
            if self.gr.ek > 0:                                    # Outer boundary rotating
                norm = self.gr.radius[0]**2/self.gr.ek
                data = (self.gr.vphi*s2D+1./self.gr.ek*s2D**2)/norm
            else:                                                # Outer boundary non-rotating, ek = -1
                norm = self.gr.omega_ic1*self.gr.radius[-1]**2
                data = (self.gr.vphi*s2D)/norm

        elif field in ('Cr', 'cr'):
            if hasattr(self.gr, 'strat'):
                temp0, rho0, beta = anelprof(self.gr.radius, self.gr.strat,
                                         self.gr.polind, self.gr.g0, self.gr.g1,
                                         self.gr.g2)
            else:
                rho0 = 1.
            vr = self.gr.vr
            vt = self.gr.vtheta
            vp = self.gr.vphi.copy()
            vp = self.gr.vphi- self.gr.vphi.mean(axis=0) # convective vp
            thlin = np.linspace(0., np.pi, self.gr.ntheta)
            th3D = np.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            vs = vr * np.sin(th3D) + vt * np.cos(th3D)
            data =  vs * vp
            denom = np.sqrt(np.mean(vs**2, axis=0)* np.mean(vp**2, axis=0))
            if labTex:
                label = r'$\langle v_s v_\phi\rangle$'
            else:
                label = 'vs vphi'
        elif field == 'vortz':
            data = self.vortz
            label = 'vortz'
        elif field == 'vopot':
            temp0, rho0, beta = anelprof(self.gr.radius, self.gr.strat,
                                         self.gr.polind, self.gr.g0, self.gr.g1,
                                         self.gr.g2)
            height = 2. * np.sqrt( self.gr.radius.max()**2-self.gr.radius**2 )
            data = (self.vortz+2./self.gr.ek)/(rho0*height)
            label = 'Pot. vort.'
        elif field in ('rhocr'):
            temp0, rho0, beta = anelprof(self.gr.radius, self.gr.strat,
                                         self.gr.polind, self.gr.g0, self.gr.g1,
                                         self.gr.g2)
            vr = self.gr.vr
            vt = self.gr.vtheta
            vp = self.gr.vphi.copy()
            vp = self.gr.vphi- self.gr.vphi.mean(axis=0) # convective vp
            thlin = np.linspace(0., np.pi, self.gr.ntheta)
            th3D = np.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            vs = vr * np.sin(th3D) + vt * np.cos(th3D)
            data =  rho0 * vs * vp
            denom = np.sqrt(np.mean(rho0*vs**2, axis=0)* np.mean(rho0*vp**2, axis=0))
            if labTex:
                label = r'$C_{s\phi}$'
            else:
                label = 'Csp'
        elif field in ('Cz', 'cz'):
            vr = self.gr.vr
            vt = self.gr.vtheta
            vp = self.gr.vphi.copy()
            vp = self.gr.vphi- self.gr.vphi.mean(axis=0) # convective vp
            thlin = np.linspace(0., np.pi, self.gr.ntheta)
            th3D = np.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            vz = vr * np.cos(th3D) - vt * np.sin(th3D)
            data =  vz * vp
            denom = np.sqrt(np.mean(vz**2, axis=0)* np.mean(vp**2, axis=0))
            if labTex:
                label = r'$\langle v_z v_\phi\rangle$'
            else:
                label = 'vz vphi'
        elif field in ('dvzdz'):
            if labTex:
                label = r'$\partial u_z/\partial z$'
            else:
                label = 'dvz/dz'
            vr = self.gr.vr
            vt = self.gr.vtheta
            thlin = np.linspace(0., np.pi, self.gr.ntheta)
            th3D = np.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            data = (vr * np.cos(th3D) - vt * np.sin(th3D))
        else:
            data, data_ic, label = selectField(self.gr, field, labTex,
                                               ic=ic)

        if field not in ('Cr', 'cr', 'ra', 'ratio', 'Cz', 'cz', 'hz', 'jphi',
                         'rhocr', 'omeffect', 'poloidal', 'flux', 'meridional'):
            phiavg = data.mean(axis=0)
            if ic and data_ic is not None: 
                phiavg_ic = data_ic.mean(axis=0)
        elif field == 'balance':
            phiavg = zderavg(data.mean(axis=0), eta=self.gr.radratio,
                             spectral=True, exclude=True)
            phiavg = phiavg + data1.mean(axis=0)
        elif field == 'dvzdz':
            phiavg = zderavg(data.mean(axis=0), eta=self.gr.radratio,
                             spectral=True, exclude=True)
        elif field in ('omeffect', 'poloidal', 'flux', 'meridional', 'jphi'):
            phiavg = data
        else:
            ro = self.gr.radius[0]
            ri = self.gr.radius[-1]
            fac = 2./(np.pi*(ro**2-ri**2))
            facOTC = ro**2.*(np.pi-2.*np.arcsin(self.gr.radratio))/2. \
                     -ri**2*np.sqrt(1.-self.gr.radratio**2)/self.gr.radratio
            facOTC = 1./facOTC
            facITC = ri**2*np.sqrt(1.-self.gr.radratio**2)/self.gr.radratio \
                     +(ro**2-ri**2)* np.arcsin(self.gr.radratio) \
                     -ri**2/2.*(np.pi - 2.*np.arcsin(self.gr.radratio))
            facITC = 1./facITC
            #mask = np.where(denom == 0, 1, 0)
            phiavg = data.mean(axis=0)

            TC = np.array([], dtype=self.precision)
            outTC = np.array([], dtype=self.precision)
            inTC = np.array([], dtype=self.precision)
            denomTC = np.array([], dtype=self.precision)
            denomoutTC = np.array([], dtype=self.precision)
            denominTC = np.array([], dtype=self.precision)
            integ = np.array([], dtype=self.precision)
            for k, th in enumerate(self.gr.colatitude):
                rr = self.gr.radius[::-1]
                dat = phiavg[k, ::-1] * rr
                dat2 = denom[k, ::-1] * rr
                corr = intcheb(dat, self.gr.nr-1, ri, ro)
                TC = np.append(TC, corr)
                corr2 = intcheb(dat2, self.gr.nr-1, ri, ro)
                denomTC = np.append(denomTC, corr2)
                if th >= np.arcsin(self.gr.radratio)  and \
                   th <= np.pi - np.arcsin(self.gr.radratio):
                    # Outside tangent cylinder
                    val = trapz(dat[rr >= ri/np.sin(th)], rr[rr >= ri/np.sin(th)])
                    outTC = np.append(outTC, val)
                    integ = np.append(integ, th)
                    val2 = trapz(dat2[rr >= ri/np.sin(th)], rr[rr >= ri/np.sin(th)])
                    denomoutTC = np.append(denomoutTC, val2)
                    # Inside tangent cylinder
                    val = trapz(dat[rr < ri/np.sin(th)], rr[rr < ri/np.sin(th)])
                    inTC = np.append(inTC, val)
                    val2 = trapz(dat2[rr < ri/np.sin(th)], rr[rr < ri/np.sin(th)])
                    denominTC = np.append(denominTC, val2)
                else:
                    val= intcheb(dat, self.gr.nr-1, ri, ro)
                    inTC = np.append(inTC, val)
                    val2= intcheb(dat2, self.gr.nr-1, ri, ro)
                    denominTC = np.append(denominTC, val2)

            num = fac*trapz(TC, self.gr.colatitude)
            den = fac*trapz(denomTC, self.gr.colatitude)
            print('Correlation', num/den)
            num = facOTC*trapz(outTC, integ)
            den = facOTC*trapz(denomoutTC, integ)
            print('Correlation out TC', num/den)
            num = facITC*trapz(inTC, self.gr.colatitude)
            den = facITC*trapz(denominTC, self.gr.colatitude)
            print('Correlation in TC', num/den)

            mask = np.where(denom == 0, 1, 0)
            phiavg /= (denom + mask)
            #phiavg /= den

        if field in ['entropy', 's', 'S', 'u2', 'b2', 'nrj', 'temperature']:
            normed = False

        fig, xx, yy, im = merContour(phiavg, self.gr.radius, label, levels, cm,
                                     normed, vmax, vmin, cbar, tit, lines=lines)
        ax = fig.get_axes()[0]

        if ic:
            th = np.linspace(0., np.pi, phiavg.shape[0])
            ri = self.gr.radratio/(1.-self.gr.radratio)
            rr, tth = np.meshgrid(self.gr.radius_ic, th)
            xx_ic = rr * np.sin(tth)
            yy_ic = rr * np.cos(tth)

            if data_ic is not None:
                if vmax is not None and vmin is not None:
                    cs = np.linspace(vmin, vmax, levels)
                else:
                    cs = levels
                im_ic = ax.contourf(xx_ic, yy_ic, phiavg_ic, cs, cmap=cm,
                                    extend='both')
                if normed:
                    im_ic.set_clim(-max(abs(phiavg.max()), abs(phiavg.min())),
                                    max(abs(phiavg.max()), abs(phiavg.min())))

            ax.plot([0, 0], [-ri, ri], 'k-')
            
        if pol:
            if ic:
                xx_big = np.concatenate((xx, xx_ic[:, 1:]), axis=-1)
                yy_big = np.concatenate((yy, yy_ic[:, 1:]), axis=-1)
                ax.contour(xx_big, yy_big, poloLines, polLevels, colors=['k'],
                           linewidths=[0.8, 0.8])
            else:
                ax.contour(xx, yy, poloLines, polLevels, colors=['k'],
                           linewidths=[0.8, 0.8])
        elif tor:
            toroLines = self.gr.Bphi.mean(axis=0)
            ax.contour(xx, yy, toroLines, polLevels, colors=['k'],
                       linewidths=[0.8])
        elif mer:
            maxMeri = abs(meriLines).max()
            minMeri = -maxMeri
            lev = np.linspace(minMeri, maxMeri, merLevels)
            im2 = ax.contour(xx, yy, meriLines, lev, colors=['k'],
                       linewidths=[0.8])

        # Variable conductivity: add a dashed line
        if hasattr(self.gr, 'nVarCond'):
            if self.gr.nVarCond == 2:
                radi = self.gr.con_RadRatio * self.gr.radius[0]
                th = np.linspace(0, np.pi, self.gr.ntheta)
                ax.plot(radi*np.sin(th), radi*np.cos(th), 'k--')

    def slice(self, field='Bphi', lon_0=0., levels=defaultLevels, cm=defaultCm,
              normed=True, vmin=None, vmax=None, cbar=True, tit=True,
              grid=False, nGridLevs=16, normRad=False, ic=False):
        """
        Plot an azimuthal slice of a given field.

           >>> s = Surf()
           >>> # vphi at 0, 30, 60 degrees in longitude
           >>> s.slice(field='vp', lon_0=[0, 30, 60], levels=65, cm='seismic')

           >>> # Minimal plot (no cbar, not title)
           >>> s.avg(field='vp', lon_0=32, tit=False, cbar=False)

           >>> # Axisymmetric Bphi + poloidal field lines
           >>> s.avg(field='Bp', pol=True, polLevels=8)

           >>> # Omega-effect, contours truncated from -1e3 to 1e3
           >>> s.avg(field='omeffect', vmax=1e3, vmin=-1e3)


        :param field: the field you want to display
        :type field: str
        :param lon_0: the longitude of the slice in degrees, or a list of longitudes
        :type lon_0: float or list
        :param levels: the number of levels in the contourf plot
        :type levels: int
        :param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
        :type cm: str
        :param tit: display the title of the figure when set to True
        :type tit: bool
        :param cbar: display the colorbar when set to True
        :type cbar: bool
        :param vmax: maximum value of the contour levels
        :type vmax: float
        :param vmin: minimum value of the contour levels
        :type vmin: float
        :param grid: display or hide the grid
        :type grid: bool
        :param nGridLevs: number of grid levels
        :type nGridLevs: int
        :param normRad: when set to True, the contour levels are normalised
                        radius by radius (default is False)
        :type normRad: bool
        :param ic: when set to True, also display the contour levels in
                   the inner core
        :type ic: bool
        """
        if field in ('Vs', 'vs'):
            if labTex:
                label = r'$v_s$'
            else:
                label = 'vs'
            vr = self.gr.vr
            vt = self.gr.vtheta
            thlin = np.linspace(0., np.pi, self.gr.ntheta)
            th3D = np.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            data = vr * np.sin(th3D) + vt * np.cos(th3D)
        elif field in ('Vz', 'vz'):
            if labTex:
                label = '$v_z$'
            else:
                label = 'vz'
            vr = self.gr.vr
            vt = self.gr.vtheta
            thlin = np.linspace(0., np.pi, self.gr.ntheta)
            th3D = np.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            data = vr * np.cos(th3D) - vt * np.sin(th3D)
        elif field in ('anel'):
            if labTex:
                label = r'$\beta v_r$'
            else:
                label = r'beta vr'
            temp0, rho0, beta = anelprof(self.gr.radius, self.gr.strat,
                                         self.gr.polind, self.gr.g0, self.gr.g1,
                                         self.gr.g2)
            data = beta * self.gr.vr
        elif field in ('dvzdz'):
            if labTex:
                label = '$\partial v_z / \partial z$'
            else:
                label = 'dvz/dz'
            vr = self.gr.vr
            vt = self.gr.vtheta
            thlin = np.linspace(0., np.pi, self.gr.ntheta)
            th3D = np.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            data = vr * np.cos(th3D) - vt * np.sin(th3D)
        elif field in ('ohm'):
            label = 'Ohmic dissipation/1e6'
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            dphi = 2.*np.pi/self.gr.nphi

            Op = (np.roll(self.gr.Btheta,-1,axis=2)-np.roll(self.gr.Btheta,1,axis=2))/\
                 (np.roll(rr3D,-1,axis=2)-np.roll(rr3D,1,axis=2)) + \
                 self.gr.Btheta/rr3D - \
                 (np.roll(self.gr.Br,-1,axis=1)-np.roll(self.gr.Br,1,axis=1))/\
                 (rr3D*(np.roll(th3D,-1,axis=1)-np.roll(th3D,1,axis=1)))
            Ot = (np.roll(self.gr.Br,-1,axis=0)-np.roll(self.gr.Br,1,axis=0))/\
                  (rr3D*np.sin(th3D)*dphi) - \
                 (np.roll(self.gr.Bphi,-1,axis=2)-np.roll(self.gr.Bphi,1,axis=2))/\
                 (np.roll(rr3D,-1,axis=2)-np.roll(rr3D,1,axis=2)) - \
                 self.gr.Bphi/rr3D
            Or = (np.roll(self.gr.Bphi,-1,axis=1)-np.roll(self.gr.Bphi,1,axis=1))/\
                 (rr3D*(np.roll(th3D,-1,axis=1)-np.roll(th3D,1,axis=1))) + \
                 np.cos(th3D)*self.gr.Bphi/(rr3D*np.sin(th3D)) - \
                 (np.roll(self.gr.Btheta,-1,axis=0)-np.roll(self.gr.Btheta,1,axis=0))/\
                 (rr3D*np.sin(th3D)*dphi)

            Or[:, 0, :] = Or[:, 1, :]
            Or[:, -1, :] = Or[:, -2, :]
            Ot[..., 0] = Ot[..., 1]
            Ot[..., -1] = Ot[..., -2]
            Op[..., 0] = Op[..., 1]
            Op[..., -1] = Op[..., -2]
            data = (Op**2+Ot**2+Or**2)/1e6
        elif field in ('flux'):
            data = rderavg(self.gr.entropy, eta=self.gr.radratio)
            label = 'flux'
        elif field in ('mag_pres_force_r'):
            data = -rderavg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, eta=self.gr.radratio)/2.0
            label = 'Rad. mag. pres. force'
        elif field in ('mag_pres_force_t'):
            rr3D = np.zeros_like(self.gr.Bphi)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            data = -thetaderavg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, order=2)/rr3D/2.0
            label = 'Lati. mag. pres. force'
        elif field in ('mag_pres_force_p'):
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = -phideravg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, self.gr.minc)/(rr3D*np.sin(th3D))/2.0
            label = 'Longi. mag. pres. force'
        elif field in ('mag_tens_force_r'):
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = self.gr.Br * rderavg(self.gr.Br, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Br, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Br, self.gr.minc) / np.sin(th3D) / rr3D - \
                   (self.gr.Btheta**2 + self.gr.Bphi**2) / rr3D
            label = 'Rad. tens. force'
        elif field in ('mag_tens_force_t'):
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = self.gr.Br * rderavg(self.gr.Btheta, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Btheta, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Btheta, self.gr.minc) / np.sin(th3D) / rr3D + \
                   self.gr.Btheta * self.gr.Br / rr3D - \
                   self.gr.Bphi**2 * np.arctan(th3D) / rr3D
            label = 'Lati. tens. force'
        elif field in ('mag_tens_force_p'):
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = self.gr.Br * rderavg(self.gr.Bphi, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Bphi, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Bphi, self.gr.minc) / np.sin(th3D) / rr3D+ \
                   self.gr.Bphi * self.gr.Br / rr3D + \
                   self.gr.Bphi * self.gr.Btheta * np.arctan(th3D) / rr3D
            label = 'Longi. tens. force'
        elif field in ('Lorentz_r'):
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = -rderavg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, eta=self.gr.radratio)/2.0 + \
                   self.gr.Br * rderavg(self.gr.Br, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Br, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Br, self.gr.minc) / np.sin(th3D) / rr3D - \
                   (self.gr.Btheta**2 + self.gr.Bphi**2) / rr3D
            label = 'Radial Lorentz force'
        elif field in ('Lorentz_t'):
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = -thetaderavg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, order=2)/rr3D/2.0 + \
                   self.gr.Br * rderavg(self.gr.Btheta, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Btheta, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Btheta, self.gr.minc) / np.sin(th3D) / rr3D + \
                   self.gr.Btheta * self.gr.Br / rr3D - \
                   self.gr.Bphi**2 * np.arctan(th3D) / rr3D
            label = 'Lati. Lorentz force'
        elif field in ('Lorentz_p'):
            th3D = np.zeros_like(self.gr.Bphi)
            rr3D = np.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = -phideravg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, self.gr.minc)/(rr3D*np.sin(th3D))/2.0 + \
                   self.gr.Br * rderavg(self.gr.Bphi, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Bphi, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Bphi, self.gr.minc) / np.sin(th3D) / rr3D+ \
                   self.gr.Bphi * self.gr.Br / rr3D + \
                   self.gr.Bphi * self.gr.Btheta * np.arctan(th3D) / rr3D
            label = 'Longi. Lorentz force'
        else:
            data, data_ic, label = selectField(self.gr, field, labTex, ic)

        data = symmetrize(data, self.gr.minc)
        if ic and data_ic is not None:
            data_ic = symmetrize(data_ic, self.gr.minc)

        th = np.linspace(np.pi/2, -np.pi/2, self.gr.ntheta)
        rr, tth = np.meshgrid(self.gr.radius, th)
        xx = rr * np.cos(tth)
        yy = rr * np.sin(tth)
        phi = np.linspace(0., 360, self.gr.nphi)
        if ic:
            ri = self.gr.radratio/(1.-self.gr.radratio)
            rr_ic, tth_ic = np.meshgrid(self.gr.radius_ic, th)
            xx_ic = rr_ic * np.cos(tth_ic)
            yy_ic = rr_ic * np.sin(tth_ic)

        lon_0 = np.asarray(lon_0)
        cmap = plt.get_cmap(cm)

        if len(lon_0) > 1:
            fig = plt.figure(figsize=(2.5*len(lon_0), 5.1))
            fig.subplots_adjust(top=0.99, bottom=0.01, right=0.99, left=0.01,
                              wspace=0.01)
            for k, lon in enumerate(lon_0):
                ind = np.nonzero(np.where(abs(phi-lon) \
                                == min(abs(phi-lon)), 1, 0))
                indPlot = ind[0][0]
                phislice = data[indPlot, ...]
                if ic and data_ic is not None:
                    phislice_ic = data_ic[indPlot, ...]
                if field == 'dvzdz':
                    phislice = zderavg(phislice, eta=self.gr.radratio,
                                       spectral=True, exclude=True)
                elif field == 'balance':
                    phislice = zderavg(phislice, eta=self.gr.radratio,
                                       spectral=True, exclude=True)
                    phislice1 = data1[indPlot, ...]
                    phislice = phislice + phislice1

                if normRad: # Normalise each radius
                    maxS = np.sqrt(np.mean(phislice**2, axis=0))
                    phislice[:, maxS!=0.] /= maxS[maxS!=0.]

                ax = fig.add_subplot(1,len(lon_0),k+1, frameon=False)
                if vmax is not None or vmin is not None:
                    normed = False
                    cs = np.linspace(vmin, vmax, levels)
                    im = ax.contourf(xx, yy, phislice, cs, cmap=cmap,
                                     extend='both')
                else:
                    cs = levels
                    im = ax.contourf(xx, yy, phislice, cs, cmap=cmap)
                ax.plot(self.gr.radius[0]*np.cos(th), self.gr.radius[0]*np.sin(th),
                   'k-')
                ax.plot(self.gr.radius[-1]*np.cos(th),
                        self.gr.radius[-1]*np.sin(th), 'k-')
                ax.plot([0., 0], [self.gr.radius[-1], self.gr.radius[0]], 'k-')
                ax.plot([0., 0], [-self.gr.radius[-1], -self.gr.radius[0]],
                        'k-')

                if ic and data_ic is not None:
                    if vmax is not None and vmin is not None:
                        cs = np.linspace(vmin, vmax, levels)
                    else:
                        cs = levels
                    im_ic = ax.contourf(xx_ic, yy_ic, phislice_ic, cs,
                                        cmap=cmap, extend='both')
                    if normed:
                        im_ic.set_clim(-max(abs(phislice.max()),
                                            abs(phislice.min())),
                                        max(abs(phislice.max()),
                                            abs(phislice.min())))
            
                    ax.plot([0, 0], [-ri, ri], 'k-') 

                ax.axis('off')

                tit1 = r'${}^\circ$'.format(lon)
                ax.text(0.9, 0.9, tit1, fontsize=18,
                      horizontalalignment='right',
                      verticalalignment='center',
                      transform = ax.transAxes)
                #fig.colorbar(im)
                if field not in ['entropy', 's', 'S'] and normed is True:
                    im.set_clim(-max(abs(phislice.max()), abs(phislice.min())),
                                 max(abs(phislice.max()), abs(phislice.min())))

                #To avoid white lines on pdfs

                for c in im.collections:
                    c.set_edgecolor("face")


        else:
            ind = np.nonzero(np.where(abs(phi-lon_0[0]) \
                            == min(abs(phi-lon_0[0])), 1, 0))
            indPlot = ind[0][0]
            phislice = data[indPlot, ...]
            if ic and data_ic is not None:
                phislice_ic = data_ic[indPlot, ...]
            if field == 'dvzdz':
                phislice = zderavg(phislice, eta=self.gr.radratio,
                                   spectral=True, exclude=True)
            elif field == 'balance':
                phislice = zderavg(phislice, eta=self.gr.radratio,
                                   spectral=True, exclude=True)
                phislice1 = data1[indPlot, ...]
                phislice = phislice + phislice1
            elif field == 'vs':
                th2D = np.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
                rr2D = np.zeros_like(th2D)
                for i in range(self.gr.ntheta):
                    th2D[i, :] = self.gr.colatitude[i]
                    rr2D[i, :] = self.gr.radius
                s2D = rr2D * np.sin(th2D)
                ro = 1./(1.-self.gr.radratio)
                coeff = -s2D/(ro**2-s2D**2)
                phislice *= coeff


            if tit:
                if cbar:
                    fig = plt.figure(figsize=(5,7.5))
                    ax = fig.add_axes([0.01, 0.01, 0.69, 0.91])
                else:
                    fig = plt.figure(figsize=(3.5,7.5))
                    ax = fig.add_axes([0.01, 0.01, 0.98, 0.91])
                ax.set_title(label, fontsize=24)
            else:
                if cbar:
                    fig = plt.figure(figsize=(5,7))
                    ax = fig.add_axes([0.01, 0.01, 0.69, 0.98])
                else:
                    fig = plt.figure(figsize=(3.5,7))
                    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])

            if vmax is not None or vmin is not None:
                normed = False
                cs = np.linspace(vmin, vmax, levels)
                im = ax.contourf(xx, yy, phislice, cs, cmap=cmap, extend='both')
            else:
                cs = levels
                im = ax.contourf(xx, yy, phislice, cs, cmap=cmap)
            ax.plot(self.gr.radius[0]*np.cos(th), self.gr.radius[0]*np.sin(th),
                   'k-', lw=1.5)
            ax.plot(self.gr.radius[-1]*np.cos(th), self.gr.radius[-1]*np.sin(th),
                   'k-', lw=1.5)
            ax.plot([0., 0], [self.gr.radius[-1], self.gr.radius[0]], 'k-', lw=1.5)
            ax.plot([0., 0], [-self.gr.radius[-1], -self.gr.radius[0]], 'k-', lw=1.5)

            if ic and data_ic is not None:
                if vmax is not None and vmin is not None:
                    cs = np.linspace(vmin, vmax, levels)
                else:
                    cs = levels
                im_ic = ax.contourf(xx_ic, yy_ic, phislice_ic, cs,
                                    cmap=cmap, extend='both')
                if normed:
                    im_ic.set_clim(-max(abs(phislice.max()),
                                        abs(phislice.min())),
                                    max(abs(phislice.max()),
                                        abs(phislice.min())))
        
                ax.plot([0, 0], [-ri, ri], 'k-') 

            if hasattr(self.gr, 'epsS'):
                if self.gr.epsS != 0:
                    rad = MagicRadial(field='anel', iplot=False)
                    idx = np.nonzero(np.where(abs(rad.dsdr)==abs(rad.dsdr).min(),
                                1, 0))[0][0]
                    ax.plot(self.gr.radius[idx]*np.cos(th),
                            self.gr.radius[idx]*np.sin(th), 'k--', lw=2)

            if grid:
                ax.contour(xx, yy, tth, nGridLevs, colors='k', linestyles='--',
                           linewidths=0.5)
            ax.axis('off')

            # Add the colorbar at the right place
            pos = ax.get_position()
            l, b, w, h = pos.bounds
            if cbar:
                if tit:
                    cax = fig.add_axes([0.82, 0.46-0.7*h/2., 0.04, 0.7*h])
                else:
                    cax = fig.add_axes([0.82, 0.5-0.7*h/2., 0.04, 0.7*h])
                mir = fig.colorbar(im, cax=cax)

            # Normalise the data
            if field not in ['entropy', 's', 'S'] and normed is True:
                im.set_clim(-max(abs(phislice.max()), abs(phislice.min())),
                             max(abs(phislice.max()), abs(phislice.min())))

            #To avoid white lines on pdfs

            for c in im.collections:
                c.set_edgecolor("face")



def report(nvar=1, levels=defaultLevels, lclean=True):
    """
    This subroutine prepares a pdf document that gather some important diagnostics

    :param lclean: clean or not the LaTeX files
    :type lclean: bool
    :param levels: number of contour levels
    :param levels: int
    :param nvar: number of graphic files
    :param nvar: int
    """
    file = open('report.tex', 'w')
    file.write("\documentclass[a4paper,10pt]{article}\n")
    file.write("\\usepackage[utf8]{inputenc}\n")
    file.write("\\usepackage{amsmath,amsfonts,amssymb}\n")
    file.write("\\usepackage[francais]{babel}\n")
    file.write("\\usepackage[T1]{fontenc}\n")
    file.write("\\usepackage[dvips]{graphicx}\n")
    file.write("\\usepackage{geometry}\n")
    file.write("\\usepackage[pdftex]{xcolor}\n")

    file.write("\geometry{hmargin=1cm,vmargin=2cm}\n")

    file.write("\\begin{document}\n")

    s = Surf(ivar=nvar)

    st = "Ek = {:.2e}, Ra = {:.2e}, Pr = {:.1f}, $N_{\\rho}$={:.2f}, $\eta$={.1f}".format(
            s.gr.ek, s.gr.ra, s.gr.pr, s.gr.strat, s.gr.radratio)
    file.write("\\begin{center}\\begin{large}\n")
    file.write(" "+st+"\n")
    file.write("\\end{large}\\end{center}\n")

    r1 = 0.98
    r3 = 1.03 * s.gr.radratio
    r2 = (r1+r3)/2.

    s.avg(field='vp', levels=levels, cm='RdYlBu_r', normed=True)
    plt.savefig('vp.png')
    plt.close()

    s.avg(field='entropy', levels=levels, cm='RdYlBu_r', normed=True)
    plt.savefig('entropy.png')
    plt.close()

    s.equat(field='entropy', levels=levels, cm='RdYlBu_r', normed=False)
    plt.savefig('equ_s.png')
    plt.close()

    s.equat(field='vr', levels=levels, cm='RdYlBu_r', normed=False)
    plt.savefig('equ_vr.png')
    plt.close()

    s.surf(field='vp', cm='RdYlBu_r', levels=levels, r=r1, proj='moll',
           normed=False)
    plt.savefig('surf_vp.png')
    plt.close()

    s.surf(field='vp', cm='RdYlBu', levels=levels, r=r2, proj='moll',
           normed=False)
    plt.savefig('surf_vp_08.png')
    plt.close()

    s.surf(field='vp', cm='RdYlBu', levels=levels, r=r3, proj='moll',
           normed=False)
    plt.savefig('surf_vp_06.png')
    plt.close()

    s.surf(field='vr', cm='RdYlBu', levels=levels, r=r1, proj='moll',
           normed=False)
    plt.savefig('surf_vr.png')
    plt.close()

    s.surf(field='vr', cm='RdYlBu', levels=levels, r=r2, proj='moll',
           normed=False)
    plt.savefig('surf_vr_08.png')
    plt.close()

    s.surf(field='vr', cm='RdYlBu', levels=levels, r=r3, proj='moll',
           normed=False)
    plt.savefig('surf_vr_06.png')
    plt.close()

    file.write("\\begin{figure}[htbp]\n")
    file.write(" \\centering\n")
    file.write(" \\includegraphics[width=9cm]{equ_s}\n")
    file.write(" \\includegraphics[width=9cm]{equ_vr}\n")
    file.write(" \\includegraphics[height=9cm]{entropy}\n")
    file.write(" \\includegraphics[height=9cm]{vp}\n")
    file.write("\\end{figure}\n")
    file.write("\\newpage\n")

    file.write("\\begin{figure}\n")
    file.write(" \\centering\n")
    file.write(" \\includegraphics[width=18cm]{surf_vr_06}\n")
    file.write(" \\includegraphics[width=18cm]{surf_vr_08}\n")
    file.write(" \\includegraphics[width=18cm]{surf_vr}\n")
    file.write("\\end{figure}\n")
    file.write("\\newpage\n")
    file.write("\\begin{figure}\n")
    file.write(" \\includegraphics[width=18cm]{surf_vp_06}\n")
    file.write(" \\includegraphics[width=18cm]{surf_vp_08}\n")
    file.write(" \\includegraphics[width=18cm]{surf_vp}\n")
    file.write("\\end{figure}\n")

    file.write("\\end{document}")


    file.close()
    os.system("pdflatex report.tex")
    if lclean:
        os.system("rm vp.png entropy.png equ_s.png equ_vr.png surf_vp.png \
                   surf_vr.png surf_vr_06.png surf_vr_08.png surf_vp_06.png\
                   surf_vp_08.png")
        os.system("rm report.log report.aux report.tex")
