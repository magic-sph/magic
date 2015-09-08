# -*- coding: utf-8 -*-
from magic import MagicGraph, MagicSetup, MagicRadial
from magic.setup import labTex
from .libmagic import *
import matplotlib.pyplot as P
import os
import numpy as N
from scipy.integrate import trapz

__author__  = "$Author$"
__date__   = "$Date$"
__version__ = "$Revision$"



class Surf:

    def __init__(self, ivar=None, datadir='.', vort=False, ave=False, tag=None,
                 precision='Float32'):
        """
        Read the graphic file, and possibly calculate vorticity if requested

        :param ivar: number of the graphic file
        :param ave: ave=True for a time-averaged graphic file
        :param tag: TAG extension of the graphic file
        :param vort: a boolean to specify whether, one wants to compute the 3-D
                     vorticiy components (take care of the memory imprint)
        :param datadir: the working directory
        :param precision: the storage precision of the graphic file (single or
                          double precision)
        """
        self.precision = precision
        self.datadir = datadir
        self.gr = MagicGraph(ivar=ivar, datadir=self.datadir, ave=ave, tag=tag,
                             precision=self.precision)

        if vort:
            thlin = self.gr.colatitude
            th3D = N.zeros_like(self.gr.vphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]

            s3D = rr3D * N.sin(th3D)
            dtheta = thetaderavg(self.gr.vphi*s3D)
            dr = rderavg(self.gr.vphi*s3D, eta=self.gr.radratio, spectral=True, 
                         exclude=False)
            ds = N.sin(th3D)*dr + N.cos(th3D)/rr3D*dtheta
            vs = self.gr.vr * N.sin(th3D) + self.gr.vtheta * N.cos(th3D) # 'vs'
            self.vortz = -1./s3D*phideravg(vs, self.gr.minc)+ds/s3D
            del dr, dtheta, ds, rr3D, th3D, s3D

    def surf(self, field='Bphi', proj='hammer', lon_0=0., r=0.85, vmax=None, 
             vmin=None, lat_0=30., levels=16, cm='RdYlBu_r', 
             normed=True, cbar=True, tit=True, lines=False):
        """
        Plot the surface distribution of a field at a given
        normalised radius.

        :param field: the name of the field you want to display
        :param proj: the type of projection. Default is Hammer, in case
                     you want to use 'ortho' or 'moll', then Basemap is
                     required.
        :param r: the radius you want to display (in normalised values)
        :param levels: the number of levels in the contour
        :param cm: the colormap name
        :param tit: display/hide the title of the figure
        :param cbar: display/hide the colorbar
        """
        r /= (1-self.gr.radratio) # as we give a normalised radius
        ind = N.nonzero(N.where(abs(self.gr.radius-r) \
                        == min(abs(self.gr.radius-r)), 1, 0))
        indPlot = ind[0][0]
        rad = self.gr.radius[indPlot] * (1.-self.gr.radratio)

        if proj != 'ortho':
            lon_0 = 0.

        if field in ('Vs', 'vs'):
            vr = self.gr.vr
            vt = self.gr.vtheta
            thlin = N.linspace(0., N.pi, self.gr.ntheta)
            th3D = N.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            data = vr * N.sin(th3D) + vt * N.cos(th3D)
            if labTex:
                label = r'$v_s$'
            else:
                label = r'vs'
        elif field in ('Vz', 'vz'):
            vr = self.gr.vr
            vt = self.gr.vtheta
            thlin = N.linspace(0., N.pi, self.gr.ntheta)
            th3D = N.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            data = vr * N.cos(th3D) - vt * N.sin(th3D)
            if labTex:
                label = r'$v_z$'
            else:
                label = r'vz'
        elif field in ('thu'):
            data = self.gr.vr*(self.gr.entropy-self.gr.entropy.mean(axis=0))
            label = 'thu'
        elif field in ('flux'):
            data = rderavg(self.gr.entropy, eta=self.gr.radratio)
            label = 'flux'
        elif field in ('mag_pres_force_r'):
            data = -rderavg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, eta=self.gr.radratio)/2.0
            label = 'Rad. mag. pres. force'
        elif field in ('mag_pres_force_t'):
            rr3D = N.zeros_like(self.gr.Bphi)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            data = -thetaderavg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, order=2)/rr3D/2.0
            label = 'Lati. mag. pres. force'
        elif field in ('mag_pres_force_p'):
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = -phideravg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, self.gr.minc)/(rr3D*N.sin(th3D))/2.0
            label = 'Longi. mag. pres. force'
        elif field in ('mag_tens_force_r'):
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = self.gr.Br * rderavg(self.gr.Br, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Br, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Br, self.gr.minc) / N.sin(th3D) / rr3D - \
                   (self.gr.Btheta**2 + self.gr.Bphi**2) / rr3D
            label = 'Rad. tens. force'
        elif field in ('mag_tens_force_t'):
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = self.gr.Br * rderavg(self.gr.Btheta, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Btheta, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Btheta, self.gr.minc) / N.sin(th3D) / rr3D + \
                   self.gr.Btheta * self.gr.Br / rr3D - \
                   self.gr.Bphi**2 * N.arctan(th3D) / rr3D
            label = 'Lati. tens. force'
        elif field in ('mag_tens_force_p'):
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = self.gr.Br * rderavg(self.gr.Bphi, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Bphi, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Bphi, self.gr.minc) / N.sin(th3D) / rr3D+ \
                   self.gr.Bphi * self.gr.Br / rr3D + \
                   self.gr.Bphi * self.gr.Btheta * N.arctan(th3D) / rr3D
            label = 'Longi. tens. force'
        elif field in ('Lorentz_r'):
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = -rderavg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, eta=self.gr.radratio)/2.0 + \
                   self.gr.Br * rderavg(self.gr.Br, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Br, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Br, self.gr.minc) / N.sin(th3D) / rr3D - \
                   (self.gr.Btheta**2 + self.gr.Bphi**2) / rr3D
            label = 'Radial Loretz force'
        elif field in ('Lorentz_t'):
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = -thetaderavg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, order=2)/rr3D/2.0 + \
                   self.gr.Br * rderavg(self.gr.Btheta, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Btheta, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Btheta, self.gr.minc) / N.sin(th3D) / rr3D + \
                   self.gr.Btheta * self.gr.Br / rr3D - \
                   self.gr.Bphi**2 * N.arctan(th3D) / rr3D
            label = 'Lati. Loretz force'
        elif field in ('Lorentz_p'):
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = -phideravg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, self.gr.minc)/(rr3D*N.sin(th3D))/2.0 + \
                   self.gr.Br * rderavg(self.gr.Bphi, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Bphi, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Bphi, self.gr.minc) / N.sin(th3D) / rr3D+ \
                   self.gr.Bphi * self.gr.Br / rr3D + \
                   self.gr.Bphi * self.gr.Btheta * N.arctan(th3D) / rr3D
            label = 'Longi. Loretz force'
        elif field in ('ohm'):
            label = 'Ohmic dissipation/1e6'
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            dphi = 2.*N.pi/self.gr.nphi

            Op = (N.roll(self.gr.Btheta,-1,axis=2)-N.roll(self.gr.Btheta,1,axis=2))/\
                 (N.roll(rr3D,-1,axis=2)-N.roll(rr3D,1,axis=2)) + \
                 self.gr.Btheta/rr3D - \
                 (N.roll(self.gr.Br,-1,axis=1)-N.roll(self.gr.Br,1,axis=1))/\
                 (rr3D*(N.roll(th3D,-1,axis=1)-N.roll(th3D,1,axis=1)))
            Ot = (N.roll(self.gr.Br,-1,axis=0)-N.roll(self.gr.Br,1,axis=0))/\
                  (rr3D*N.sin(th3D)*dphi) - \
                 (N.roll(self.gr.Bphi,-1,axis=2)-N.roll(self.gr.Bphi,1,axis=2))/\
                 (N.roll(rr3D,-1,axis=2)-N.roll(rr3D,1,axis=2)) - \
                 self.gr.Bphi/rr3D
            Or = (N.roll(self.gr.Bphi,-1,axis=1)-N.roll(self.gr.Bphi,1,axis=1))/\
                 (rr3D*(N.roll(th3D,-1,axis=1)-N.roll(th3D,1,axis=1))) + \
                 N.cos(th3D)*self.gr.Bphi/(rr3D*N.sin(th3D)) - \
                 (N.roll(self.gr.Btheta,-1,axis=0)-N.roll(self.gr.Btheta,1,axis=0))/\
                 (rr3D*N.sin(th3D)*dphi)

            Or[:, 0, :] = Or[:, 1, :] 
            Or[:, -1, :] = Or[:, -2, :]
            Ot[..., 0] = Ot[..., 1]
            Ot[..., -1] = Ot[..., -2]
            Op[..., 0] = Op[..., 1]
            Op[..., -1] = Op[..., -2]
            data = (Op**2+Ot**2+Or**2)/1e6
        elif field in ('vortzfluct'):
            th3D = N.zeros_like(self.gr.vphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            s3D = rr3D*N.sin(th3D)
            dth = thetaderavg((self.gr.vphi-self.gr.vphi.mean(axis=0))*rr3D*N.sin(th3D))
            dr = rderavg((self.gr.vphi-self.gr.vphi.mean(axis=0))*rr3D*N.sin(th3D), \
                         eta=self.gr.radratio, spectral=True, exclude=False)
            ds = N.sin(th3D)*dr + N.cos(th3D)/rr3D*dth
            data = -1./(rr3D*N.sin(th3D))*phideravg(self.gr.vr*N.sin(th3D)+self.gr.vtheta*N.cos(th3D), self.gr.minc)+ds/(rr3D*N.sin(th3D))

            del dr, dth, ds, rr3D, th3D

            if labTex:
                label = r"$\omega_z'$"
            else:
                label = 'vortzfluct'
        elif field in ('vortz'):
            th3D = N.zeros_like(self.gr.vphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            s3D = rr3D*N.sin(th3D)
            dth = thetaderavg(self.gr.vphi*rr3D*N.sin(th3D))
            dr = rderavg(self.gr.vphi*rr3D*N.sin(th3D), eta=self.gr.radratio, 
                         spectral=True, exclude=False)
            ds = N.sin(th3D)*dr + N.cos(th3D)/rr3D*dth
            data = -1./(rr3D*N.sin(th3D))*phideravg(self.gr.vr*N.sin(th3D)+self.gr.vtheta*N.cos(th3D), self.gr.minc)+ds/(rr3D*N.sin(th3D))

            del dr, dth, ds, rr3D, th3D

            if labTex:
                label = r'$\omega_z$'
            else:
                label = 'vortz'
        else:
            data, label = selectField(self.gr, field, labTex)

        phi = N.linspace(-N.pi, N.pi, self.gr.nphi)
        theta = N.linspace(N.pi/2, -N.pi/2, self.gr.ntheta)
        pphi, ttheta = N.mgrid[-N.pi:N.pi:self.gr.nphi*1j,
                            N.pi/2.:-N.pi/2.:self.gr.ntheta*1j]
        lon2 = pphi * 180./N.pi
        lat2 = ttheta * 180./N.pi

        delat = 30. ; delon = 60.
        circles = N.arange(delat, 90.+delat, delat).tolist()+\
                  N.arange(-delat, -90.-delat, -delat).tolist()
        meridians = N.arange(-180+delon, 180, delon)

        if proj == 'moll' or proj == 'hammer':
            if tit:
                if cbar:
                    fig = P.figure(figsize=(9,4.5))
                    ax = fig.add_axes([0.01, 0.01, 0.87, 0.87])
                else:
                    fig = P.figure(figsize=(8,4.5))
                    ax = fig.add_axes([0.01, 0.01, 0.98, 0.87])
                ax.set_title('%s: r/ro = %.3f' % (label, rad), 
                             fontsize=24)
            else:
                if cbar:
                    fig = P.figure(figsize=(9,4))
                    ax = fig.add_axes([0.01, 0.01, 0.87, 0.98])
                else:
                    fig = P.figure(figsize=(8,4))
                    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])
                tit1 = r'%.2f Ro' % rad
                ax.text(0.12, 0.9, tit1, fontsize=16,
                      horizontalalignment='right',
                      verticalalignment='center',
                      transform = ax.transAxes)
        else:
            if tit:
                if cbar:
                    fig = P.figure(figsize=(6,5.5))
                    ax = fig.add_axes([0.01, 0.01, 0.82, 0.9])
                else:
                    fig = P.figure(figsize=(5,5.5))
                    ax = fig.add_axes([0.01, 0.01, 0.98, 0.9])
                ax.set_title('%s: r/ro = %.3f' % (label, rad), 
                             fontsize=24)
            else:
                if cbar:
                    fig = P.figure(figsize=(6,5))
                    ax = fig.add_axes([0.01, 0.01, 0.82, 0.98])
                else:
                    fig = P.figure(figsize=(5,5))
                    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])
                tit1 = r'%.2f Ro' % rad
                ax.text(0.12, 0.9, tit1, fontsize=16,
                      horizontalalignment='right',
                      verticalalignment='center',
                      transform = ax.transAxes)
            
        if proj != 'hammer':
            from mpl_toolkits.basemap import Basemap
            map = Basemap(projection=proj, lon_0=lon_0, lat_0=lat_0,
                          resolution='c')
            map.drawparallels([0.], dashes=[2, 3], linewidth=0.5)
            map.drawparallels(circles, dashes=[2,3], linewidth=0.5)
            map.drawmeridians(meridians, dashes=[2,3], linewidth=0.5)
            map.drawmeridians([-180], dashes=[20,0], linewidth=0.5)
            map.drawmeridians([180], dashes=[20,0], linewidth=0.5)
            x, y = list(map(lon2, lat2))
        else:
            x, y = hammer2cart(ttheta, pphi)
            for lat0 in circles:
                x0, y0 = hammer2cart(lat0*N.pi/180., phi)
                ax.plot(x0, y0, 'k:', linewidth=0.7)
            for lon0 in meridians:
                x0, y0 = hammer2cart(theta, lon0*N.pi/180.)
                ax.plot(x0, y0, 'k:', linewidth=0.7)
            xxout, yyout  = hammer2cart(theta, -N.pi)
            xxin, yyin  = hammer2cart(theta, N.pi)
            ax.plot(xxin, yyin, 'k-')
            ax.plot(xxout, yyout, 'k-')
            ax.axis('off')

        rprof = data[..., indPlot]
        rprof = symmetrize(rprof, self.gr.minc)

        cmap = P.get_cmap(cm)

        if proj == 'ortho': 
            lats = N.linspace(-90., 90., self.gr.ntheta)
            dat = map.transform_scalar(rprof.T, phi*180/N.pi, lats, 
                                       self.gr.nphi, self.gr.ntheta, masked=True)
            im = map.imshow(dat, cmap=cmap)
        else:
            if vmax is not None or vmin is not None:
                normed = False
                cs = N.linspace(vmin, vmax, levels)
                im = ax.contourf(x, y, rprof, cs, cmap=cmap, extend='both')
                if lines:
                    ax.contour(x, y, rprof, cs, colors='k', linewidths=0.5, extend='both')
                    ax.contour(x, y, rprof, 1, colors=['k'])
                #im = ax.pcolormesh(x, y, rprof, cmap=cmap, antialiased=True)
            else:
                cs = levels
                im = ax.contourf(x, y, rprof, cs, cmap=cmap)
                if lines:
                    ax.contour(x, y, rprof, cs, colors='k', linewidths=0.5)
                    ax.contour(x, y, rprof, 1, colors=['k'])
                #im = ax.pcolormesh(x, y, rprof, cmap=cmap, antialiased=True)

        # Add the colorbar at the right place
        pos = ax.get_position()
        l, b, w, h = pos.bounds
        if cbar:
            if tit:
                cax = fig.add_axes([0.9, 0.46-0.7*h/2., 0.03, 0.7*h])
            else:
                cax = fig.add_axes([0.9, 0.51-0.7*h/2., 0.03, 0.7*h])
            mir = fig.colorbar(im, cax=cax)

        # Normalise around zero
        if field not in ['entropy', 's', 'S', 'u2', 'b2', 'nrj'] \
            and normed is True:
            im.set_clim(-max(abs(rprof.max()), abs(rprof.min())), 
                         max(abs(rprof.max()), abs(rprof.min())))

    def equat(self, field='vr', levels=16, cm='RdYlBu_r', normed=True,
              vmax=None, vmin=None, cbar=True, tit=True, avg=False, normRad=False):
        """
        Plot the equatorial cut of a given field

        :param field: the name of the field
        :param levels: the number of contour levels
        :param cm: the name of the colormap
        :param cbar: display/hide the colorbar
        :param tit: display/hide the title
        :param avg: display also the radial profile of the quantity
                    (average in azimuth)
        :pram normRad: normalise each radius
        """
        phi = N.linspace(0., 2.*N.pi, self.gr.nphi)
        rr, pphi = N.meshgrid(self.gr.radius, phi)
        xx = rr * N.cos(pphi)
        yy = rr * N.sin(pphi)

        if field in ('vortzfluct'):
            philoc = N.linspace(0., 2.*N.pi/self.gr.minc, self.gr.npI)
            rrloc, pphiloc = N.meshgrid(self.gr.radius, philoc)
            vpfluct = self.gr.vphi-self.gr.vphi.mean(axis=0)
            vrfluct = self.gr.vr-self.gr.vr.mean(axis=0)
            dr = rderavg(rrloc*vpfluct[:,self.gr.ntheta/2,:], spectral=False,
                         eta=self.gr.radratio, exclude=True)
            equator = 1./rrloc*(dr-phideravg(vrfluct[:,self.gr.ntheta/2,:], self.gr.minc))
            if labTex:
                label = r"$\omega_z'$"
            else:
                label = 'vortzfluct'
        elif field in ('vortz'):
            philoc = N.linspace(0., 2.*N.pi/self.gr.minc, self.gr.npI)
            rrloc, pphiloc = N.meshgrid(self.gr.radius, philoc)
            dr = rderavg(rrloc*self.gr.vphi[:,self.gr.ntheta/2,:], spectral=False,
                         eta=self.gr.radratio, exclude=True)
            equator = 1./rrloc*(dr - phideravg(self.gr.vr[:,self.gr.ntheta/2,:], self.gr.minc))
            if labTex:
                label = r'$\omega_z$'
            else:
                label = 'vortz'
        elif field in ('jz'):
            philoc = N.linspace(0., 2.*N.pi/self.gr.minc, self.gr.npI)
            rrloc, pphiloc = N.meshgrid(self.gr.radius, philoc)
            dr = rderavg(rrloc*self.gr.Bphi[:,self.gr.ntheta/2,:], spectral=False,
                         eta=self.gr.radratio, exclude=True)
            equator = 1./rrloc*(dr - phideravg(self.gr.Br[:,self.gr.ntheta/2,:], self.gr.minc))
            if labTex:
                label = r'$j_z$'
            else:
                label = 'jz'
        elif field in ('vopot'):
            philoc = N.linspace(0., 2.*N.pi/self.gr.minc, self.gr.npI)
            rrloc, pphiloc = N.meshgrid(self.gr.radius, philoc)
            dr = rderavg(rrloc*self.gr.vphi[:,self.gr.ntheta/2,:], spectral=False,
                         eta=self.gr.radratio, exclude=True)
            wz = 1./rrloc*(dr - phideravg(self.gr.vr[:,self.gr.ntheta/2,:], self.gr.minc))
            temp0, rho0, beta = anelprof(self.gr.radius, self.gr.strat,
                                         self.gr.polind, self.gr.g0, self.gr.g1,
                                         self.gr.g2)
            #equator = (wz + 2./(self.gr.ek))/(rho0)
            height = 2. * N.sqrt( self.gr.radius.max()**2-self.gr.radius**2 )
            equator = (wz + 2./(self.gr.ek))/(rho0*height)
            #equator = wz - 2./(self.gr.ek)*N.log(height)
            label = r'PV'
            fig1 = P.figure()
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
            data, label = selectField(self.gr, field, labTex)

        if field not in ('vortz', 'vopot', 'jz', 'vortzfluct'):
            equator = data[:, self.gr.ntheta/2,:]

        equator = symmetrize(equator, self.gr.minc)

        if normRad: # Normalise each radius
            maxS = N.sqrt(N.mean(equator**2, axis=0))
            equator[:, maxS!=0.] /= maxS[maxS!=0.]

        if tit:
            if cbar:
                fig = P.figure(figsize=(6.5,5.5))
                ax = fig.add_axes([0.01, 0.01, 0.76, 0.9])
            else:
                fig = P.figure(figsize=(5,5.5))
                ax = fig.add_axes([0.01, 0.01, 0.98, 0.9])
            ax.set_title(label, fontsize=24)
        else:
            if cbar:
                fig = P.figure(figsize=(6.5,5))
                ax = fig.add_axes([0.01, 0.01, 0.76, 0.98])
            else:
                fig = P.figure(figsize=(5, 5))
                ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])

        cmap = P.get_cmap(cm)
        if vmax is not None or vmin is not None:
            normed = False
            cs = N.linspace(vmin, vmax, levels)
            im = ax.contourf(xx, yy, equator, cs, cmap=cmap, extend='both')
        else:
            cs = levels
            im = ax.contourf(xx, yy, equator, cs, cmap=cmap)
            #im = ax.pcolormesh(xx, yy, equator, cmap=cmap, antialiased=True)
        ax.plot(self.gr.radius[0] * N.cos(phi), self.gr.radius[0]*N.sin(phi),
               'k-', lw=1.5)
        ax.plot(self.gr.radius[-1] * N.cos(phi), self.gr.radius[-1]*N.sin(phi),
               'k-', lw=1.5)

        ax.axis('off')

        # Variable conductivity: add a dashed line
        if hasattr(self.gr, 'nVarCond'):
            if self.gr.nVarCond == 2:
                radi = self.gr.con_RadRatio * self.gr.radius[0]
                ax.plot(radi*N.cos(phi), radi*N.sin(phi), 'k--', lw=1.5)

        #if hasattr(self.gr, 'epsS'):
        #    rad = MagicRadial(field='anel', iplot=False)
        #    idx = N.nonzero(N.where(abs(rad.dsdr)==abs(rad.dsdr).min(), 1, 0))[0][0]
        #    ax.plot(self.gr.radius[idx]*N.cos(phi), self.gr.radius[idx]*N.sin(phi), 
        #            'k--', lw=2)

        # Add the colorbar at the right place
        pos = ax.get_position()
        l, b, w, h = pos.bounds
        if cbar:
            if tit:
                cax = fig.add_axes([0.85, 0.46-0.7*h/2., 0.03, 0.7*h])
            else:
                cax = fig.add_axes([0.85, 0.5-0.7*h/2., 0.03, 0.7*h])
            mir = fig.colorbar(im, cax=cax)

        # Normalise data 
        if field not in ['entropy', 's', 'S', 'u2', 'b2', 'nrj'] \
            and normed is True:
            im.set_clim(-max(abs(equator.max()), abs(equator.min())), 
                         max(abs(equator.max()), abs(equator.min())))

        # If avg is requested, then display an additional figure
        # with azimutal average 
        if avg:
            fig1 = P.figure()
            ax1 = fig1.add_subplot(111)
            ax1.plot(self.gr.radius, equator.mean(axis=0))
            ax1.set_xlabel('Radius')
            ax1.set_ylabel(label)
            ax1.set_xlim(self.gr.radius.min(), self.gr.radius.max())

    def avg(self, field='vphi', levels=16, cm='RdYlBu_r', normed=True,
            vmax=None, vmin=None, cbar=True, tit=True,
            pol=False, tor=False, mer=False,
            merLevels=16, polLevels=16):
        """
        Plot the azimutal average of a given field.

        :param field: the field you want to display
        :param levels: the number of levels of the colormap
        :param cm: the name of the colormap
        :param cbar: display/hide the colorbar
        :param tit: display/hide the title
        :param pol: poloidal field lines
        :param mer: meridional circulation
        :param merLevels: number of levels to display meridional circulation
        :param polLevels: number of levels to display poloidal field lines
        """
        if pol:
            rr2D = N.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
            th2D = N.zeros_like(rr2D)
            data = N.zeros_like(rr2D)
            brm = self.gr.Br.mean(axis=0)
            for i in range(self.gr.ntheta):
                rr2D[i, :] = self.gr.radius
                th2D[i, :] = self.gr.colatitude[i]+N.pi/2.
            s2D = rr2D * N.abs(N.cos(th2D))
            data[0, :] = -0.5*s2D[0, :]*brm[0, :]*self.gr.colatitude[0]

            for i in range(1, self.gr.ntheta):
                data[i, :] = data[i-1, :] \
                            -(s2D[i, :]*brm[i, :]+s2D[i-1,:]*brm[i-1,:])* \
                             (th2D[i,:]-th2D[i-1,:])
            dataerr = data[-1, :]-0.5*(s2D[-1,:]*brm[-1,:])*\
                      (N.pi-self.gr.colatitude[-1])
            for i in range(self.gr.ntheta):
                data[i, :] = data[i, :] - dataerr*self.gr.colatitude[i]/N.pi
            poloLines = 0.5*data/N.cos(th2D)

        if mer:
            rr2D = N.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
            th2D = N.zeros_like(rr2D)
            data = N.zeros_like(rr2D)
            if hasattr(self.gr, 'strat'):
                temp, rho, beta = anelprof(self.gr.radius, self.gr.strat, 
                           self.gr.polind,
                           g0=self.gr.g0, g1=self.gr.g1, g2=self.gr.g2)
            else:
                rho = 1.
            vrm = self.gr.vr.mean(axis=0)*rho
            for i in range(self.gr.ntheta):
                rr2D[i, :] = self.gr.radius
                th2D[i, :] = self.gr.colatitude[i]+N.pi/2.
            s2D = rr2D * N.abs(N.cos(th2D))
            data[0, :] = -0.5*s2D[0, :]*vrm[0, :]*self.gr.colatitude[0]

            for i in range(1, self.gr.ntheta):
                data[i, :] = data[i-1, :] \
                            -(s2D[i, :]*vrm[i, :]+s2D[i-1,:]*vrm[i-1,:])* \
                             (th2D[i,:]-th2D[i-1,:])
            dataerr = data[-1, :]-0.5*(s2D[-1,:]*vrm[-1,:])*\
                      (N.pi-self.gr.colatitude[-1])
            for i in range(self.gr.ntheta):
                data[i, :] = data[i, :] - dataerr*self.gr.colatitude[i]/N.pi
            meriLines = 0.5*data/N.cos(th2D)

        if field in ('Vs', 'vs'):
            vr = self.gr.vr
            vt = self.gr.vtheta
            thlin = N.linspace(0., N.pi, self.gr.ntheta)
            th3D = N.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr), 
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            data = vr * N.sin(th3D) + vt * N.cos(th3D)
            label = 'Vs'
        elif field in ('Vz', 'vz'):
            vr = self.gr.vr
            vt = self.gr.vtheta
            thlin = N.linspace(0., N.pi, self.gr.ntheta)
            th3D = N.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr), 
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            data = vr * N.cos(th3D) - vt * N.sin(th3D)
            label = 'Vz'
        elif field in ('Omega'):
            if labTex:
                label = r'$\Omega$'
            else:
                label = 'omega'
            th2D = N.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
            rr2D = N.zeros_like(th2D)
            for i in range(self.gr.ntheta):
                th2D[i, :] = self.gr.colatitude[i]
                rr2D[i, :] = self.gr.radius
            s2D = rr2D * N.sin(th2D)
            data = self.gr.vphi/s2D + 1./self.gr.ek
        elif field in ('jphi'):
            if labTex:
                label = r'$j_\phi$'
            else:
                label = 'jphi'
            th2D = N.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
            rr2D = N.zeros_like(th2D)
            for i in range(self.gr.ntheta):
                th2D[i, :] = self.gr.colatitude[i]
                rr2D[i, :] = self.gr.radius
            Brm = self.gr.Br.mean(axis=0)
            Btm = self.gr.Btheta.mean(axis=0)
            data = 1./rr2D*( rderavg(rr2D*Btm, eta=self.gr.radratio) - \
                             thetaderavg(Brm) )
        elif field in ('ohm'):
            if labTex:
                label = r'$\lambda\,j^2$'
            else:
                label = 'Ohmic dissipation'
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            dphi = 2.*N.pi/self.gr.nphi

            Op = (N.roll(self.gr.Btheta,-1,axis=2)-N.roll(self.gr.Btheta,1,axis=2))/\
                 (N.roll(rr3D,-1,axis=2)-N.roll(rr3D,1,axis=2)) + \
                 self.gr.Btheta/rr3D - \
                 (N.roll(self.gr.Br,-1,axis=1)-N.roll(self.gr.Br,1,axis=1))/\
                 (rr3D*(N.roll(th3D,-1,axis=1)-N.roll(th3D,1,axis=1)))
            Ot = (N.roll(self.gr.Br,-1,axis=0)-N.roll(self.gr.Br,1,axis=0))/\
                  (rr3D*N.sin(th3D)*dphi) - \
                 (N.roll(self.gr.Bphi,-1,axis=2)-N.roll(self.gr.Bphi,1,axis=2))/\
                 (N.roll(rr3D,-1,axis=2)-N.roll(rr3D,1,axis=2)) - \
                 self.gr.Bphi/rr3D
            Or = (N.roll(self.gr.Bphi,-1,axis=1)-N.roll(self.gr.Bphi,1,axis=1))/\
                 (rr3D*(N.roll(th3D,-1,axis=1)-N.roll(th3D,1,axis=1))) + \
                 N.cos(th3D)*self.gr.Bphi/(rr3D*N.sin(th3D)) - \
                 (N.roll(self.gr.Btheta,-1,axis=0)-N.roll(self.gr.Btheta,1,axis=0))/\
                 (rr3D*N.sin(th3D)*dphi)

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
            rr2D = N.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
            th2D = N.zeros_like(rr2D)
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
                   #vtm*bpm*N.cos(th2D)/(N.sin(th2D)*rr2D)
            # M. Schrinner and U. Christensen
            # Phi component of curl <Vphi> x <B> 
            data = brm*dvpdr+btm/rr2D*dvpdt-vpm*brm/rr2D-\
                   vpm*btm*N.cos(th2D)/(N.sin(th2D)*rr2D)
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
            th3D = N.zeros_like(self.gr.vphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            dphi = 2.*N.pi/self.gr.nphi

            vt = self.gr.vtheta - self.gr.vtheta.mean(axis=0)
            vr = self.gr.vr-self.gr.vr.mean(axis=0)
            vp = self.gr.vphi-self.gr.vphi.mean(axis=0)
            wp = (N.roll(vt,-1,axis=2)-N.roll(vt,1,axis=2))/\
                 (N.roll(rr3D,-1,axis=2)-N.roll(rr3D,1,axis=2)) + \
                 vt/rr3D - \
                 (N.roll(vr,-1,axis=1)-N.roll(vr,1,axis=1))/\
                 (rr3D*(N.roll(th3D,-1,axis=1)-N.roll(th3D,1,axis=1)))
            wt = (N.roll(vr,-1,axis=0)-N.roll(vr,1,axis=0))/\
                  (rr3D*N.sin(th3D)*dphi) - \
                 (N.roll(vp,-1,axis=2)-N.roll(vp,1,axis=2))/\
                 (N.roll(rr3D,-1,axis=2)-N.roll(rr3D,1,axis=2)) - \
                 vp/rr3D
            wr = (N.roll(vp,-1,axis=1)-N.roll(vp,1,axis=1))/\
                 (rr3D*(N.roll(th3D,-1,axis=1)-N.roll(th3D,1,axis=1))) + \
                 N.cos(th3D)*vp/(rr3D*N.sin(th3D)) - \
                 (N.roll(vt,-1,axis=0)-N.roll(vt,1,axis=0))/\
                 (rr3D*N.sin(th3D)*dphi)

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
            th3D = N.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            vz = vr * N.cos(th3D) - vt * N.sin(th3D)
            data = self.vortz * vz
            denom = N.sqrt(N.mean(vz**2, axis=0)* N.mean(self.vortz**2, axis=0))
        elif field in ('enstrophy'):
            label = 'Enstrophy'
            normed = False
            data = self.vortz**2
        elif field in ('helicity'):
            label = 'Helicity'
            th3D = N.zeros_like(self.gr.vphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            dphi = 2.*N.pi/self.gr.nphi

            wp = (N.roll(self.gr.vtheta,-1,axis=2)-N.roll(self.gr.vtheta,1,axis=2))/\
                 (N.roll(rr3D,-1,axis=2)-N.roll(rr3D,1,axis=2)) + \
                 self.gr.vtheta/rr3D - \
                 (N.roll(self.gr.vr,-1,axis=1)-N.roll(self.gr.vr,1,axis=1))/\
                 (rr3D*(N.roll(th3D,-1,axis=1)-N.roll(th3D,1,axis=1)))
            wt = (N.roll(self.gr.vr,-1,axis=0)-N.roll(self.gr.vr,1,axis=0))/\
                  (rr3D*N.sin(th3D)*dphi) - \
                 (N.roll(self.gr.vphi,-1,axis=2)-N.roll(self.gr.vphi,1,axis=2))/\
                 (N.roll(rr3D,-1,axis=2)-N.roll(rr3D,1,axis=2)) - \
                 self.gr.vphi/rr3D
            wr = (N.roll(self.gr.vphi,-1,axis=1)-N.roll(self.gr.vphi,1,axis=1))/\
                 (rr3D*(N.roll(th3D,-1,axis=1)-N.roll(th3D,1,axis=1))) + \
                 N.cos(th3D)*self.gr.vphi/(rr3D*N.sin(th3D)) - \
                 (N.roll(self.gr.vtheta,-1,axis=0)-N.roll(self.gr.vtheta,1,axis=0))/\
                 (rr3D*N.sin(th3D)*dphi)

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
            rr2D = N.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
            th2D = N.zeros_like(rr2D)
            data = N.zeros_like(rr2D)
            brm = self.gr.Br.mean(axis=0)
            for i in range(self.gr.ntheta):
                rr2D[i, :] = self.gr.radius
                th2D[i, :] = self.gr.colatitude[i]+N.pi/2.
            s2D = rr2D * N.abs(N.cos(th2D))
            data[0, :] = -0.5*s2D[0, :]*brm[0, :]*self.gr.colatitude[0]

            for i in range(1, self.gr.ntheta):
                data[i, :] = data[i-1, :] \
                            -(s2D[i, :]*brm[i, :]+s2D[i-1,:]*brm[i-1,:])* \
                             (th2D[i,:]-th2D[i-1,:])
            dataerr = data[-1, :]-0.5*(s2D[-1,:]*brm[-1,:])*\
                      (N.pi-self.gr.colatitude[-1])
            for i in range(self.gr.ntheta):
                data[i, :] = data[i, :] - dataerr*self.gr.colatitude[i]/N.pi
            data = 0.5*data/N.cos(th2D)
        elif field in ('meridional'):
            label = "meridional circulation"
            rr2D = N.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
            th2D = N.zeros_like(rr2D)
            data = N.zeros_like(rr2D)
            temp, rho, beta = anelprof(self.gr.radius, self.gr.strat, 
                           self.gr.polind,
                           g0=self.gr.g0, g1=self.gr.g1, g2=self.gr.g2)
            vrm = self.gr.vr.mean(axis=0)*rho
            for i in range(self.gr.ntheta):
                rr2D[i, :] = self.gr.radius
                th2D[i, :] = self.gr.colatitude[i]+N.pi/2.
            s2D = rr2D * N.abs(N.cos(th2D))
            data[0, :] = -0.5*s2D[0, :]*vrm[0, :]*self.gr.colatitude[0]

            for i in range(1, self.gr.ntheta):
                data[i, :] = data[i-1, :] \
                            -(s2D[i, :]*vrm[i, :]+s2D[i-1,:]*vrm[i-1,:])* \
                             (th2D[i,:]-th2D[i-1,:])
            dataerr = data[-1, :]-0.5*(s2D[-1,:]*vrm[-1,:])*\
                      (N.pi-self.gr.colatitude[-1])
            for i in range(self.gr.ntheta):
                data[i, :] = data[i, :] - dataerr*self.gr.colatitude[i]/N.pi
            data = 0.5*data/N.cos(th2D)
        elif field in ('ra', 'ratio'):
            label = 'Ratio'
            data = self.gr.vphi**2#/(self.gr.vphi**2+\
                   #self.gr.vtheta**2+self.gr.vr**2)
            denom = N.mean(self.gr.vphi**2+ self.gr.vtheta**2+self.gr.vr**2, 
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
            data = beta * N.ones_like(self.gr.vr)#* self.gr.vr
        elif field in ('angular'):
            label = 'Angular momentum'
            th2D = N.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
            rr2D = N.zeros_like(th2D)
            rho2D = N.zeros_like(th2D)
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
            s2D = rr2D * N.sin(th2D)
            norm = self.gr.radius[0]**2/self.gr.ek
            data = (self.gr.vphi*s2D+1./self.gr.ek*s2D**2)/norm
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
            thlin = N.linspace(0., N.pi, self.gr.ntheta)
            th3D = N.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            vs = vr * N.sin(th3D) + vt * N.cos(th3D)
            data =  vs * vp
            denom = N.sqrt(N.mean(vs**2, axis=0)* N.mean(vp**2, axis=0))
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
            height = 2. * N.sqrt( self.gr.radius.max()**2-self.gr.radius**2 )
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
            thlin = N.linspace(0., N.pi, self.gr.ntheta)
            th3D = N.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            vs = vr * N.sin(th3D) + vt * N.cos(th3D)
            data =  rho0 * vs * vp
            denom = N.sqrt(N.mean(rho0*vs**2, axis=0)* N.mean(rho0*vp**2, axis=0))
            if labTex:
                label = r'$C_{s\phi}$'
            else:
                label = 'Csp'
        elif field in ('Cz', 'cz'):
            vr = self.gr.vr
            vt = self.gr.vtheta
            vp = self.gr.vphi.copy()
            vp = self.gr.vphi- self.gr.vphi.mean(axis=0) # convective vp
            thlin = N.linspace(0., N.pi, self.gr.ntheta)
            th3D = N.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            vz = vr * N.cos(th3D) - vt * N.sin(th3D)
            data =  vz * vp
            denom = N.sqrt(N.mean(vz**2, axis=0)* N.mean(vp**2, axis=0))
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
            thlin = N.linspace(0., N.pi, self.gr.ntheta)
            th3D = N.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            data = (vr * N.cos(th3D) - vt * N.sin(th3D))
        else:
            data, label = selectField(self.gr, field, labTex)

        if field not in ('Cr', 'cr', 'ra', 'ratio', 'Cz', 'cz', 'hz', 'jphi',
                         'rhocr', 'omeffect', 'poloidal', 'flux', 'meridional'):
            phiavg = data.mean(axis=0)
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
            fac = 2./(N.pi*(ro**2-ri**2))
            facOTC = ro**2.*(N.pi-2.*N.arcsin(self.gr.radratio))/2. \
                     -ri**2*N.sqrt(1.-self.gr.radratio**2)/self.gr.radratio
            facOTC = 1./facOTC
            facITC = ri**2*N.sqrt(1.-self.gr.radratio**2)/self.gr.radratio \
                     +(ro**2-ri**2)* N.arcsin(self.gr.radratio) \
                     -ri**2/2.*(N.pi - 2.*N.arcsin(self.gr.radratio))
            facITC = 1./facITC
            #mask = N.where(denom == 0, 1, 0)
            phiavg = data.mean(axis=0)

            TC = N.array([], dtype=self.precision)
            outTC = N.array([], dtype=self.precision)
            inTC = N.array([], dtype=self.precision)
            denomTC = N.array([], dtype=self.precision)
            denomoutTC = N.array([], dtype=self.precision)
            denominTC = N.array([], dtype=self.precision)
            integ = N.array([], dtype=self.precision)
            for k, th in enumerate(self.gr.colatitude):
                rr = self.gr.radius[::-1]
                dat = phiavg[k, ::-1] * rr
                dat2 = denom[k, ::-1] * rr
                corr = intcheb(dat, self.gr.nr-1, ri, ro)
                TC = N.append(TC, corr)
                corr2 = intcheb(dat2, self.gr.nr-1, ri, ro)
                denomTC = N.append(denomTC, corr2)
                if th >= N.arcsin(self.gr.radratio)  and \
                   th <= N.pi - N.arcsin(self.gr.radratio):
                    # Outside tangent cylinder
                    val = trapz(dat[rr >= ri/N.sin(th)], rr[rr >= ri/N.sin(th)])
                    outTC = N.append(outTC, val)
                    integ = N.append(integ, th)
                    val2 = trapz(dat2[rr >= ri/N.sin(th)], rr[rr >= ri/N.sin(th)])
                    denomoutTC = N.append(denomoutTC, val2)
                    # Inside tangent cylinder
                    val = trapz(dat[rr < ri/N.sin(th)], rr[rr < ri/N.sin(th)])
                    inTC = N.append(inTC, val)
                    val2 = trapz(dat2[rr < ri/N.sin(th)], rr[rr < ri/N.sin(th)])
                    denominTC = N.append(denominTC, val2)
                else:
                    val= intcheb(dat, self.gr.nr-1, ri, ro)
                    inTC = N.append(inTC, val)
                    val2= intcheb(dat2, self.gr.nr-1, ri, ro)
                    denominTC = N.append(denominTC, val2)

            num = fac*trapz(TC, self.gr.colatitude)
            den = fac*trapz(denomTC, self.gr.colatitude)
            print('Correlation', num/den)
            num = facOTC*trapz(outTC, integ)
            den = facOTC*trapz(denomoutTC, integ)
            print('Correlation out TC', num/den)
            num = facITC*trapz(inTC, self.gr.colatitude)
            den = facITC*trapz(denominTC, self.gr.colatitude)
            print('Correlation in TC', num/den)

            mask = N.where(denom == 0, 1, 0)
            phiavg /= (denom + mask)
            #phiavg /= den

        th = N.linspace(0, N.pi, self.gr.ntheta)
        rr, tth = N.meshgrid(self.gr.radius, th)
        xx = rr * N.sin(tth)
        yy = rr * N.cos(tth)

        if tit:
            if cbar:
                fig = P.figure(figsize=(5,7.5))
                ax = fig.add_axes([0.01, 0.01, 0.69, 0.91])
            else:
                fig = P.figure(figsize=(3.5,7.5))
                ax = fig.add_axes([0.01, 0.01, 0.98, 0.91])
            ax.set_title(label, fontsize=24)
        else:
            if cbar:
                fig = P.figure(figsize=(5,7))
                ax = fig.add_axes([0.01, 0.01, 0.69, 0.98])
            else:
                fig = P.figure(figsize=(3.5,7))
                ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])

        cmap = P.get_cmap(cm)
        if vmax is not None and vmin is not None:
            normed = False
            cs = N.linspace(vmin, vmax, levels)
            im = ax.contourf(xx, yy, phiavg, cs, cmap=cmap, extend='both')
        else:
            cs = levels
            im = ax.contourf(xx, yy, phiavg, cs, cmap=cmap)
            #im = ax.pcolormesh(xx, yy, phiavg, cmap=cmap, antialiased=True)

        if pol:
            ax.contour(xx, yy, poloLines, polLevels, colors=['k'],
                       linewidths=[0.8, 0.8])
        elif tor:
            toroLines = self.gr.Bphi.mean(axis=0)
            ax.contour(xx, yy, toroLines, polLevels, colors=['k'], 
                       linewidths=[0.8])
        elif mer:
            maxMeri = abs(meriLines).max()
            minMeri = -maxMeri
            lev = N.linspace(minMeri, maxMeri, merLevels)
            im2 = ax.contour(xx, yy, meriLines, lev, colors=['k'],
                       linewidths=[0.8])
        ax.plot(self.gr.radius[0]*N.sin(th), self.gr.radius[0]*N.cos(th), 'k-')
        ax.plot(self.gr.radius[-1]*N.sin(th), self.gr.radius[-1]*N.cos(th), 'k-')
        ax.plot([0., 0], [self.gr.radius[-1], self.gr.radius[0]], 'k-')
        ax.plot([0., 0], [-self.gr.radius[-1], -self.gr.radius[0]], 'k-')
        #ax.axvline(0., color='k', ymin=self.gr.radius[-1])
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
        
        if field == 'b2':
            fig1 = P.figure()
            ax1 = fig1.add_subplot(111)
            ax1.plot(self.gr.radius, phiavg.mean(axis=0), 'k-')
            ax1.set_xlim(self.gr.radius.min(), self.gr.radius.max())
            ax1.set_xlabel('Radius')
            ax1.set_ylabel('Amplitude of B')
            if hasattr(self.gr, 'nVarCond'):
                if self.gr.nVarCond == 2:
                    ax1.axvline(self.gr.con_RadRatio*self.gr.radius[0],
                              color='k', linestyle='--')


        # Variable conductivity: add a dashed line
        if hasattr(self.gr, 'nVarCond'):
            if self.gr.nVarCond == 2:
                radi = self.gr.con_RadRatio * self.gr.radius[0]
                ax.plot(radi*N.sin(th), radi*N.cos(th), 'k--')

        # Normalisation of the contours around zero
        if field not in ['entropy', 's', 'S', 'u2', 'b2', 'nrj'] \
            and normed is True:
            im.set_clim(-max(abs(phiavg.max()), abs(phiavg.min())), 
                         max(abs(phiavg.max()), abs(phiavg.min())))

    def slice(self, field='Bphi', lon_0=0., levels=12, cm='RdYlBu_r', 
              normed=True, vmin=None, vmax=None, cbar=True, tit=True,
              grid=False, nGridLevs=16):
        """
        Plot an azimuthal slice of a given field.

        :param field: the field to display
        :param lon_0: the longitude of the slice, or an array of longitudes
        :param levels: the number of contour levels
        :param cm: the name of the colormap
        :param cbar: display/hide the colorbar
        :param tit: display/hide the title
        :param grid: display/hide the grid
        :param nGridLevs: number of grid levels
        """
        if field in ('Vs', 'vs'):
            if labTex:
                label = r'$v_s$'
            else:
                label = 'vs'
            vr = self.gr.vr
            vt = self.gr.vtheta
            thlin = N.linspace(0., N.pi, self.gr.ntheta)
            th3D = N.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            data = vr * N.sin(th3D) + vt * N.cos(th3D)
        elif field in ('Vz', 'vz'):
            if labTex:
                label = '$v_z$'
            else:
                label = 'vz'
            vr = self.gr.vr
            vt = self.gr.vtheta
            thlin = N.linspace(0., N.pi, self.gr.ntheta)
            th3D = N.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            data = vr * N.cos(th3D) - vt * N.sin(th3D)
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
            thlin = N.linspace(0., N.pi, self.gr.ntheta)
            th3D = N.zeros((self.gr.npI, self.gr.ntheta, self.gr.nr),
                           dtype=self.precision)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = thlin[i]
            data = vr * N.cos(th3D) - vt * N.sin(th3D)
        elif field in ('ohm'):
            label = 'Ohmic dissipation/1e6'
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            dphi = 2.*N.pi/self.gr.nphi

            Op = (N.roll(self.gr.Btheta,-1,axis=2)-N.roll(self.gr.Btheta,1,axis=2))/\
                 (N.roll(rr3D,-1,axis=2)-N.roll(rr3D,1,axis=2)) + \
                 self.gr.Btheta/rr3D - \
                 (N.roll(self.gr.Br,-1,axis=1)-N.roll(self.gr.Br,1,axis=1))/\
                 (rr3D*(N.roll(th3D,-1,axis=1)-N.roll(th3D,1,axis=1)))
            Ot = (N.roll(self.gr.Br,-1,axis=0)-N.roll(self.gr.Br,1,axis=0))/\
                  (rr3D*N.sin(th3D)*dphi) - \
                 (N.roll(self.gr.Bphi,-1,axis=2)-N.roll(self.gr.Bphi,1,axis=2))/\
                 (N.roll(rr3D,-1,axis=2)-N.roll(rr3D,1,axis=2)) - \
                 self.gr.Bphi/rr3D
            Or = (N.roll(self.gr.Bphi,-1,axis=1)-N.roll(self.gr.Bphi,1,axis=1))/\
                 (rr3D*(N.roll(th3D,-1,axis=1)-N.roll(th3D,1,axis=1))) + \
                 N.cos(th3D)*self.gr.Bphi/(rr3D*N.sin(th3D)) - \
                 (N.roll(self.gr.Btheta,-1,axis=0)-N.roll(self.gr.Btheta,1,axis=0))/\
                 (rr3D*N.sin(th3D)*dphi)

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
            rr3D = N.zeros_like(self.gr.Bphi)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            data = -thetaderavg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, order=2)/rr3D/2.0
            label = 'Lati. mag. pres. force'
        elif field in ('mag_pres_force_p'):
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = -phideravg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, self.gr.minc)/(rr3D*N.sin(th3D))/2.0
            label = 'Longi. mag. pres. force'
        elif field in ('mag_tens_force_r'):
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = self.gr.Br * rderavg(self.gr.Br, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Br, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Br, self.gr.minc) / N.sin(th3D) / rr3D - \
                   (self.gr.Btheta**2 + self.gr.Bphi**2) / rr3D
            label = 'Rad. tens. force'
        elif field in ('mag_tens_force_t'):
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = self.gr.Br * rderavg(self.gr.Btheta, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Btheta, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Btheta, self.gr.minc) / N.sin(th3D) / rr3D + \
                   self.gr.Btheta * self.gr.Br / rr3D - \
                   self.gr.Bphi**2 * N.arctan(th3D) / rr3D
            label = 'Lati. tens. force'
        elif field in ('mag_tens_force_p'):
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = self.gr.Br * rderavg(self.gr.Bphi, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Bphi, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Bphi, self.gr.minc) / N.sin(th3D) / rr3D+ \
                   self.gr.Bphi * self.gr.Br / rr3D + \
                   self.gr.Bphi * self.gr.Btheta * N.arctan(th3D) / rr3D
            label = 'Longi. tens. force'
        elif field in ('Lorentz_r'):
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = -rderavg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, eta=self.gr.radratio)/2.0 + \
                   self.gr.Br * rderavg(self.gr.Br, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Br, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Br, self.gr.minc) / N.sin(th3D) / rr3D - \
                   (self.gr.Btheta**2 + self.gr.Bphi**2) / rr3D
            label = 'Radial Loretz force'
        elif field in ('Lorentz_t'):
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = -thetaderavg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, order=2)/rr3D/2.0 + \
                   self.gr.Br * rderavg(self.gr.Btheta, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Btheta, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Btheta, self.gr.minc) / N.sin(th3D) / rr3D + \
                   self.gr.Btheta * self.gr.Br / rr3D - \
                   self.gr.Bphi**2 * N.arctan(th3D) / rr3D
            label = 'Lati. Loretz force'
        elif field in ('Lorentz_p'):
            th3D = N.zeros_like(self.gr.Bphi)
            rr3D = N.zeros_like(th3D)
            for i in range(self.gr.nr):
                rr3D[:, :, i] = self.gr.radius[i]
            for i in range(self.gr.ntheta):
                th3D[:, i, :] = self.gr.colatitude[i]
            data = -phideravg(self.gr.Br**2+self.gr.Btheta**2+self.gr.Bphi**2, self.gr.minc)/(rr3D*N.sin(th3D))/2.0 + \
                   self.gr.Br * rderavg(self.gr.Bphi, eta=self.gr.radratio) + \
                   self.gr.Btheta * thetaderavg(self.gr.Bphi, order=2) / rr3D + \
                   self.gr.Bphi * phideravg(self.gr.Bphi, self.gr.minc) / N.sin(th3D) / rr3D+ \
                   self.gr.Bphi * self.gr.Br / rr3D + \
                   self.gr.Bphi * self.gr.Btheta * N.arctan(th3D) / rr3D
            label = 'Longi. Loretz force'
        else:
            data, label = selectField(self.gr, field, labTex)

        data = symmetrize(data, self.gr.minc)

        th = N.linspace(N.pi/2, -N.pi/2, self.gr.ntheta)
        rr, tth = N.meshgrid(self.gr.radius, th)
        xx = rr * N.cos(tth)
        yy = rr * N.sin(tth)
        phi = N.linspace(0., 360, self.gr.nphi)

        lon_0 = N.asarray(lon_0)
        cmap = P.get_cmap(cm)

        if len(lon_0) > 1:
            fig = P.figure(figsize=(2.5*len(lon_0), 5.1))
            fig.subplots_adjust(top=0.99, bottom=0.01, right=0.99, left=0.01,
                              wspace=0.01)
            for k, lon in enumerate(lon_0):
                ind = N.nonzero(N.where(abs(phi-lon) \
                                == min(abs(phi-lon)), 1, 0))
                indPlot = ind[0][0]
                phislice = data[indPlot, ...]
                if field == 'dvzdz':
                    phislice = zderavg(phislice, eta=self.gr.radratio,
                                       spectral=True, exclude=True)
                elif field == 'balance':
                    phislice = zderavg(phislice, eta=self.gr.radratio,
                                       spectral=True, exclude=True)
                    phislice1 = data1[indPlot, ...]
                    phislice = phislice + phislice1

                ax = fig.add_subplot(1,len(lon_0),k+1, frameon=False)
                if vmax is not None or vmin is not None:
                    normed = False
                    cs = N.linspace(vmin, vmax, levels)
                    im = ax.contourf(xx, yy, phislice, cs, cmap=cmap, 
                                     extend='both')
                else:
                    cs = levels
                    im = ax.contourf(xx, yy, phislice, cs, cmap=cmap)
                ax.plot(self.gr.radius[0]*N.cos(th), self.gr.radius[0]*N.sin(th),
                   'k-')
                ax.plot(self.gr.radius[-1]*N.cos(th), 
                        self.gr.radius[-1]*N.sin(th), 'k-')
                ax.plot([0., 0], [self.gr.radius[-1], self.gr.radius[0]], 'k-')
                ax.plot([0., 0], [-self.gr.radius[-1], -self.gr.radius[0]], 'k-')
                ax.axis('off')

                tit1 = r'$%i^\circ$' % lon
                ax.text(0.9, 0.9, tit1, fontsize=18,
                      horizontalalignment='right',
                      verticalalignment='center',
                      transform = ax.transAxes)
                #fig.colorbar(im)
                if field not in ['entropy', 's', 'S'] and normed is True:
                    im.set_clim(-max(abs(phislice.max()), abs(phislice.min())), 
                                 max(abs(phislice.max()), abs(phislice.min())))

        else:
            ind = N.nonzero(N.where(abs(phi-lon_0[0]) \
                            == min(abs(phi-lon_0[0])), 1, 0))
            indPlot = ind[0][0]
            phislice = data[indPlot, ...]
            if field == 'dvzdz':
                phislice = zderavg(phislice, eta=self.gr.radratio,
                                   spectral=True, exclude=True)
            elif field == 'balance':
                phislice = zderavg(phislice, eta=self.gr.radratio,
                                   spectral=True, exclude=True)
                phislice1 = data1[indPlot, ...]
                phislice = phislice + phislice1
            elif field == 'vs':
                th2D = N.zeros((self.gr.ntheta, self.gr.nr), dtype=self.precision)
                rr2D = N.zeros_like(th2D)
                for i in range(self.gr.ntheta):
                    th2D[i, :] = self.gr.colatitude[i]
                    rr2D[i, :] = self.gr.radius
                s2D = rr2D * N.sin(th2D)
                ro = 1./(1.-self.gr.radratio)
                coeff = -s2D/(ro**2-s2D**2)
                phislice *= coeff


            if tit:
                if cbar:
                    fig = P.figure(figsize=(5,7.5))
                    ax = fig.add_axes([0.01, 0.01, 0.69, 0.91])
                else:
                    fig = P.figure(figsize=(3.5,7.5))
                    ax = fig.add_axes([0.01, 0.01, 0.98, 0.91])
                ax.set_title(label, fontsize=24)
            else:
                if cbar:
                    fig = P.figure(figsize=(5,7))
                    ax = fig.add_axes([0.01, 0.01, 0.69, 0.98])
                else:
                    fig = P.figure(figsize=(3.5,7))
                    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])

            if vmax is not None or vmin is not None:
                normed = False
                cs = N.linspace(vmin, vmax, levels)
                im = ax.contourf(xx, yy, phislice, cs, cmap=cmap, extend='both')
            else:
                cs = levels
                im = ax.contourf(xx, yy, phislice, cs, cmap=cmap)
            ax.plot(self.gr.radius[0]*N.cos(th), self.gr.radius[0]*N.sin(th),
                   'k-', lw=1.5)
            ax.plot(self.gr.radius[-1]*N.cos(th), self.gr.radius[-1]*N.sin(th),
                   'k-', lw=1.5)
            ax.plot([0., 0], [self.gr.radius[-1], self.gr.radius[0]], 'k-', lw=1.5)
            ax.plot([0., 0], [-self.gr.radius[-1], -self.gr.radius[0]], 'k-', lw=1.5)

            if hasattr(self.gr, 'epsS'):
                rad = MagicRadial(field='anel', iplot=False)
                idx = N.nonzero(N.where(abs(rad.dsdr)==abs(rad.dsdr).min(), 1, 0))[0][0]
                ax.plot(self.gr.radius[idx]*N.cos(th), self.gr.radius[idx]*N.sin(th), 
                        'k--', lw=2)

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


def report(nvar=1, levels=16, lclean=True):
    """
    This subroutine prepares a pdf document that gather some important diagnostics

    :param lclean: clean or not the LaTeX files
    :param levels: number of contour levels
    :param nvar: number of graphic files
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

    st = "Ek = %.2e, Ra = %.2e, Pr = %.1f, $N_{\\rho}$=%.2f, $\eta$=%.1f" % \
            (s.gr.ek, s.gr.ra, s.gr.pr, s.gr.strat, s.gr.radratio)
    file.write("\\begin{center}\\begin{large}\n")
    file.write(" "+st+"\n")
    file.write("\\end{large}\\end{center}\n")

    r1 = 0.98
    r3 = 1.03 * s.gr.radratio
    r2 = (r1+r3)/2.

    s.avg(field='vp', levels=levels, cm='RdYlBu_r', normed=True)
    P.savefig('vp.png')
    P.close()

    s.avg(field='entropy', levels=levels, cm='RdYlBu_r', normed=True)
    P.savefig('entropy.png')
    P.close()

    s.equat(field='entropy', levels=levels, cm='RdYlBu_r', normed=False)
    P.savefig('equ_s.png')
    P.close()

    s.equat(field='vr', levels=levels, cm='RdYlBu_r', normed=False)
    P.savefig('equ_vr.png')
    P.close()

    s.surf(field='vp', cm='RdYlBu_r', levels=levels, r=r1, proj='moll', 
           normed=False)
    P.savefig('surf_vp.png')
    P.close()

    s.surf(field='vp', cm='RdYlBu', levels=levels, r=r2, proj='moll', 
           normed=False)
    P.savefig('surf_vp_08.png')
    P.close()

    s.surf(field='vp', cm='RdYlBu', levels=levels, r=r3, proj='moll', 
           normed=False)
    P.savefig('surf_vp_06.png')
    P.close()

    s.surf(field='vr', cm='RdYlBu', levels=levels, r=r1, proj='moll', 
           normed=False)
    P.savefig('surf_vr.png')
    P.close()

    s.surf(field='vr', cm='RdYlBu', levels=levels, r=r2, proj='moll', 
           normed=False)
    P.savefig('surf_vr_08.png')
    P.close()

    s.surf(field='vr', cm='RdYlBu', levels=levels, r=r3, proj='moll', 
           normed=False)
    P.savefig('surf_vr_06.png')
    P.close()

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
