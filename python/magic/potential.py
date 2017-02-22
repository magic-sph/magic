# -*- coding: utf-8 -*-
from magic import npfile, MagicSetup, scanDir
from .setup import labTex, defaultCm, defaultLevels, labTex
from .libmagic import *
from .plotlib import radialContour, merContour, equatContour
import os, re, sys
if sys.version_info.major == 3:
    from legendre3 import *
elif  sys.version_info.major == 2:
    from legendre2 import *

import numpy as np
import matplotlib.pyplot as plt


class MagicPotential(MagicSetup):

    def __init__(self, field='V', datadir='.', tag=None, ave=False, ipot=None, 
                 precision='Float32'):

        if ave:
            self.name = '%s_lmr_ave' % field
        else:
            self.name = '%s_lmr_' % field

        if tag is not None:
            if ipot is not None:
                file = '%s%i.%s' % (self.name, ipot, tag)
                filename = os.path.join(datadir, file)
            else:
                pattern = os.path.join(datadir, '%s*%s' % (self.name, tag))
                files = scanDir(pattern)
                if len(files) != 0:
                    filename = files[-1]
                else:
                    print('No such tag... try again')
                    return

            if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
        else:
            if ipot is not None:
                pattern = os.path.join(datadir, '%s%i*' % (self.name, ipot))
                files = scanDir(pattern)
                filename = files[-1]
            else:
                pattern = os.path.join(datadir, '%s*' % self.name)
                files = scanDir(pattern)
                filename = files[-1]
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists(os.path.join(datadir, 'log.%s' % ending)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % ending)

        infile = npfile(filename, endian='B')

        # Read header
        self.l_max, self.n_r_max, self.n_r_ic_max, self.minc, \
                                          self.lm_max = infile.fort_read('i4')
        self.m_max = (self.l_max/self.minc)*self.minc
        self.n_m_max = self.m_max/self.minc+1
        self.ra, self.ek, self.pr, self.prmag, self.radratio, \
           self.sigma_ratio, self.omega_ma, self.omega_ic = infile.fort_read(precision)
        self.time = infile.fort_read(precision)
        dat = infile.fort_read(precision)
        # Read radius and density
        self.radius = dat[:self.n_r_max]
        self.rho0 = dat[self.n_r_max:]

        # Read field in the outer core
        self.pol = infile.fort_read('Complex32')
        self.pol = self.pol.reshape((self.n_r_max, self.lm_max))
        self.pol = self.pol.T
        if ( field != 'T' and field != 'Xi' ):
            self.tor = infile.fort_read('Complex32')
            self.tor = self.tor.reshape((self.n_r_max, self.lm_max))
            self.tor = self.tor.T


        self.n_theta_max = 3*self.l_max/2
        self.n_phi_max = 2*self.n_theta_max
        self.nphi = self.n_phi_max+1
        self.sp = legendre
        self.sp.inittransform(self.l_max, self.minc, self.lm_max, 
                              self.n_theta_max)
        self.colat = self.sp.sinth

        self.idx = np.zeros((self.l_max+1, self.m_max+1), 'i')
        self.ell = np.zeros(self.lm_max, 'i')
        self.idx[0:self.l_max+2, 0] = np.arange(self.l_max+1)
        self.ell[0:self.l_max+2] = np.arange(self.l_max+2)
        k = self.l_max+1
        for m in range(self.minc, self.l_max+1, self.minc):
            for l in range(m, self.l_max+1):
                self.idx[l, m] = k
                self.ell[self.idx[l,m]] = l
                k +=1

    def avg(self, field='vphi', levels=defaultLevels, cm=defaultCm, normed=True,
            vmax=None, vmin=None, cbar=True, tit=True):

        phiavg = np.zeros((self.n_theta_max, self.n_r_max), 'f')
        if field in ('T', 'temp', 'S', 'entropy'):
            for i in range(self.n_r_max):
                out = self.sp.specspat_scal(1, self.n_theta_max, 1, self.pol[:, i],
                                         self.lm_max)
                phiavg[:, i] = out[0,:].real

            label = 'T'
        elif field in ('vr', 'Vr', 'ur', 'Ur'):
            for i in range(self.n_r_max):
                field = self.pol[:,i]*self.ell*(self.ell+1)/self.radius[i]**2\
                        /self.rho0[i]

                out = self.sp.specspat_scal(1, self.n_theta_max, 1, field, self.lm_max)
                phiavg[:, i] = out[0,:].real
            if labTex:
                label = r'$v_r$'
            else:
                label = 'vr'
        elif field in ('vt', 'Vt', 'ut', 'Ut', 'utheta', 'vtheta'):
            field = rderavg(self.pol, self.radratio)
            for i in range(self.n_r_max):
                vt, vp = self.sp.specspat_vec(1, self.n_theta_max, 1,
                                       field[:, i], self.tor[:, i], self.lm_max)
                phiavg[:, i] = vt[0,:].real/self.radius[i]/self.rho0[i]
            if labTex:
                label = r'$v_\theta$'
            else:
                label = 'vtheta'
        elif field in ('vp', 'Vp', 'up', 'Up', 'uphi', 'vphi'):
            field = rderavg(self.pol, self.radratio)
            for i in range(self.n_r_max):
                vt, vp = self.sp.specspat_vec(1, self.n_theta_max, 1,
                                       field[:, i], self.tor[:, i], self.lm_max)
                phiavg[:, i] = vp[0,:].real/self.radius[i]/self.rho0[i]
            if labTex:
                label = r'$v_\phi$'
            else:
                label = 'vphi'
        elif field in ('br', 'Br'):
            for i in range(self.n_r_max):
                field = self.pol[:,i]*self.ell*(self.ell+1)/self.radius[i]**2

                out = self.sp.specspat_scal(1, self.n_theta_max, 1, field, self.lm_max)
                phiavg[:, i] = out[0,:].real
            if labTex:
                label = r'$B_r$'
            else:
                label = 'br'
        elif field in ('bt', 'Bt', 'Btheta', 'btheta'):
            field = rderavg(self.pol, self.radratio)
            for i in range(self.n_r_max):
                bt, bp = self.sp.specspat_vec(1, self.n_theta_max, 1,
                                       field[:, i], self.tor[:, i], self.lm_max)
                phiavg[:, i] = bt[0,:].real/self.radius[i]
            if labTex:
                label = r'$B_\theta$'
            else:
                label = 'Btheta'
        elif field in ('bp', 'Bp', 'bphi', 'Bphi'):
            field = rderavg(self.pol, self.radratio)
            for i in range(self.n_r_max):
                bt, bp = self.sp.specspat_vec(1, self.n_theta_max, 1,
                                       field[:, i], self.tor[:, i], self.lm_max)
                phiavg[:, i] = bp[0,:].real/self.radius[i]
            if labTex:
                label = r'$B_\phi$'
            else:
                label = 'Bphi'
        
        if field in ('temperature', 't', 'T', 'entropy', 's', 'S', 'u2', 'b2', 'nrj'):
            normed = False

        fig, xx, yy = merContour(phiavg, self.radius, label, levels, cm, normed,
                                 vmax, vmin, cbar, tit)

    def equat(self, field='vr', levels=defaultLevels, cm=defaultCm, normed=True,
              vmax=None, vmin=None, cbar=True, tit=True, normRad=False):

        equator = np.zeros((self.n_phi_max, self.n_r_max), 'f')

        if field in ('temperature', 't', 'T', 'entropy', 's', 'S'):
            for i in range(self.n_r_max):
                out = self.sp.specspat_equat_scal(self.n_phi_max, self.pol[:, i],
                                                  self.lm_max)
                out = np.fft.ifft(out)*self.n_phi_max
                equator[:, i] = out.real
            label = 'T'
        elif field in ('vr', 'Vr', 'ur', 'Ur'):
            for i in range(self.n_r_max):
                field = self.ell*(self.ell+1)/self.radius[i]**2/self.rho0[i]*\
                        self.pol[:, i]
                out = self.sp.specspat_equat_scal(self.n_phi_max, field, self.lm_max)
                out = np.fft.ifft(out)*self.n_phi_max
                equator[:, i] = out.real
            if labTex:
                label = r'$v_r$'
            else:
                label = 'vr'
        elif field in ('vt', 'Vt', 'ut', 'Ut', 'vtheta', 'utheta'):
            field = rderavg(self.pol, self.radratio)
            for i in range(self.n_r_max):
                vt, vp = self.sp.specspat_equat_vec(self.n_phi_max, field[:, i],
                                                    self.tor[:, i], self.lm_max)
                rprof = np.fft.ifft(vt, axis=0)*self.n_phi_max
                equator[:, i] = rprof.real/self.radius[i]/self.rho0[i]
            if labTex:
                label = r'$v_\theta$'
            else:
                label = 'vtheta'
        elif field in ('vp', 'Vp', 'up', 'Up', 'vphi', 'uphi'):
            field = rderavg(self.pol, self.radratio)
            for i in range(self.n_r_max):
                vt, vp = self.sp.specspat_equat_vec(self.n_phi_max, field[:, i],
                                                    self.tor[:, i], self.lm_max)
                rprof = np.fft.ifft(vp, axis=0)*self.n_phi_max
                equator[:, i] = rprof.real/self.radius[i]/self.rho0[i]
            if labTex:
                label = r'$v_\phi$'
            else:
                label = 'vphi'
        elif field in ('Br', 'br'):
            for i in range(self.n_r_max):
                field = self.ell*(self.ell+1)/self.radius[i]**2*self.pol[:, i]
                out = self.sp.specspat_equat_scal(self.n_phi_max, field, self.lm_max)
                out = np.fft.ifft(out)*self.n_phi_max
                equator[:, i] = out.real
            if labTex:
                label = r'$B_r$'
            else:
                label = 'Br'
        elif field in ('bt', 'Bt', 'Btheta', 'btheta'):
            field = rderavg(self.pol, self.radratio)
            for i in range(self.n_r_max):
                bt, bp = self.sp.specspat_equat_vec(self.n_phi_max, field[:, i],
                                                    self.tor[:, i], self.lm_max)
                rprof = np.fft.ifft(bt, axis=0)*self.n_phi_max
                equator[:, i] = rprof.real/self.radius[i]
            if labTex:
                label = r'$B_\theta$'
            else:
                label = 'Btheta'
        elif field in ('bp', 'Bp', 'Bphi', 'bphi'):
            field = rderavg(self.pol, self.radratio)
            for i in range(self.n_r_max):
                bt, bp = self.sp.specspat_equat_vec(self.n_phi_max, field[:, i],
                                                    self.tor[:, i], self.lm_max)
                rprof = np.fft.ifft(bp, axis=0)*self.n_phi_max
                equator[:, i] = rprof.real/self.radius[i]
            if labTex:
                label = r'$B_\phi$'
            else:
                label = 'Bphi'

        equator = symmetrize(equator, self.minc)

        if field in ('temperature', 't', 'T', 'entropy', 's', 'S', 'u2', 'b2', 'nrj'):
            normed = False

        fig, xx, yy = equatContour( equator, self.radius, label, levels, cm,
                                    normed, vmax, vmin, cbar, tit, normRad)
        ax = fig.get_axes()[0]

        # Variable conductivity: add a dashed line
        if hasattr(self, 'nVarCond'):
            if self.nVarCond == 2:
                radi = self.con_RadRatio * self.radius[0]
                ax.plot(radi*np.cos(phi), radi*np.sin(phi), 'k--', lw=1.5)

    def surf(self, field='vr', proj='hammer', lon_0=0., r=0.85, vmax=None,
             vmin=None, lat_0=30., levels=defaultLevels, cm=defaultCm, lon_shift=0,
             normed=True, cbar=True, tit=True, lines=False):

        r /= (1-self.radratio) # as we give a normalised radius
        ind = np.nonzero(np.where(abs(self.radius-r) \
                        == min(abs(self.radius-r)), 1, 0))
        indPlot = ind[0][0]
        rad = self.radius[indPlot] * (1.-self.radratio)

        if field in ('T', 'temp', 'S', 'entropy'):
            rprof = self.sp.specspat_scal(self.n_phi_max, self.n_theta_max,
                                          self.n_m_max, self.pol[:, indPlot],
                                          self.lm_max)
            rprof = np.fft.ifft(rprof, axis=0)*self.n_phi_max
            rprof = rprof.real
            label = 'T'
        elif field in ('vr', 'Vr', 'ur', 'Ur'):
            field = self.ell*(self.ell+1)/self.radius[indPlot]**2/self.rho0[indPlot]*\
                    self.pol[:, indPlot]
            rprof = self.sp.specspat_scal(self.n_phi_max, self.n_theta_max,
                                          self.n_m_max, field, self.lm_max)
            rprof = np.fft.ifft(rprof, axis=0)*self.n_phi_max
            rprof = rprof.real
            if labTex:
                label = r'$v_r$'
            else:
                label = 'vr'
        elif field in ('vt', 'Vt', 'ut', 'Ut', 'vtheta', 'utheta'):
            field = rderavg(self.pol, self.radratio)
            field = field[:, indPlot]
            vt, vp = self.sp.specspat_vec(self.n_phi_max, self.n_theta_max,
                                          self.n_m_max, field, self.tor[:, indPlot],
                                          self.lm_max)
            rprof = np.fft.ifft(vt, axis=0)*self.n_phi_max
            rprof = rprof.real/self.radius[indPlot]/self.rho0[indPlot]
            if labTex:
                label = r'$v_\theta$'
            else:
                label = 'vtheta'
        elif field in ('vp', 'Vp', 'up', 'Up', 'uphi', 'vphi'):
            field = rderavg(self.pol, self.radratio)
            field = field[:, indPlot]
            vt, vp = self.sp.specspat_vec(self.n_phi_max, self.n_theta_max,
                                          self.n_m_max, field, self.tor[:, indPlot],
                                          self.lm_max)
            rprof = np.fft.ifft(vp, axis=0)*self.n_phi_max
            rprof = rprof.real/self.radius[indPlot]/self.rho0[indPlot]
            if labTex:
                label = r'$v_\phi$'
            else:
                label = 'vphi'
        elif field in ('br', 'Br'):
            field = self.ell*(self.ell+1)/self.radius[indPlot]**2*\
                    self.pol[:, indPlot]
            rprof = self.sp.specspat_scal(self.n_phi_max, self.n_theta_max,
                                          self.n_m_max, field, self.lm_max)
            rprof = np.fft.ifft(rprof, axis=0)*self.n_phi_max
            rprof = rprof.real
            if labTex:
                label = r'$B_r$'
            else:
                label = 'Br'
        elif field in ('bt', 'Bt', 'Btheta', 'btheta'):
            field = rderavg(self.pol, self.radratio)
            field = field[:, indPlot]
            bt, bp = self.sp.specspat_vec(self.n_phi_max, self.n_theta_max,
                                          self.n_m_max, field, self.tor[:, indPlot],
                                          self.lm_max)
            rprof = np.fft.ifft(bt, axis=0)*self.n_phi_max
            rprof = rprof.real/self.radius[indPlot]
            if labTex:
                label = r'$B_\theta$'
            else:
                label = 'Btheta'
        elif field in ('bp', 'Bp', 'bphi', 'Bphi'):
            field = rderavg(self.pol, self.radratio)
            field = field[:, indPlot]
            bt, bp = self.sp.specspat_vec(self.n_phi_max, self.n_theta_max,
                                          self.n_m_max, field, self.tor[:, indPlot],
                                          self.lm_max)
            rprof = np.fft.ifft(bp, axis=0)*self.n_phi_max
            rprof = rprof.real/self.radius[indPlot]
            if labTex:
                label = r'$B_\phi$'
            else:
                label = 'Bphi'


        if field in ('temperature', 't', 'T', 'entropy', 's', 'S', 'u2', 'b2', 'nrj'):
            normed = False

        fig = radialContour(rprof, rad, label, proj, lon_0, vmax, vmin, lat_0, levels,
                            cm, normed, cbar, tit, lines)


if __name__ == '__main__':

    p = MagicPotential(field='B', ave=True)
    p.surf(field='br', r=0.9)
    p.avg(field='br')
    p.equat(field='br')

    #p = MagicPotential(field='T', ave=True)
    #p.surf(field='T', r=0.63)
    #p.avg(field='T')
    #p.equat(field='T')

    plt.show()
