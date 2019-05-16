# -*- coding: utf-8 -*-
from magic import MagicSetup, scanDir
from .setup import labTex, defaultCm, defaultLevels, labTex, buildSo
from .libmagic import *
from .spectralTransforms import SpectralTransforms
from .plotlib import radialContour, merContour, equatContour
import os, re, sys, time
import numpy as np
import matplotlib.pyplot as plt
from .npfile import *

if buildSo:
    if sys.version_info.major == 3:
        import magic.lmrreader_single3 as Psngl
    elif  sys.version_info.major == 2:
        import magic.lmrreader_single2 as Psngl

    readingMode = 'f2py'
else:
    readingMode = 'python'


def getPotEndianness(filename):
    """
    This function determines the endianness of the potential files

    :param filename: input of the filename
    :type filename: str
    :returns: the endianness of the file ('B'='big_endian' or
              'l'='little_endian')
    :rtype: str
    """
    f = npfile(filename, endian='B')
    try:
        f.fort_read('i4')
        endian = 'B'
    except TypeError:
        endian = 'l'
    f.close()

    return endian


class MagicPotential(MagicSetup):
    """
    This class allows to load and display the content of the potential
    files: :ref:`V_lmr.TAG <secVpotFile>`, :ref:`B_lmr.TAG <secBpotFile>`
    and :ref:`T_lmr.TAG <secTpotFile>`. This class allows to transform
    the poloidal/toroidal potential in spectral space to the physical
    quantities in the physical space. It allows to plot radial and
    equatorial cuts as well as phi-averages.

    >>> # To read T_lmr.test
    >>> p = MagicPotential(field='T', ipot=1, tag='test')
    >>> # To read the latest V_lmr file in the working directory
    >>> p = MagicPotential(field='V')
    >>> # Get the poloidal potential (lm, nR)
    >>> wlm = p.pol
    >>> # Obtain the value of w(l=12, m=12, nR=33)
    >>> print( p.pol[p.idx[12,12], 32] )

    >>> # Possible plots
    >>> p.equat(field='vr')
    >>> p.avg(field='vp')
    >>> p.surf(field='vt', r=0.8)
    """

    def __init__(self, field='V', datadir='.', tag=None, ave=False, ipot=None,
                 precision=np.float32):
        """
        :param field: 'B', 'V' or 'T' (magnetic field, velocity field or
                      temperature)
        :type field: str
        :param datadir: the working directory
        :type datadir: str
        :param tag: if you specify a pattern, it tries to read the
                    corresponding files
        :type tag: str
        :param ave: plot a time-averaged spectrum when set to True
        :type ave: bool
        :param ipot: the number of the lmr file you want to plot
        :type ipot: int
        :param precision: single or double precision
        :type precision: str
        """

        if hasattr(self, 'radial_scheme'):
            if self.radial_scheme == 'FD':
                self.rcheb = False
            else:
                self.rcheb = True
        else:
            self.rcheb = True

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


        # Determine file endianness
        endian = getPotEndianness(filename)

        t1 = time.time()
        if readingMode  == 'python':

            infile = npfile(filename, endian=endian)

            # Read header
            self.l_max, self.n_r_max, self.n_r_ic_max, self.minc, \
                                          self.lm_max = infile.fort_read('i4')
            self.m_max = int((self.l_max/self.minc)*self.minc)
            self.n_m_max = int(self.m_max/self.minc+1)
            self.ra, self.ek, self.pr, self.prmag, self.radratio, self.sigma_ratio, \
                         self.omega_ma, self.omega_ic = infile.fort_read(precision)
            self.time = infile.fort_read(precision)
            dat = infile.fort_read(precision)

            # Read radius and density
            self.radius = dat[:self.n_r_max]
            self.rho0 = dat[self.n_r_max:]

            # Read field in the outer core
            self.pol = infile.fort_read(np.complex64)
            self.pol = self.pol.reshape((self.n_r_max, self.lm_max))
            self.pol = self.pol.T
            if ( field != 'T' and field != 'Xi' ):
                self.tor = infile.fort_read(np.complex64)
                self.tor = self.tor.reshape((self.n_r_max, self.lm_max))
                self.tor = self.tor.T

            infile.close()

        else: # F2py reader

            if ( field != 'T' and field != 'Xi' ):
                l_read_tor = True
            else:
                l_read_tor = False

            Prd = Psngl.potreader_single
            Prd.readpot(filename, endian, l_read_tor)
            self.n_r_max = Prd.n_r_max
            self.l_max = Prd.l_max
            self.n_r_ic_max = Prd.n_r_ic_max
            self.minc = Prd.minc
            self.lm_max = Prd.lm_max
            self.m_max = int((self.l_max/self.minc)*self.minc)
            self.n_m_max = int(self.m_max/self.minc+1)
            self.ra = Prd.ra
            self.ek = Prd.ek
            self.radratio = Prd.radratio
            self.sigma_ratio = Prd.sigma_ratio
            self.prmag = Prd.prmag
            self.pr = Prd.pr
            self.omega_ic = Prd.omega_ic
            self.omega_ma = Prd.omega_ma
            self.radius = Prd.radius
            self.rho0 = Prd.rho0

            self.pol = Prd.pol
            if ( field != 'T' and field != 'Xi' ):
                self.tor = Prd.tor
        t2 = time.time()
        print('Time to read %s: %.2f' % (filename, t2-t1))

        self.n_theta_max = int(3*self.l_max/2)
        if self.n_theta_max % 2: # odd number
            self.n_theta_max += 1
        self.n_phi_max = int(2*self.n_theta_max/self.minc)
        t1 = time.time()
        self.sh = SpectralTransforms(l_max=self.l_max, minc=self.minc,
                                     lm_max=self.lm_max, n_theta_max=self.n_theta_max)
        t2 = time.time()
        print('Time to set up the spectral transforms: %.2f' % (t2-t1))
        self.colat = self.sh.colat

        self.idx = self.sh.idx
        self.ell = self.sh.ell

    def avg(self, field='vphi', levels=defaultLevels, cm=defaultCm, normed=True,
            vmax=None, vmin=None, cbar=True, tit=True):
        """
        Plot the azimutal average of a given field.

           >>> p = MagicPotential(field='V')
           >>> # Axisymmetric zonal flows, 65 contour levels
           >>> p.avg(field='vp', levels=65, cm='seismic')

           >>> # Minimal plot (no cbar, not title)
           >>> p.avg(field='vr', tit=False, cbar=False)

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
        """

        phiavg = np.zeros((self.n_theta_max, self.n_r_max), np.float32)
        t1 = time.time()
        if field in ('T', 'temp', 'S', 'entropy'):
            for i in range(self.n_r_max):
                phiavg[:, i] = self.sh.spec_spat(self.pol[:, i], l_axi=True)

            label = 'T'
        elif field in ('vr', 'Vr', 'ur', 'Ur'):
            for i in range(self.n_r_max):
                field = self.pol[:,i]*self.ell*(self.ell+1)/self.radius[i]**2\
                        /self.rho0[i]
                phiavg[:, i] = self.sh.spec_spat(field, l_axi=True)
            if labTex:
                label = r'$v_r$'
            else:
                label = 'vr'
        elif field in ('vt', 'Vt', 'ut', 'Ut', 'utheta', 'vtheta'):
            field = rderavg(self.pol, self.radratio, spectral=self.rcheb)
            for i in range(self.n_r_max):
                vt, vp = self.sh.spec_spat( field[:, i], self.tor[:, i], l_axi=True)
                phiavg[:, i] = vt/self.radius[i]/self.rho0[i]
            if labTex:
                label = r'$v_\theta$'
            else:
                label = 'vtheta'
        elif field in ('vp', 'Vp', 'up', 'Up', 'uphi', 'vphi'):
            field = rderavg(self.pol, self.radratio, spectral=self.rcheb)
            for i in range(self.n_r_max):
                vt, vp = self.sh.spec_spat( field[:, i], self.tor[:, i], l_axi=True)
                phiavg[:, i] = vp/self.radius[i]/self.rho0[i]
            if labTex:
                label = r'$v_\phi$'
            else:
                label = 'vphi'
        elif field in ('br', 'Br'):
            for i in range(self.n_r_max):
                field = self.pol[:,i]*self.ell*(self.ell+1)/self.radius[i]**2

                phiavg[:, i] = self.sh.spec_spat(field, l_axi=True)
            if labTex:
                label = r'$B_r$'
            else:
                label = 'br'
        elif field in ('bt', 'Bt', 'Btheta', 'btheta'):
            field = rderavg(self.pol, self.radratio, spectral=self.rcheb)
            for i in range(self.n_r_max):
                bt, bp = self.sh.spec_spat( field[:, i], self.tor[:, i], l_axi=True)
                phiavg[:, i] = bt.real/self.radius[i]
            if labTex:
                label = r'$B_\theta$'
            else:
                label = 'Btheta'
        elif field in ('bp', 'Bp', 'bphi', 'Bphi'):
            field = rderavg(self.pol, self.radratio, spectral=self.rcheb)
            for i in range(self.n_r_max):
                bt, bp = self.sh.spec_spat( field[:, i], self.tor[:, i], l_axi=True)
                phiavg[:, i] = bp.real/self.radius[i]
            if labTex:
                label = r'$B_\phi$'
            else:
                label = 'Bphi'
        t2 = time.time()
        print('Transform time (avg): %.2f' % (t2-t1))

        if field in ('temperature', 'entropy', 's', 'S', 'u2', 'b2', 'nrj'):
            normed = False

        fig, xx, yy = merContour(phiavg, self.radius, label, levels, cm, normed,
                                 vmax, vmin, cbar, tit)

    def equat(self, field='vr', levels=defaultLevels, cm=defaultCm, normed=True,
              vmax=None, vmin=None, cbar=True, tit=True, normRad=False):
        """
        Plot the equatorial cut of a given field

           >>> p = MagicPotential(field='B')
           >>> # Equatorial cut of the Br
           >>> p.equat(field='Br')

           >>> # Normalise the contour levels radius by radius
           >>> p.equat(field='Bphi', normRad=True)

        :param field: the name of the input physical quantity you want to
                      display
        :type field: str
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
        """

        equator = np.zeros((self.n_phi_max, self.n_r_max), np.float32)
        t1 = time.time()
        if field in ('temperature', 't', 'T', 'entropy', 's', 'S'):
            for i in range(self.n_r_max):
                equator[:, i] = self.sh.spec_spat_equat(self.pol[:, i])
            label = 'T'
        elif field in ('vr', 'Vr', 'ur', 'Ur'):
            for i in range(self.n_r_max):
                field = self.ell*(self.ell+1)/self.radius[i]**2/self.rho0[i]*\
                        self.pol[:, i]
                equator[:, i] = self.sh.spec_spat_equat(field)
            if labTex:
                label = r'$v_r$'
            else:
                label = 'vr'
        elif field in ('vt', 'Vt', 'ut', 'Ut', 'vtheta', 'utheta'):
            field = rderavg(self.pol, self.radratio, spectral=self.rcheb)
            for i in range(self.n_r_max):
                vt, vp = self.sh.spec_spat_equat(field[:, i], self.tor[:, i])
                equator[:, i] = vt/self.radius[i]/self.rho0[i]
            if labTex:
                label = r'$v_\theta$'
            else:
                label = 'vtheta'
        elif field in ('vp', 'Vp', 'up', 'Up', 'vphi', 'uphi'):
            field = rderavg(self.pol, self.radratio, spectral=self.rcheb)
            for i in range(self.n_r_max):
                vt, vp = self.sh.spec_spat_equat(field[:, i], self.tor[:, i])
                equator[:, i] = vp/self.radius[i]/self.rho0[i]
            if labTex:
                label = r'$v_\phi$'
            else:
                label = 'vphi'
        elif field in ('Br', 'br'):
            for i in range(self.n_r_max):
                field = self.ell*(self.ell+1)/self.radius[i]**2*self.pol[:, i]
                equator[:, i] = self.sh.spec_spat_equat(field)
            if labTex:
                label = r'$B_r$'
            else:
                label = 'Br'
        elif field in ('bt', 'Bt', 'Btheta', 'btheta'):
            field = rderavg(self.pol, self.radratio, spectral=self.rcheb)
            for i in range(self.n_r_max):
                bt, bp = self.sh.spec_spat_equat(field[:, i], self.tor[:, i])
                equator[:, i] = bt/self.radius[i]
            if labTex:
                label = r'$B_\theta$'
            else:
                label = 'Btheta'
        elif field in ('bp', 'Bp', 'Bphi', 'bphi'):
            field = rderavg(self.pol, self.radratio, spectral=self.rcheb)
            for i in range(self.n_r_max):
                bt, bp = self.sh.spec_spat_equat(field[:, i], self.tor[:, i])
                equator[:, i] = bp/self.radius[i]
            if labTex:
                label = r'$B_\phi$'
            else:
                label = 'Bphi'
        t2 = time.time()
        print('Transform time (equat): %.2f' % (t2-t1))

        equator = symmetrize(equator, self.minc)

        if field in ('temperature', 't', 'T', 'entropy', 's', 'S', 'u2', 'b2', 'nrj'):
            normed = False

        fig, xx, yy = equatContour( equator, self.radius, self.minc, label,
                                    levels, cm, normed, vmax, vmin, cbar, tit,
                                    normRad )
        ax = fig.get_axes()[0]

        # Variable conductivity: add a dashed line
        if hasattr(self, 'nVarCond'):
            if self.nVarCond == 2:
                radi = self.con_RadRatio * self.radius[0]
                ax.plot(radi*np.cos(phi), radi*np.sin(phi), 'k--', lw=1.5)

    def surf(self, field='vr', proj='hammer', lon_0=0., r=0.85, vmax=None,
             vmin=None, lat_0=30., levels=defaultLevels, cm=defaultCm, lon_shift=0,
             normed=True, cbar=True, tit=True, lines=False):
        """
        Plot the surface distribution of an input field at a given
        input radius (normalised by the outer boundary radius).

           >>> p = MagicPotential(field='V')
           >>> # Radial flow component at ``r=0.95 r_o``, 65 contour levels
           >>> p.surf(field='vr', r=0.95, levels=65, cm='seismic')

           >>> # Control the limit of the colormap from -1e3 to 1e3
           >>> p.surf(field='vp', r=1., vmin=-1e3, vmax=1e3, levels=33)

        :param field: the name of the field you want to display
        :type field: str
        :param proj: the type of projection. Default is Hammer, in case
                     you want to use 'ortho' or 'moll', then Basemap is
                     required.
        :type proj: str
        :param lon_0: central azimuth (only used with Basemap)
        :type lon_0: float
        :param lat_0: central latitude (only used with Basemap)
        :type lat_0: float
        :param r: the radius at which you want to display the input
                  data (in normalised units with the radius of the outer boundary)
        :type r: float
        :param levels: the number of levels in the contour
        :type levels: int
        :param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
        :type cm: str
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
        :type normed: bool
        """

        r /= (1-self.radratio) # as we give a normalised radius
        ind = np.nonzero(np.where(abs(self.radius-r) \
                        == min(abs(self.radius-r)), 1, 0))
        indPlot = ind[0][0]
        rad = self.radius[indPlot] * (1.-self.radratio)

        t1 = time.time()
        rprof = np.zeros((self.n_phi_max, self.n_theta_max), np.complex64)
        if field in ('T', 'temp', 'S', 'entropy'):
            rprof = self.sh.spec_spat(self.pol[:, indPlot])
            label = 'T'
        elif field in ('vr', 'Vr', 'ur', 'Ur'):
            field = self.ell*(self.ell+1)/self.radius[indPlot]**2/self.rho0[indPlot]*\
                    self.pol[:, indPlot]
            rprof = self.sh.spec_spat(field)
            if labTex:
                label = r'$v_r$'
            else:
                label = 'vr'
        elif field in ('vt', 'Vt', 'ut', 'Ut', 'vtheta', 'utheta'):
            field = rderavg(self.pol, self.radratio, spectral=self.rcheb)
            field = field[:, indPlot]
            vt, vp = self.sh.spec_spat(field, self.tor[:, indPlot])
            rprof = vt/self.radius[indPlot]/self.rho0[indPlot]
            if labTex:
                label = r'$v_\theta$'
            else:
                label = 'vtheta'
        elif field in ('vp', 'Vp', 'up', 'Up', 'uphi', 'vphi'):
            field = rderavg(self.pol, self.radratio, spectral=self.rcheb)
            field = field[:, indPlot]
            vt, vp  = self.sh.spec_spat(field, self.tor[:, indPlot])
            rprof = vp/self.radius[indPlot]/self.rho0[indPlot]
            if labTex:
                label = r'$v_\phi$'
            else:
                label = 'vphi'
        elif field in ('br', 'Br'):
            field = self.ell*(self.ell+1)/self.radius[indPlot]**2*\
                    self.pol[:, indPlot]
            rprof = self.sh.spec_spat(field)
            if labTex:
                label = r'$B_r$'
            else:
                label = 'Br'
        elif field in ('bt', 'Bt', 'Btheta', 'btheta'):
            field = rderavg(self.pol, self.radratio, spectral=self.rcheb)
            field = field[:, indPlot]
            bt, bp = self.sh.spec_spat(field, self.tor[:, indPlot])
            rprof = bt/self.radius[indPlot]
            if labTex:
                label = r'$B_\theta$'
            else:
                label = 'Btheta'
        elif field in ('bp', 'Bp', 'bphi', 'Bphi'):
            field = rderavg(self.pol, self.radratio, spectral=self.rcheb)
            field = field[:, indPlot]
            bt, bp = self.sh.spec_spat(field, self.tor[:, indPlot])
            rprof = bp/self.radius[indPlot]
            if labTex:
                label = r'$B_\phi$'
            else:
                label = 'Bphi'
        t2 = time.time()
        print('Transform time (surf): %.2f' % (t2-t1))

        rprof = symmetrize(rprof, self.minc)

        if field in ('temperature', 't', 'T', 'entropy', 's', 'S', 'u2', 'b2', 'nrj'):
            normed = False

        fig = radialContour(rprof, rad, label, proj, lon_0, vmax, vmin, lat_0, levels,
                            cm, normed, cbar, tit, lines)


if __name__ == '__main__':
    p = MagicPotential(field='B', ave=True)
    p.surf(field='br', r=0.9)
    p.avg(field='br')
    p.equat(field='br')

    plt.show()
