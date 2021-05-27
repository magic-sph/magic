# -*- coding: utf-8 -*-
from magic import npfile, scanDir, MagicSetup, hammer2cart, symmetrize, progressbar
import os, re
import numpy as np
import matplotlib.pyplot as plt
from .spectralTransforms import SpectralTransforms
from magic.setup import labTex
import copy


def deriv(x, y, axis=0):
    """
    This function is a simple second order derivative

    :param x: input x-axis
    :type x: numpy.ndarray
    :param y: input array
    :type y: numpy.ndarray
    :returns: an array that contains the derivatives
    :rtype: numpy.ndarray
    """
    if len(x) < 3:
        exit("Paramaters must have at least three points")
    #if len(x) != len(y):
        #exit("Vectors must have the same size")
    d = (np.roll(y, -1, axis=axis)-np.roll(y, 1, axis=axis))/ \
        (np.roll(x, -1)-np.roll(x, 1))
    d[..., 0] = (-3.*y[..., 0]+4.*y[..., 1]-y[..., 2])/(x[..., 2]-x[..., 0])
    d[..., -1] = (3*y[..., -1]-4*y[..., -2]+y[..., -3])/(x[..., -1]-x[..., -3])
    return d

def getGauss(alm, blm, ell, m, scale_b, ratio_cmb_surface, rcmb):
    """
    Get the Gauss coefficients from the real and imaginary parts
    of the poloidal potential

    :param alm: real part of the poloidal potential
    :type alm: numpy.ndarray
    :param blm: imaginary part of the poloidal potential
    :type blm: numpy.ndarray
    :param ell: spherical harmonic degree \ell
    :type ell: numpy.ndarray
    :param scale_b: magnetic field unit (default is 1)
    :type scale_b: float
    :param ratio_cmb_surface: ratio of surface ratio to CMB radius (default is 1)
    :type ratio_cmb_surface: float
    :param rcmb: radius of the outer boundary
    :type rcmb: float
    """
    fac = (-1)**m*ell*np.sqrt((2*ell+1.)/(4.*np.pi))
    fac[m > 0 ] *= np.sqrt(2)
    glm = scale_b*ratio_cmb_surface**(ell+2.)/rcmb**2*fac*alm
    hlm = -scale_b*ratio_cmb_surface**(ell+2.)/rcmb**2*fac*blm

    return glm, hlm

def rearangeLat(field):
    """
    This function is used to unfold the colatitudes

    :param field: input array with MagIC ordering of colatitudes (i.e.
                  successively Northern Hemisphere and Southern
                  Hemisphere)
    :type field: numpy.ndarray
    :return: an array with the regular ordering of the colatitudes
    :rtype: numpy.ndarray
    """
    even = field[:, ::2]
    odd = field[:, 1::2]
    return np.concatenate((even, odd[:, ::-1]), axis=1)


class MagicCoeffCmb(MagicSetup):
    """
    This class allows to read the :ref:`B_coeff_cmb.TAG <secCoeffFiles>` files.
    It first read the poloidal potential at the CMB and then transform
    it to the Gauss coefficients :math:`g_{\ell m}` and :math:`h_{\ell m}`
    using the getGauss function.

    >>> # Reads the files B_coeff_cmb.testa, B_coeff_cmb.testb
    >>> # and B_coeff_cmb.testc and stack them in one single time series
    >>> cmb = MagicCoeffCmb(tag='test[a-c]')
    >>> print(cmb.ell, cmb.glm) # print \ell and g_{\ell m}
    >>> print(cmb.glm[:, cmb.idx[1, 0]]) # time-series of the axisymmetric dipole
    >>> plot(cmb.time, cmb.dglmdt[:, cmb.idx[2, 0]]) # Secular variation of the quadrupole
    >>> # Display the time-evolution of the CMB field
    >>> cmb.movieCmb(levels=12, cm='seismic')
    >>> # Save the time-evolution of the CMB field
    >>> cmb.movieCmb(levels=12, cm='seismic', png=True)
    """

    def __init__(self, tag=None, datadir='.', ratio_cmb_surface=1, scale_b=1,
                 iplot=True, lCut=None, precision=np.float64, ave=False, sv=False,
                 quiet=False):
        """
        A class to read the B_coeff_cmb files

        :param tag: if you specify a pattern, it tries to read the corresponding files
        :type tag: str
        :param ratio_cmb_surface: ratio of surface ratio to CMB radius (default is 1)
        :type ratio_cmb_surface: float
        :param scale_b: magnetic field unit (default is 1)
        :type scale_b: float
        :param iplot: a logical to toggle the plot (default is True)
        :type iplot: int
        :param precision: single or double precision
        :type precision: char
        :param ave: load a time-averaged CMB file when set to True
        :type ave: bool
        :param sv: load a dt_b CMB file when set to True
        :type sv: bool
        :param quiet: verbose when toggled to True (default is True)
        :type quiet: bool
        :param lCut: reduce the spherical harmonic truncation to l <= lCut
        :type lCut: int
        :param datadir: working directory
        :type datadir: str
        """
        pattern = os.path.join(datadir, 'log.*')
        logFiles = scanDir(pattern)

        if ave:
            self.name = 'B_coeff_cmb_ave'
        elif sv:
            self.name = 'B_coeff_dt_cmb'
        else:
            self.name = 'B_coeff_cmb'

        if tag is not None:
            pattern = os.path.join(datadir,  '{}.{}'.format(self.name, tag))
            files = scanDir(pattern)

            # Either the log.tag directly exists and the setup is easy to obtain
            if os.path.exists(os.path.join(datadir, 'log.{}'.format(tag))):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.{}'.format(tag))
            # Or the tag is a bit more complicated and we need to find
            # the corresponding log file
            else:
                pattern = os.path.join(datadir, '{}'.format(self.name))
                mask = re.compile(r'{}\.(.*)'.format(pattern))
                if mask.match(files[-1]):
                    ending = mask.search(files[-1]).groups(0)[0]
                    pattern = os.path.join(datadir, 'log.{}'.format(ending))
                    if logFiles.__contains__(pattern):
                        MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml='log.{}'.format(ending))
        else:
            pattern = os.path.join(datadir, '{}.*'.format(self.name))
            files = scanDir(pattern)
            filename = files[-1]
            # Determine the setup
            mask = re.compile(r'{}\.(.*)'.format(self.name))
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists('log.{}'.format(ending)):
                try:
                    MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.{}'.format(ending))
                except AttributeError:
                    pass

        self.rcmb = 1./(1.-self.radratio)
        ricb = self.radratio/(1.-self.radratio)

        # Read the B_coeff files (by stacking the different tags)
        data = []
        for k, file in enumerate(files):
            if not quiet: print('Reading {}'.format(file))
            f = npfile(file, endian='B')
            self.l_max_cmb, self.minc, n_data = f.fort_read('i')
            self.m_max_cmb = int((self.l_max_cmb/self.minc)*self.minc)

            while 1:
                try:
                    data.append(f.fort_read(precision))
                except TypeError:
                    break
        self.lm_max_cmb = self.m_max_cmb*(self.l_max_cmb+1)//self.minc - \
                          self.m_max_cmb*(self.m_max_cmb-self.minc)//(2*self.minc) + \
                          self.l_max_cmb-self.m_max_cmb+1

        # Get indices location
        self.idx = np.zeros((self.l_max_cmb+1, self.m_max_cmb+1), 'i')
        self.ell = np.zeros(self.lm_max_cmb, 'i')
        self.ms = np.zeros(self.lm_max_cmb, 'i')
        self.idx[0:self.l_max_cmb+2, 0] = np.arange(self.l_max_cmb+1)
        self.ell[0:self.l_max_cmb+2] = np.arange(self.l_max_cmb+2)
        k = self.l_max_cmb+1
        for m in range(self.minc, self.l_max_cmb+1, self.minc):
            for l in range(m, self.l_max_cmb+1):
                self.idx[l, m] = k
                self.ell[self.idx[l,m]] = l
                self.ms[self.idx[l,m]] = m
                k +=1

        # Rearange data
        data = np.array(data, dtype=precision)
        self.nstep = data.shape[0]
        self.blm = np.zeros((self.nstep, self.lm_max_cmb), np.complex128)
        self.blm[:, 1:self.l_max_cmb+1] = data[:, 1:self.l_max_cmb+1]
        self.blm[:, self.l_max_cmb+1:] = data[:, self.l_max_cmb+1::2]+\
                                         1j*data[:, self.l_max_cmb+2::2]

        # Truncate!
        if lCut is not None:
            if lCut < self.l_max_cmb:
                self.truncate(lCut)

        # Get time
        self.time = np.zeros(self.nstep, precision)
        self.time = data[:, 0]

        # Get Gauss coefficients
        self.glm = np.zeros((self.nstep, self.lm_max_cmb), precision)
        self.hlm = np.zeros((self.nstep, self.lm_max_cmb), precision)

        self.glm, self.hlm = getGauss(self.blm.real, self.blm.imag,
                                      self.ell, self.ms, scale_b,
                                      ratio_cmb_surface, self.rcmb)

        # Time-averaged Gauss coefficient
        if not ave:
            facT = 1./(self.time[-1]-self.time[0])
            self.glmM = facT * np.trapz(self.glm, self.time, axis=0)
            self.hlmM = facT * np.trapz(self.hlm, self.time, axis=0)

            if len(self.time) > 3:
                self.dglmdt = deriv(self.time, self.glm.T, axis=1)
                self.dhlmdt = deriv(self.time, self.hlm.T, axis=1)
                self.dglmdt = self.dglmdt.T
                self.dhlmdt = self.dhlmdt.T

            else:
                self.dglmdt = np.zeros_like(self.glm)
                self.dhlmdt = np.zeros_like(self.hlm)

        # Magnetic energy (Lowes)
        self.El = np.zeros((self.nstep, self.l_max_cmb+1), precision)
        self.Em = np.zeros((self.nstep, self.m_max_cmb+1), precision)
        self.ESVl = np.zeros((self.nstep, self.l_max_cmb+1), precision)
        E = 0.
        for l in range(1, self.l_max_cmb+1):
            self.El[:, l] = 0.
            self.ESVl[:, l] = 0.
            for m in range(0, l+1, self.minc):
                lm = self.idx[l, m]
                self.El[:, l] += (self.ell[lm]+1)*\
                                (self.glm[:,lm]**2+self.hlm[:,lm]**2)
                self.Em[:, m] += (self.ell[lm]+1)*\
                                (self.glm[:,lm]**2+self.hlm[:,lm]**2)
                if not ave:
                    self.ESVl[:, l] += (self.ell[lm]+1)*\
                                      (self.dglmdt[:, lm]**2+self.dhlmdt[:, lm]**2)

        if not ave:
            # Time-averaged energy
            self.ElM = facT * np.trapz(self.El, self.time, axis=0)
            self.EmM = facT * np.trapz(self.Em, self.time, axis=0)

            # Secular variation
            self.ESVlM = facT * np.trapz(self.ESVl, self.time, axis=0)
            if abs(self.ESVlM[1:]).min() > 0.:
                self.taul = np.sqrt(self.ElM[1:]/self.ESVlM[1:])
            else:
                self.taul = np.zeros_like(self.ElM[1:])

        if iplot and not ave:
            self.plot()

    def __add__(self, new):
        """
        Built-in function to sum two cmb files

        .. note:: So far this function only works for two cmb files with the same
                  grid sizes. At some point, we might introduce grid extrapolation
                  to allow any summation/
        """

        out = copy.deepcopy(new)
        out.nstep = new.nstep+self.nstep
        out.time = np.concatenate((self.time, new.time), axis=0)
        out.blm = np.concatenate((self.blm, new.blm), axis=0)
        out.glm = np.concatenate((self.glm, new.glm), axis=0)
        out.hlm = np.concatenate((self.hlm, new.hlm), axis=0)
        out.El = np.concatenate((self.El, new.El), axis=0)
        out.Em = np.concatenate((self.Em, new.Em), axis=0)
        out.ESVl = np.concatenate((self.ESVl, new.ESVl), axis=0)

        return out

    def truncate(self, lCut):
        """
        :param lCut: truncate to spherical harmonic degree lCut
        :type lCut: int
        """
        self.l_max_cmb = lCut
        self.m_max_cmb = int((self.l_max_cmb/self.minc)*self.minc)
        self.lm_max_cmb = self.m_max_cmb*(self.l_max_cmb+1)/self.minc - \
                        self.m_max_cmb*(self.m_max_cmb-self.minc)/(2*self.minc) + \
                        self.l_max_cmb-self.m_max_cmb+1

        # Get indices location
        idx_new = np.zeros((self.l_max_cmb+1, self.m_max_cmb+1), 'i')
        ell_new = np.zeros(self.lm_max_cmb, 'i')
        ms_new = np.zeros(self.lm_max_cmb, 'i')
        idx_new[0:self.l_max_cmb+2, 0] = np.arange(self.l_max_cmb+1)
        ell_new[0:self.l_max_cmb+2] = np.arange(self.l_max_cmb+2)
        k = self.l_max_cmb+1
        for m in range(self.minc, self.l_max_cmb+1, self.minc):
            for l in range(m, self.l_max_cmb+1):
                idx_new[l, m] = k
                ell_new[idx_new[l,m]] = l
                ms_new[idx_new[l,m]] = m
                k +=1

        blm_new = np.zeros((self.nstep, self.lm_max_cmb), np.complex128)
        for l in range(1, self.l_max_cmb+1):
            for m in range(0, l+1, self.minc):
                lm = idx_new[l, m]
                blm_new[:, lm] = self.blm[:, self.idx[l,m]]

        self.idx = idx_new
        self.ell = ell_new
        self.ms = ms_new
        self.blm = blm_new

    def plot(self):
        """
        Display some results when iplot is set to True
        """
        ell = np.arange(self.l_max_cmb+1)
        fig = plt.figure()
        ax = fig.add_subplot(211)
        ax.semilogy(ell[1:], self.ElM[1:], marker='o')
        if labTex:
            ax.set_xlabel(r'$\ell$')
        else:
            ax.set_xlabel('Degree l')
        ax.set_ylabel('Magnetic energy')
        ax.set_xlim(1., self.l_max_cmb)

        ax1 = fig.add_subplot(212)
        ax1.semilogy(ell[0:self.m_max_cmb+1:self.minc], self.EmM[::self.minc],
                     marker='o')
        if labTex:
            ax1.set_xlabel(r'$m$')
        else:
            ax1.set_xlabel('Order m')
        ax1.set_ylabel('Magnetic energy')

        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
        ax.loglog(ell[1:], self.taul, marker='o')
        if labTex:
            ax.set_xlabel(r'$\ell$')
            ax.set_ylabel(r'$\tau_\ell$')
        else:
            ax.set_xlabel('Degree l')
            ax.set_ylabel('tau  l')
        ax.set_xlim(1, self.l_max_cmb)

        fig2 = plt.figure()
        ax = fig2.add_subplot(111)
        ax.plot(self.time, self.glm[:, self.idx[1,0]], label='g10')
        ax.plot(self.time, self.glm[:, self.idx[2,0]], label='g20')
        ax.plot(self.time, self.glm[:, self.idx[3,0]], label='g30')
        ax.set_xlabel('Time')
        ax.set_ylabel('Gauss coefficients')

    def timeLongitude(self, removeMean=True, lat0=0., levels=12, cm='RdYlBu_r',
                      deminc=True, shtns_lib='shtns'):
        """
        Plot the time-longitude diagram of Br (input latitude can be chosen)

        .. warning:: the python bindings of `SHTns <https://bitbucket.org/bputigny/shtns-magic>`_ are mandatory to use this plotting function!

        :param lat0: value of the latitude
        :type lat0: float
        :param levels: number of contour levels
        :type levels: int
        :param cm: name of the colormap
        :type cm: str
        :param deminc: a logical to indicate if one wants do get rid of the
                       possible azimuthal symmetry
        :type deminc: bool
        :param shtns_lib: version of shtns library used: can be either 'shtns'
                          or 'shtns-magic'
        :type shtns_lib: char
        :param removeMean: remove the time-averaged part when set to True
        :type removeMean: bool
        """
        # The python bindings of shtns are mandatory to use this function !!!
        import shtns

        if removeMean:
            blmCut = self.blm-self.blm.mean(axis=0)
        else:
            blmCut = self.blm

        # Define shtns setup
        sh = shtns.sht(int(self.l_max_cmb), int(self.m_max_cmb/self.minc),
                       mres=int(self.minc),
                       norm=shtns.sht_orthonormal | shtns.SHT_NO_CS_PHASE)

        polar_opt_threshold = 1e-10
        nlat = max(int(self.l_max_cmb*(3./2./2.)*2.),192)
        nphi = 2*nlat/self.minc
        nlat, nphi = sh.set_grid(nlat, nphi, polar_opt=polar_opt_threshold)

        th = np.linspace(np.pi/2., -np.pi/2., nlat)
        lat0 *= np.pi/180.
        mask = np.where(abs(th-lat0) == abs(th-lat0).min(), 1, 0)
        idx = np.nonzero(mask)[0][0]

        # Transform data on grid space
        BrCMB = np.zeros((self.nstep, nphi, nlat), np.float64)
        if deminc:
            dat = np.zeros((self.nstep, self.minc*nphi+1), np.float64)
        else:
            dat = np.zeros((self.nstep, nphi), np.float64)
        for k in range(self.nstep):
            tmp = sh.synth(blmCut[k, :]*sh.l*(sh.l+1)/self.rcmb**2)
            tmp = tmp.T # Longitude, Latitude

            if shtns_lib == 'shtns-magic':
                BrCMB[k, ...] = rearangeLat(tmp)
            else:
                BrCMB[k, ...] = tmp

            if deminc:
                dat[k, :] = symmetrize(BrCMB[k, :, idx], self.minc)
            else:
                dat[k, :] = BrCMB[k, :, idx]


        th = np.linspace(np.pi/2., -np.pi/2., nlat)
        if deminc:
            phi = np.linspace(-np.pi, np.pi, self.minc*nphi+1)
        else:
            phi = np.linspace(-np.pi/self.minc, np.pi/self.minc, nphi)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        vmin = -max(abs(dat.max()), abs(dat.min()))
        vmax = -vmin
        cs = np.linspace(vmin, vmax, levels)
        ax.contourf(phi, self.time, dat, cs, cmap=plt.get_cmap(cm))

        ax.set_xlabel('Longitude')
        ax.set_ylabel('Time')

        w2 = np.fft.fft2(dat)
        w2 = abs(w2[1:self.nstep/2+1, 0:self.m_max_cmb+1])

        dw = 2.*np.pi/(self.time[-1]-self.time[0])
        omega = dw*np.arange(self.nstep)
        omega = omega[1:self.nstep/2+1]
        ms = np.arange(self.m_max_cmb+1)

        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax1.contourf(ms, omega, w2, 17, cmap=plt.get_cmap('jet'))
        ax1.set_yscale('log')
        ax1.set_xlim(0,13)
        ax1.set_xlabel(r'Azimuthal wavenumber')
        ax1.set_ylabel(r'Frequency')

    def movieCmb(self, cut=0.5, levels=12, cm='RdYlBu_r', png=False, step=1,
                 normed=False, dpi=80, bgcolor=None, deminc=True, removeMean=False,
                 precision=np.float64, contour=False, mer=False):
        """
        Plotting function (it can also write the png files)

        :param levels: number of contour levels
        :type levels: int
        :param cm: name of the colormap
        :type cm: str
        :param cut: adjust the contour extrema to max(abs(data))*cut
        :type cut: float
        :param png: save the movie as a series of png files when
                    set to True
        :type png: bool
        :param dpi: dot per inch when saving PNGs
        :type dpi: int
        :param bgcolor: background color of the figure
        :type bgcolor: str
        :param normed: the colormap is rescaled every timestep when set to True,
                       otherwise it is calculated from the global extrema
        :type normed: bool
        :param step: the stepping between two timesteps
        :type step: int
        :param deminc: a logical to indicate if one wants do get rid of the
                       possible azimuthal symmetry
        :type deminc: bool
        :param precision: single or double precision
        :type precision: char
        :param contour: also display the solid contour levels when set to True
        :type contour: bool
        :param mer: display meridians and circles when set to True
        :type mer: bool
        :param removeMean: remove the time-averaged part when set to True
        :type removeMean: bool
        """

        if removeMean:
            blmCut = self.blm-self.blm.mean(axis=0)
        else:
            blmCut = self.blm

        nlat = int(max(int(self.l_max_cmb*(3./2./2.)*2.),192))
        if np.mod(nlat, 2) == 1:
            nlat += 1
        nphi = int(2*nlat/self.minc)

        # Define spectral transform setup
        sh = SpectralTransforms(l_max=self.l_max_cmb, minc=self.minc,
                                lm_max=self.lm_max_cmb,
                                n_theta_max=nlat)

        # Transform data on grid space
        BrCMB = np.zeros((self.nstep, nphi, nlat), precision)
        print('Spectral -> Spatial transform')
        for k in progressbar(range(self.nstep)):
            BrCMB[k, ...] = sh.spec_spat(blmCut[k, :]*self.ell*(sh.ell+1)/self.rcmb**2)
        print('Done')

        if png:
            plt.ioff()
            if not os.path.exists('movie'):
                os.mkdir('movie')
        else:
            plt.ion()

        if not normed:
            vmin = - max(abs(BrCMB.max()), abs(BrCMB.min()))
            vmin = cut * vmin
            vmax = -vmin
            cs = np.linspace(vmin, vmax, levels)

        th = np.linspace(np.pi/2., -np.pi/2., nlat)
        if deminc:
            phi = np.linspace(-np.pi, np.pi, self.minc*nphi+1)
            xxout, yyout = hammer2cart(th, -np.pi)
            xxin, yyin = hammer2cart(th, np.pi)
        else:
            phi = np.linspace(-np.pi/self.minc, np.pi/self.minc, nphi)
            xxout, yyout = hammer2cart(th, -np.pi/self.minc)
            xxin, yyin = hammer2cart(th, np.pi/self.minc)
        ttheta, pphi = np.meshgrid(th, phi)
        xx, yy = hammer2cart(ttheta, pphi)
        if deminc:
            fig = plt.figure(figsize=(8, 4))
        else:
            fig = plt.figure(figsize=(8/self.minc, 4))
        fig.subplots_adjust(top=0.99, right=0.99, bottom=0.01, left=0.01)
        ax = fig.add_subplot(111, frameon=False)

        if mer:
            theta = np.linspace(np.pi/2, -np.pi/2, nlat)
            meridians = np.r_[-120, -60, 0, 60, 120]
            circles = np.r_[ 60, 30, 0, -30, -60]

        for k in range(self.nstep):
            if k == 0:
                if normed:
                    vmin = - max(abs(BrCMB[k, ...].max()), abs(BrCMB[k, ...].min()))
                    vmin = cut * vmin
                    vmax = -vmin
                    cs = np.linspace(vmin, vmax, levels)
                if deminc:
                    dat = symmetrize(BrCMB[k, ...], self.minc)
                else:
                    dat = BrCMB[k, ...]
                im = ax.contourf(xx, yy, dat, cs, cmap=plt.get_cmap(cm), extend='both')
                if contour:
                    ax.contour(xx, yy, dat, cs, linestyles=['-', '-'],
                               colors=['k', 'k'], linewidths=[0.7, 0.7])
                ax.plot(xxout, yyout, 'k-', lw=1.5)
                ax.plot(xxin, yyin, 'k-', lw=1.5)
                #ax.text(0.12, 0.9, 't={:.6f}'.format(self.time[0]), fontsize=16,
                        #horizontalalignment='right',
                        #verticalalignment='center', transform = ax.transAxes)

                if mer:
                    for lat0 in circles:
                        x0, y0 = hammer2cart(lat0*np.pi/180., phi)
                        ax.plot(x0, y0, 'k:', linewidth=0.7)
                    for lon0 in meridians:
                        x0, y0 = hammer2cart(theta, lon0*np.pi/180.)
                        ax.plot(x0, y0, 'k:', linewidth=0.7)

                ax.axis('off')
                man = plt.get_current_fig_manager()
                man.canvas.draw()
            if k != 0 and k % step == 0:
                if not png:
                    print(k)
                plt.cla()
                if normed:
                    vmin = - max(abs(BrCMB[k, ...].max()), abs(BrCMB[k, ...].min()))
                    vmin = cut * vmin
                    vmax = -vmin
                    cs = np.linspace(vmin, vmax, levels)
                if deminc:
                    dat = symmetrize(BrCMB[k, ...], self.minc)
                else:
                    dat = BrCMB[k, ...]
                im = ax.contourf(xx, yy, dat, cs, cmap=plt.get_cmap(cm), extend='both')
                if contour:
                    ax.contour(xx, yy, dat, cs, colors=['k'],
                               linestyles=['-', '-'], linewidths=[0.7, 0.7])
                ax.plot(xxout, yyout, 'k-', lw=1.5)
                ax.plot(xxin, yyin, 'k-', lw=1.5)
                #ax.text(0.12, 0.9, 't={:.6f}'.format(self.time[k]), fontsize=16,
                        #horizontalalignment='right',
                        #verticalalignment='center', transform = ax.transAxes)

                if mer:
                    for lat0 in circles:
                        x0, y0 = hammer2cart(lat0*np.pi/180., phi)
                        ax.plot(x0, y0, 'k:', linewidth=0.7)
                    for lon0 in meridians:
                        x0, y0 = hammer2cart(theta, lon0*np.pi/180.)
                        ax.plot(x0, y0, 'k:', linewidth=0.7)

                ax.axis('off')
                man.canvas.draw()
            if png:
                filename = 'movie/img{:05d}.png'.format(k)
                print('write {}'.format(filename))
                #st = 'echo {}'.format(ivar) + ' > movie/imgmax'
                if bgcolor is not None:
                    fig.savefig(filename, facecolor=bgcolor, dpi=dpi)
                else:
                    fig.savefig(filename, dpi=dpi)



class MagicCoeffR(MagicSetup):
    """
    This class allows to read the :ref:`B_coeff_r#.TAG <secBcoeffrFile>`
    and :ref:`V_coeff_r#.TAG <secVcoeffrFile>` files.
    It reads the poloidal and toroidal potentials and reconstruct the time
    series (or the energy) contained in any given mode.

    >>> # Reads the files V_coeff_r2.test*
    >>> cr = MagicCoeffR(tag='test*', field='V', r=2)
    >>> print(cr.ell, cr.wlm) # print \ell and w_{\ell m}
    >>> # Time-evolution of the poloidal energy in the (\ell=10, m=10) mode
    >>> plot(cr.time, cr.epolLM[:, cr.idx[10, 10]])
    """

    def __init__(self, tag, ratio_cmb_surface=1, scale_b=1, iplot=True,
                 field='B', r=1, precision=np.float64, lCut=None, quiet=False):
        """
        :param tag: if you specify a pattern, it tries to read the corresponding files
        :type tag: str
        :param ratio_cmb_surface: ratio of surface ratio to CMB radius (default is 1)
        :type ratio_cmb_surface: float
        :param scale_b: magnetic field unit (default is 1)
        :type scale_b: float
        :param iplot: a logical to toggle the plot (default is True)
        :type iplot: bool
        :param field: 'B', 'V', 'T' or 'Xi' (magnetic field, velocity field, temperature or composition)
        :type field: str
        :param r: an integer to characterise which file we want to plot
        :type r: int
        :param precision: single or double precision
        :type precision: str
        :param lCut: reduce the spherical harmonic truncation to l <= lCut
        :type lCut: int
        :param quiet: verbose when toggled to True (default is True)
        :type quiet: bool
        """

        logFiles = scanDir('log.*')
        if len(logFiles) != 0:
            MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])
        else:
            str1 = 'Aspect ratio ?\n'
            self.radratio = float(input(str1))

        self.rcmb = 1./(1.-self.radratio)
        ricb = self.radratio/(1.-self.radratio)

        files = scanDir('{}_coeff_r{}.{}'.format(field,r,tag))

        # Read the B_coeff files (by stacking the different tags)
        data = []
        for k, file in enumerate(files):
            if not quiet: print('Reading {}'.format(file))
            f = npfile(file, endian='B')
            if precision == np.float32:
                out = f.fort_read('3i4,f4')[0]
            else:
                out = f.fort_read('3i4,f8')[0]
            self.l_max_r, self.minc, n_data = out[0]
            self.m_max_r = int((self.l_max_r//self.minc)*self.minc)
            self.radius = out[1]

            while 1:
                try:
                    data.append(f.fort_read(precision))
                except TypeError:
                    break
        self.lm_max_r = self.m_max_r*(self.l_max_r+1)//self.minc - \
                        self.m_max_r*(self.m_max_r-self.minc)//(2*self.minc) + \
                        self.l_max_r-self.m_max_r+1

        # Get indices location
        self.idx = np.zeros((self.l_max_r+1, self.m_max_r+1), 'i')
        self.ell = np.zeros(self.lm_max_r, 'i')
        self.ms = np.zeros(self.lm_max_r, 'i')
        self.idx[0:self.l_max_r+2, 0] = np.arange(self.l_max_r+1)
        self.ell[0:self.l_max_r+2] = np.arange(self.l_max_r+2)
        k = self.l_max_r+1
        for m in range(self.minc, self.l_max_r+1, self.minc):
            for l in range(m, self.l_max_r+1):
                self.idx[l, m] = k
                self.ell[self.idx[l,m]] = l
                self.ms[self.idx[l,m]] = m
                k +=1


        # Rearange data
        data = np.array(data, dtype=precision)
        self.nstep = data.shape[0]
        self.wlm = np.zeros((self.nstep, self.lm_max_r), np.complex128)
        if field == 'V' or field == 'B':
            self.dwlm = np.zeros((self.nstep, self.lm_max_r), np.complex128)
            self.zlm = np.zeros((self.nstep, self.lm_max_r), np.complex128)

        # Get time
        self.time = np.zeros(self.nstep, dtype=precision)
        self.time = data[:, 0]

        # wlm
        if field == 'T' or field == 'Xi': #T or Xi contains l = m = 0
            self.wlm[:, 0:self.l_max_r+1] = data[:, 1:self.l_max_r+2]
        else:
            self.wlm[:, 1:self.l_max_r+1] = data[:, 1:self.l_max_r+1]
        k = self.l_max_r+1
        for m in range(self.minc, self.l_max_r+1, self.minc):
            for l in range(m, self.l_max_r+1):
                self.wlm[:, self.idx[l, m]] = data[:, k]+1j*data[:, k+1]
                k += 2

        if field == 'V' or field == 'B':
            # dwlm
            self.dwlm[:, 1:self.l_max_r+1] = data[:, k:k+self.l_max_r]
            k += self.l_max_r
            for m in range(self.minc, self.l_max_r+1, self.minc):
                for l in range(m, self.l_max_r+1):
                    self.dwlm[:, self.idx[l, m]] = data[:, k]+1j*data[:, k+1]
                    k += 2
            # zlm
            self.zlm[:, 1:self.l_max_r+1] = data[:, k:k+self.l_max_r]
            k += self.l_max_r
            for m in range(self.minc, self.l_max_r+1, self.minc):
                for l in range(m, self.l_max_r+1):
                    self.zlm[:, self.idx[l, m]] = data[:, k]+1j*data[:, k+1]
                    k += 2

        # ddw in case B is stored
        if field == 'B':
            self.ddwlm = np.zeros((self.nstep, self.lm_max_r), np.complex128)
            self.ddwlm[:, 1:self.l_max_r+1] = data[:, k:k+self.l_max_r]
            k += self.l_max_r
            for m in range(self.minc, self.l_max_r+1, self.minc):
                for l in range(m, self.l_max_r+1):
                    self.ddwlm[:, self.idx[l, m]] = data[:, k]+1j*data[:, k+1]
                    k += 2

        # Truncate!
        if lCut is not None:
            if lCut < self.l_max_r:
                self.truncate(lCut, field=field)

        if field == 'V' or field == 'B':
            self.e_pol_axi_l = np.zeros((self.nstep, self.l_max_r+1), precision)
            self.e_tor_axi_l = np.zeros((self.nstep, self.l_max_r+1), precision)
            self.e_pol_l = np.zeros((self.nstep, self.l_max_r+1), precision)
            self.e_tor_l = np.zeros((self.nstep, self.l_max_r+1), precision)

            for l in range(1, self.l_max_r+1):
                self.e_pol_l[:, l] = 0.
                self.e_tor_l[:, l] = 0.
                self.e_pol_axi_l[:, l] = 0.
                self.e_tor_axi_l[:, l] = 0.
                for m in range(0, l+1, self.minc):
                    lm = self.idx[l, m]

                    if m == 0:
                        epol = 0.5*self.ell[lm]*(self.ell[lm]+1)*( \
                               self.ell[lm]*(self.ell[lm]+1)/self.radius**2* \
                               abs(self.wlm[:,lm])**2+ abs(self.dwlm[:,lm])**2 )
                        etor = 0.5*self.ell[lm]*(self.ell[lm]+1)*abs(self.zlm[:, lm])**2

                        self.e_pol_axi_l[:, l] += epol
                        self.e_tor_axi_l[:, l] += etor
                    else:
                        epol = self.ell[lm]*(self.ell[lm]+1)*( \
                               self.ell[lm]*(self.ell[lm]+1)/self.radius**2* \
                               abs(self.wlm[:,lm])**2+ abs(self.dwlm[:,lm])**2 )
                        etor = self.ell[lm]*(self.ell[lm]+1)*abs(self.zlm[:, lm])**2

                    self.e_pol_l[:, l] += epol
                    self.e_tor_l[:, l] += etor


            # Time-averaged energy
            facT = 1./(self.time[-1]-self.time[0])

            self.e_pol_lM = facT * np.trapz(self.e_pol_l, self.time, axis=0)
            self.e_tor_lM = facT * np.trapz(self.e_tor_l, self.time, axis=0)
            self.e_pol_axi_lM = facT * np.trapz(self.e_pol_axi_l, self.time, axis=0)
            self.e_tor_axi_lM = facT * np.trapz(self.e_tor_axi_l, self.time, axis=0)

    def truncate(self, lCut, field='B'):
        """
        :param lCut: truncate to spherical harmonic degree lCut
        :type lCut: int
        :param field: name of the field ('V', 'B' or 'T')
        :type field: char
        """
        self.l_max_r = lCut
        self.m_max_r = int((self.l_max_r/self.minc)*self.minc)

        self.lm_max_r = int(self.m_max_r*(self.l_max_r+1)/self.minc - \
                        self.m_max_r*(self.m_max_r-self.minc)/(2*self.minc) + \
                        self.l_max_r-self.m_max_r+1)

        # Get indices location
        idx_new = np.zeros((self.l_max_r+1, self.m_max_r+1), 'i')
        ell_new = np.zeros(self.lm_max_r, 'i')
        ms_new = np.zeros(self.lm_max_r, 'i')
        idx_new[0:self.l_max_r+2, 0] = np.arange(self.l_max_r+1)
        ell_new[0:self.l_max_r+2] = np.arange(self.l_max_r+2)
        k = self.l_max_r+1
        for m in range(self.minc, self.l_max_r+1, self.minc):
            for l in range(m, self.l_max_r+1):
                idx_new[l, m] = k
                ell_new[idx_new[l,m]] = l
                ms_new[idx_new[l,m]] = m
                k +=1

        wlm_new = np.zeros((self.nstep, self.lm_max_r), np.complex128)
        if field == 'V' or field == 'B':
            dwlm_new = np.zeros((self.nstep, self.lm_max_r), np.complex128)
            zlm_new = np.zeros((self.nstep, self.lm_max_r), np.complex128)
            if field == 'B':
                ddwlm_new = np.zeros((self.nstep, self.lm_max_r), np.complex128)

            for l in range(1, self.l_max_r+1):
                for m in range(0, l+1, self.minc):
                    lm = idx_new[l, m]
                    wlm_new[:, lm] = self.wlm[:, self.idx[l,m]]
                    zlm_new[:, lm] = self.zlm[:, self.idx[l,m]]
                    dwlm_new[:, lm] = self.dwlm[:, self.idx[l,m]]
                    if field == 'B':
                        ddwlm_new[:, lm] = self.ddwlm[:, self.idx[l,m]]
        else:
            for l in range(1, self.l_max_r+1):
                for m in range(0, l+1, self.minc):
                    lm = idx_new[l, m]
                    wlm_new[:, lm] = self.wlm[:, self.idx[l,m]]

        #for m in range(self.minc, self.l_max_r+1, self.minc):
            #for l in range(m, self.l_max_r+1):
                #wlm_new[:, idx_new[l, m]] = self.wlm[:, self.idx[l,m]]
                #dwlm_new[:, idx_new[l, m]] = self.dwlm[:, self.idx[l,m]]
                #zlm_new[:, idx_new[l, m]] = self.zlm[:, self.idx[l,m]]
                #if field == 'B':
                    #ddwlm_new[:, idx_new[l, m]] = self.ddwlm[:, self.idx[l,m]]
        self.idx = idx_new
        self.ell = ell_new
        self.ms = ms_new

        self.wlm = wlm_new
        if field == 'V' or field == 'B':
            self.dwlm = dwlm_new
            self.zlm = zlm_new
            if field == 'B':
                self.ddwlm = ddwlm_new

    def movieRad(self, cut=0.5, levels=12, cm='RdYlBu_r', png=False, step=1,
                 normed=False, dpi=80, bgcolor=None, deminc=True, removeMean=False,
                 precision=np.float64, contour=False, mer=False):
        """
        Plotting function (it can also write the png files)

        .. warning:: the python bindings of `SHTns <https://bitbucket.org/bputigny/shtns-magic>`_ are mandatory to use this plotting function!

        :param levels: number of contour levels
        :type levels: int
        :param cm: name of the colormap
        :type cm: str
        :param cut: adjust the contour extrema to max(abs(data))*cut
        :type cut: float
        :param png: save the movie as a series of png files when
                    set to True
        :type png: bool
        :param dpi: dot per inch when saving PNGs
        :type dpi: int
        :param bgcolor: background color of the figure
        :type bgcolor: str
        :param normed: the colormap is rescaled every timestep when set to True,
                       otherwise it is calculated from the global extrema
        :type normed: bool
        :param step: the stepping between two timesteps
        :type step: int
        :param deminc: a logical to indicate if one wants do get rid of the
                       possible azimuthal symmetry
        :type deminc: bool
        :param precision: single or double precision
        :type precision: char
        :param contour: also display the solid contour levels when set to True
        :type contour: bool
        :param mer: display meridians and circles when set to True
        :type mer: bool
        :param removeMean: remove the time-averaged part when set to True
        :type removeMean: bool
        """

        if removeMean:
            dataCut = self.wlm-self.wlm.mean(axis=0)
        else:
            dataCut = self.wlm

        nlat = max(int(self.l_max_r*(3./2./2.)*2.),192)
        nphi = 2*nlat/self.minc

        # Define spectral transform setup
        sh = SpectralTransforms(l_max=self.l_max_r, minc=self.minc,
                                lm_max=self.lm_max_r,
                                n_theta_max=nlat)

        """
        # The python bindings of shtns are mandatory to use this function !!!
        import shtns

        # Define shtns setup
        sh = shtns.sht(int(self.l_max_r), int(self.m_max_r/self.minc),
                       mres=int(self.minc),
                       norm=shtns.sht_orthonormal | shtns.SHT_NO_CS_PHASE)
        """

        # Transform data on grid space
        data = np.zeros((self.nstep, nphi, nlat), precision)
        print('Spectral -> Spatial transform')
        for k in progressbar(range(self.nstep)):
            data[k, ...] = sh.spec_spat(dataCut[k, :]*self.ell*(self.ell+1)/self.radius**2)
        print('Done')

        if png:
            plt.ioff()
            if not os.path.exists('movie'):
                os.mkdir('movie')
        else:
            plt.ion()

        if not normed:
            vmin = - max(abs(data.max()), abs(data.min()))
            vmin = cut * vmin
            vmax = -vmin
            cs = np.linspace(vmin, vmax, levels)

        th = np.linspace(np.pi/2., -np.pi/2., nlat)
        if deminc:
            phi = np.linspace(-np.pi, np.pi, self.minc*nphi+1)
            xxout, yyout = hammer2cart(th, -np.pi)
            xxin, yyin = hammer2cart(th, np.pi)
        else:
            phi = np.linspace(-np.pi/self.minc, np.pi/self.minc, nphi)
            xxout, yyout = hammer2cart(th, -np.pi/self.minc)
            xxin, yyin = hammer2cart(th, np.pi/self.minc)
        ttheta, pphi = np.meshgrid(th, phi)
        xx, yy = hammer2cart(ttheta, pphi)
        if deminc:
            fig = plt.figure(figsize=(8, 4))
        else:
            fig = plt.figure(figsize=(8/self.minc, 4))
        fig.subplots_adjust(top=0.99, right=0.99, bottom=0.01, left=0.01)
        ax = fig.add_subplot(111, frameon=False)

        if mer:
            theta = np.linspace(np.pi/2, -np.pi/2, nlat)
            meridians = np.r_[-120, -60, 0, 60, 120]
            circles = np.r_[ 60, 30, 0, -30, -60]

        for k in range(self.nstep):
            if k == 0:
                if normed:
                    vmin = - max(abs(data[k, ...].max()), abs(data[k, ...].min()))
                    vmin = cut * vmin
                    vmax = -vmin
                    cs = np.linspace(vmin, vmax, levels)
                if deminc:
                    dat = symmetrize(data[k, ...], self.minc)
                else:
                    dat = data[k, ...]
                im = ax.contourf(xx, yy, dat, cs, cmap=plt.get_cmap(cm), extend='both')
                if contour:
                    ax.contour(xx, yy, dat, cs, linestyles=['-', '-'],
                               colors=['k', 'k'], linewidths=[0.7, 0.7])
                ax.plot(xxout, yyout, 'k-', lw=1.5)
                ax.plot(xxin, yyin, 'k-', lw=1.5)

                if mer:
                    for lat0 in circles:
                        x0, y0 = hammer2cart(lat0*np.pi/180., phi)
                        ax.plot(x0, y0, 'k:', linewidth=0.7)
                    for lon0 in meridians:
                        x0, y0 = hammer2cart(theta, lon0*np.pi/180.)
                        ax.plot(x0, y0, 'k:', linewidth=0.7)

                ax.axis('off')
                man = plt.get_current_fig_manager()
                man.canvas.draw()

                if png:
                    filename = 'movie/img{:05d}.png'.format(k)
                    print('write {}'.format(filename))
                    #st = 'echo {}'.format(ivar) + ' > movie/imgmax'
                    if bgcolor is not None:
                        fig.savefig(filename, facecolor=bgcolor, dpi=dpi)
                    else:
                        fig.savefig(filename, dpi=dpi)

            elif k != 0 and k % step == 0:
                if not png:
                    print(k)
                plt.cla()
                if normed:
                    vmin = - max(abs(data[k, ...].max()), abs(data[k, ...].min()))
                    vmin = cut * vmin
                    vmax = -vmin
                    cs = np.linspace(vmin, vmax, levels)
                if deminc:
                    dat = symmetrize(data[k, ...], self.minc)
                else:
                    dat = data[k, ...]
                im = ax.contourf(xx, yy, dat, cs, cmap=plt.get_cmap(cm), extend='both')
                if contour:
                    ax.contour(xx, yy, dat, cs, colors=['k'],
                               linestyles=['-', '-'], linewidths=[0.7, 0.7])
                ax.plot(xxout, yyout, 'k-', lw=1.5)
                ax.plot(xxin, yyin, 'k-', lw=1.5)

                if mer:
                    for lat0 in circles:
                        x0, y0 = hammer2cart(lat0*np.pi/180., phi)
                        ax.plot(x0, y0, 'k:', linewidth=0.7)
                    for lon0 in meridians:
                        x0, y0 = hammer2cart(theta, lon0*np.pi/180.)
                        ax.plot(x0, y0, 'k:', linewidth=0.7)

                ax.axis('off')
                man.canvas.draw()
                if png:
                    filename = 'movie/img{:05d}.png'.format(k)
                    print('write {}'.format(filename))
                    #st = 'echo {}'.format(ivar) + ' > movie/imgmax'
                    if bgcolor is not None:
                        fig.savefig(filename, facecolor=bgcolor, dpi=dpi)
                    else:
                        fig.savefig(filename, dpi=dpi)

    def fft(self):
        """
        Fourier transform of the poloidal energy
        """
        w2 = np.fft.fft(self.e_pol_l, axis=0)
        w2 = abs(w2[1:self.nstep/2+1,1:])
        dw = 2.*np.pi/(self.time[-1]-self.time[0])
        omega = dw*np.arange(self.nstep)
        omega = omega[1:self.nstep/2+1]
        ls = np.arange(self.l_max_r+1)
        ls = ls[1:]

        dat = np.log10(w2)
        vmax = dat.max()-1
        vmin = dat.min()+2
        levs = np.linspace(vmin, vmax, 65)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        im = ax.contourf(ls, omega, np.log10(w2), levs, cmap=plt.get_cmap('jet'),
                         extend='both')

        cbar = fig.colorbar(im)

        ax.set_xlabel(r'Spherical harmonic degree')
        ax.set_ylabel(r'Frequency')
