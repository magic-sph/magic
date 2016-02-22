# -*- coding: utf-8 -*-
from magic import npfile, scanDir, MagicSetup, hammer2cart, symmetrize
import os
import numpy as N
import matplotlib.pyplot as P
from magic.setup import labTex


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
    d = (N.roll(y, -1, axis=axis)-N.roll(y, 1, axis=axis))/ \
        (N.roll(x, -1)-N.roll(x, 1))
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
    fac = (-1)**m*ell*N.sqrt((2*ell+1.)/(4.*N.pi))
    fac[m > 0 ] *= N.sqrt(2)
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
    return N.concatenate((even, odd[:, ::-1]), axis=1)


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
    >>> # Display the time-evolution of the CMB field (requires shtns)
    >>> cmb.movieCmb(levels=12, cm='seismic')
    >>> # Save the time-evolution of the CMB field (requires shtns)
    >>> cmb.movieCmb(levels=12, cm='seismic', png=True)
    """
    
    def __init__(self, tag, ratio_cmb_surface=1, scale_b=1, iplot=True,
                 precision='Float64', ave=False, sv=False):
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
        """

        logFiles = scanDir('log.*')
        if len(logFiles) != 0:
            MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])
        else:
            str1 = 'Aspect ratio ?\n'
            self.radratio = float(input(str1))

        self.rcmb = 1./(1.-self.radratio)
        ricb = self.radratio/(1.-self.radratio)

        if ave:
            files = scanDir('B_coeff_cmb_ave.%s' % tag)
        elif sv:
            files = scanDir('B_coeff_dt_cmb.%s' % tag)
        else:
            files = scanDir('B_coeff_cmb.%s' % tag)

        # Read the B_coeff files (by stacking the different tags)
        data = []
        for k, file in enumerate(files):
            print('Reading %s' % file)
            f = npfile(file, endian='B')
            self.l_max_cmb, self.minc, n_data = f.fort_read('i')
            self.m_max_cmb = int((self.l_max_cmb/self.minc)*self.minc)

            while 1:
                try:
                    data.append(f.fort_read(precision))
                except TypeError:
                    break
        self.lm_max_cmb = self.m_max_cmb*(self.l_max_cmb+1)/self.minc - \
                          self.m_max_cmb*(self.m_max_cmb-self.minc)/(2*self.minc) + \
                          self.l_max_cmb-self.m_max_cmb+1

        # Get indices location
        self.idx = N.zeros((self.l_max_cmb+1, self.m_max_cmb+1), 'i')
        self.ell = N.zeros(self.lm_max_cmb, 'i')
        self.ms = N.zeros(self.lm_max_cmb, 'i')
        self.idx[0:self.l_max_cmb+2, 0] = N.arange(self.l_max_cmb+1)
        self.ell[0:self.l_max_cmb+2] = N.arange(self.l_max_cmb+2)
        k = self.l_max_cmb+1
        for m in range(self.minc, self.l_max_cmb+1, self.minc):
            for l in range(m, self.l_max_cmb+1):
                self.idx[l, m] = k
                self.ell[self.idx[l,m]] = l
                self.ms[self.idx[l,m]] = m
                k +=1

        # Rearange data
        data = N.array(data, dtype=precision)
        self.nstep = data.shape[0]
        self.blm = N.zeros((self.nstep, self.lm_max_cmb), 'Complex64')
        self.blm[:, 1:self.l_max_cmb+1] = data[:, 1:self.l_max_cmb+1]
        self.blm[:, self.l_max_cmb+1:] = data[:, self.l_max_cmb+1::2]+\
                                         1j*data[:, self.l_max_cmb+2::2]

        # Get time
        self.time = N.zeros(self.nstep, precision)
        self.time = data[:, 0]

        # Get Gauss coefficients
        self.glm = N.zeros((self.nstep, self.lm_max_cmb), precision)
        self.hlm = N.zeros((self.nstep, self.lm_max_cmb), precision)

        self.glm, self.hlm = getGauss(self.blm.real, self.blm.imag, 
                                      self.ell, self.ms, scale_b, 
                                      ratio_cmb_surface, self.rcmb)

        # Time-averaged Gauss coefficient
        facT = 1./(self.time[-1]-self.time[0])
        self.glmM = facT * N.trapz(self.glm, self.time, axis=0)
        self.hlmM = facT * N.trapz(self.hlm, self.time, axis=0)

        if len(self.time) > 3:
            self.dglmdt = deriv(self.time, self.glm.T, axis=1)
            self.dhlmdt = deriv(self.time, self.hlm.T, axis=1)
            self.dglmdt = self.dglmdt.T
            self.dhlmdt = self.dhlmdt.T

        else:
            self.dglmdt = N.zeros_like(self.glm)
            self.dhlmdt = N.zeros_like(self.hlm)

        # Magnetic energy (Lowes)
        self.El = N.zeros((self.nstep, self.l_max_cmb+1), precision)
        self.Em = N.zeros((self.nstep, self.m_max_cmb+1), precision)
        self.ESVl = N.zeros((self.nstep, self.l_max_cmb+1), precision)
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
                self.ESVl[:, l] += (self.ell[lm]+1)*\
                                  (self.dglmdt[:, lm]**2+self.dhlmdt[:, lm]**2)

        # Time-averaged energy
        self.ElM = facT * N.trapz(self.El, self.time, axis=0)
        self.EmM = facT * N.trapz(self.Em, self.time, axis=0)

        # Secular variation
        self.ESVlM = facT * N.trapz(self.ESVl, self.time, axis=0)
        self.taul = N.sqrt(self.ElM[1:]/self.ESVlM[1:])

        if iplot:
            self.plot()


    def plot(self):
        """
        Display some results when iplot is set to True
        """
        ell = N.arange(self.l_max_cmb+1)
        fig = P.figure()
        ax = fig.add_subplot(211)
        ax.semilogy(ell[1:], self.ElM[1:], 'b-o')
        if labTex:
            ax.set_xlabel(r'$\ell$')
        else:
            ax.set_xlabel('Degree l')
        ax.set_ylabel('Magnetic energy')
        ax.set_xlim(1., self.l_max_cmb)

        ax1 = fig.add_subplot(212)
        ax1.semilogy(ell[0:self.m_max_cmb+1:self.minc], self.EmM[::self.minc], 
                     'b-o')   
        if labTex:
            ax1.set_xlabel(r'$m$')
        else:
            ax1.set_xlabel('Order m')
        ax1.set_ylabel('Magnetic energy')

        fig1 = P.figure()
        ax = fig1.add_subplot(111)
        ax.loglog(ell[1:], self.taul, 'b-o')
        if labTex:
            ax.set_xlabel(r'$\ell$')
            ax.set_ylabel(r'$\tau_\ell$')
        else:
            ax.set_xlabel('Degree l')
            ax.set_ylabel('tau  l')
        ax.set_xlim(1, self.l_max_cmb)

        fig2 = P.figure()
        ax = fig2.add_subplot(111)
        ax.plot(self.time, self.glm[:, self.idx[1,0]], label='g10')
        ax.plot(self.time, self.glm[:, self.idx[2,0]], label='g20')
        ax.plot(self.time, self.glm[:, self.idx[3,0]], label='g30')
        ax.set_xlabel('Time')
        ax.set_ylabel('Gauss coefficients')

    def movieCmb(self, cut=0.5, levels=12, cmap='RdYlBu_r', png=False, step=1,
                 normed=False, dpi=80, bgcolor=None, deminc=True, 
                 precision='Float64', shtns_lib='shtns'):
        """
        Plotting function (it can also write the png files)

        .. warning:: the python bindings of `SHTns <https://bitbucket.org/bputigny/shtns-magic>`_ are mandatory to use this plotting function!

        :param levels: number of contour levels
        :type levels: int
        :param cmap: name of the colormap
        :type cmap: str
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
        :param shtns_lib: version of shtns library used: can be either 'shtns'
                          or 'shtns-magic'
        :type shtns_lib: char
        """

        # The python bindings of shtns are mandatory to use this function !!!
        import shtns

        # Define shtns setup
        sh = shtns.sht(int(self.l_max_cmb), int(self.m_max_cmb/self.minc), 
                       mres=int(self.minc), 
                       norm=shtns.sht_orthonormal | shtns.SHT_NO_CS_PHASE)

        polar_opt_threshold = 1e-10
        nlat = max((self.l_max_cmb*(3/2/2)),192)
        nphi = 2*nlat/self.minc
        nlat, nphi = sh.set_grid(nlat, nphi, polar_opt=polar_opt_threshold)

        # Transform data on grid space
        BrCMB = N.zeros((self.nstep, nphi, nlat), precision)
        for k in range(self.nstep):
            tmp = sh.synth(self.blm[k, :]*sh.l*(sh.l+1)/self.rcmb**2)
            tmp = tmp.T # Longitude, Latitude

            if shtns_lib == 'shtns-magic':
                BrCMB[k, ...] = rearangeLat
            BrCMB[k, ...] = tmp

        if png:
            P.ioff()
            if not os.path.exists('movie'):
                os.mkdir('movie')
        else:
            P.ion()

        if not normed:
            vmin = - max(abs(BrCMB.max()), abs(BrCMB.min()))
            vmin = cut * vmin
            vmax = -vmin
            cs = N.linspace(vmin, vmax, levels)

        th = N.linspace(N.pi/2., -N.pi/2., nlat)
        if deminc:
            phi = N.linspace(-N.pi, N.pi, self.minc*nphi+1)
            xxout, yyout = hammer2cart(th, -N.pi)
            xxin, yyin = hammer2cart(th, N.pi)
        else:
            phi = N.linspace(-N.pi/self.minc, N.pi/self.minc, nphi)
            xxout, yyout = hammer2cart(th, -N.pi/self.minc)
            xxin, yyin = hammer2cart(th, N.pi/self.minc)
        ttheta, pphi = N.meshgrid(th, phi)
        xx, yy = hammer2cart(ttheta, pphi)
        if deminc:
            fig = P.figure(figsize=(8, 4))
        else:
            fig = P.figure(figsize=(8/self.minc, 4))
        fig.subplots_adjust(top=0.99, right=0.99, bottom=0.01, left=0.01)
        ax = fig.add_subplot(111, frameon=False)


        for k in range(self.nstep):
            if k == 0:
                if normed:
                    vmin = - max(abs(BrCMB[k, ...].max()), abs(BrCMB[k, ...].min()))
                    vmin = cut * vmin
                    vmax = -vmin
                    cs = N.linspace(vmin, vmax, levels)
                if deminc:
                    dat = symmetrize(BrCMB[k, ...], self.minc)
                else:
                    dat = BrCMB[k, ...]
                im = ax.contourf(xx, yy, dat, cs, cmap=cmap, extend='both')
                ax.plot(xxout, yyout, 'k-', lw=1.5)
                ax.plot(xxin, yyin, 'k-', lw=1.5)
                ax.axis('off')
                man = P.get_current_fig_manager()
                man.canvas.draw()
            if k != 0 and k % step == 0:
                if not png:
                    print(k)
                P.cla()
                if normed:
                    vmin = - max(abs(BrCMB[k, ...].max()), abs(BrCMB[k, ...].min()))
                    vmin = cut * vmin
                    vmax = -vmin
                    cs = N.linspace(vmin, vmax, levels)
                if deminc:
                    dat = symmetrize(BrCMB[k, ...], self.minc)
                else:
                    dat = BrCMB[k, ...]
                im = ax.contourf(xx, yy, dat, cs, cmap=cmap, extend='both')
                ax.plot(xxout, yyout, 'k-', lw=1.5)
                ax.plot(xxin, yyin, 'k-', lw=1.5)
                ax.axis('off')
                man.canvas.draw()
            if png:
                filename = 'movie/img%05d.png' % k
                print('write %s' % filename)
                #st = 'echo %i' % ivar + ' > movie/imgmax'
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
    >>> plot(cr.time, cr.epolLM[:, 10, 10]) 
    """
    
    def __init__(self, tag, ratio_cmb_surface=1, scale_b=1, iplot=True,
                 field='B', r=1, precision='Float64'):
        """
        :param tag: if you specify a pattern, it tries to read the corresponding files
        :type tag: str
        :param ratio_cmb_surface: ratio of surface ratio to CMB radius (default is 1)
        :type ratio_cmb_surface: float
        :param scale_b: magnetic field unit (default is 1)
        :type scale_b: float
        :param iplot: a logical to toggle the plot (default is True)
        :type iplot: bool
        :param field: 'B', 'V' or 'T' (magnetic field, velocity field or temperature)
        :type field: str
        :param r: an integer to characterise which file we want to plot
        :type r: int
        :param precision: single or double precision
        :type precision: str
        """

        logFiles = scanDir('log.*')
        if len(logFiles) != 0:
            MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])
        else:
            str1 = 'Aspect ratio ?\n'
            self.radratio = float(input(str1))

        self.rcmb = 1./(1.-self.radratio)
        ricb = self.radratio/(1.-self.radratio)

        files = scanDir('%s_coeff_r%i.%s' % (field,r,tag))

        # Read the B_coeff files (by stacking the different tags)
        data = []
        for k, file in enumerate(files):
            print('Reading %s' % file)
            f = npfile(file, endian='B')
            out = f.fort_read('3i4,%s' % precision)[0]
            self.l_max_cmb, self.minc, n_data = out[0]
            self.m_max_cmb = int((self.l_max_cmb/self.minc)*self.minc)
            self.radius = out[1]

            while 1:
                try:
                    data.append(f.fort_read(precision))
                except TypeError:
                    break
        data = N.array(data, dtype=precision)

        self.ell = N.arange(self.l_max_cmb+1)
        self.nstep = data.shape[0]

        # Get time
        self.time = N.zeros(self.nstep, dtype=precision)
        self.time = data[:, 0]

        # Rearange and get Gauss coefficients
        self.wlm = N.zeros((self.nstep, self.l_max_cmb+1, self.m_max_cmb+1), 'Complex64')
        self.dwlm = N.zeros((self.nstep, self.l_max_cmb+1, self.m_max_cmb+1), 'Complex64')
        self.zlm = N.zeros((self.nstep, self.l_max_cmb+1, self.m_max_cmb+1), 'Complex64')

        # wlm
        # Axisymmetric coefficients (m=0)
        self.wlm[:, 1:, 0] = data[:, 1:self.l_max_cmb+1]
        # Other coefficients (m!=0)
        k = self.l_max_cmb+1
        for m in range(self.minc, self.l_max_cmb+1, self.minc):
            for l in range(m, self.l_max_cmb+1):
                self.wlm[:, l, m] = data[:, k]+1j*data[:, k+1]
                k += 2

        # dwlm
        self.dwlm[:, 1:, 0] = data[:, k:k+self.l_max_cmb]
        k += self.l_max_cmb
        for m in range(self.minc, self.l_max_cmb+1, self.minc):
            for l in range(m, self.l_max_cmb+1):
                self.dwlm[:, l, m] = data[:, k]+1j*data[:, k+1]
                k += 2
        # zlm
        self.zlm[:, 1:, 0] = data[:, k:k+self.l_max_cmb]
        k += self.l_max_cmb
        for m in range(self.minc, self.l_max_cmb+1, self.minc):
            for l in range(m, self.l_max_cmb+1):
                self.zlm[:, l, m] = data[:, k]+1j*data[:, k+1]
                k += 2

        # ddw in case B is stored
        if field == 'B':
            self.ddwlm = N.zeros((self.nstep, self.l_max_cmb+1, self.m_max_cmb+1), 'Complex64')
            self.ddwlm[:, 1:, 0] = data[:, k:k+self.l_max_cmb]
            k += self.l_max_cmb
            for m in range(self.minc, self.l_max_cmb+1, self.minc):
                for l in range(m, self.l_max_cmb+1):
                    self.ddwlm[:, l, m] = data[:, k]+1j*data[:, k+1]
                    k += 2

        #self.epolLM = 0.5*self.ell*(self.ell+1)* (self.ell*(self.ell+1)* \
        #                            abs(self.wlm)**2+abs(self.dwlm)**2)
        self.epolAxiL = 0.5*self.ell*(self.ell+1)*(self.ell*(self.ell+1)* \
                                    abs(self.wlm[:,:,0])**2+abs(self.dwlm[:,:,0])**2)
        #self.etorLM = 0.5*self.ell*(self.ell+1)*abs(self.zlm)**2
        self.etorAxiL = 0.5*self.ell*(self.ell+1)*abs(self.zlm[:,:,0])**2

        #epolTot = self.epolLM.sum(axis=1).sum(axis=1)
        #etorTot = self.etorLM.sum(axis=1).sum(axis=1)
        etorAxiTot = self.etorAxiL.sum(axis=1)
        epolAxiTot = self.epolAxiL.sum(axis=1)
        
        if iplot:
            P.plot(self.time, epolTot)
            P.plot(self.time, etorTot)
            P.plot(self.time, epolAxiTot)
            P.plot(self.time, etorAxiTot)
        
        """
        self.dglmdt = deriv(self.time, self.glm.T, axis=2)
        self.dhlmdt = deriv(self.time, self.hlm.T, axis=2)
        self.dglmdt = self.dglmdt.T
        self.dhlmdt = self.dhlmdt.T

        #print(self.glm[-1, 1, 0], self.glm[-1, 2, 0], self.glm[-1, 2, 2])

        # Time-averaged Gauss coefficient
        facT = 1./(self.time[-1]-self.time[0])
        self.glmM = facT * N.trapz(self.glm, self.time, axis=0)
        self.hlmM = facT * N.trapz(self.hlm, self.time, axis=0)

        # Magnetic energy (Lowes)
        self.El = (self.ell+1)*(self.glm**2+self.hlm**2).sum(axis=2)
        self.Em = N.zeros((self.nstep, self.m_max_cmb+1), dtype=precision)
        # For m, we need to unfold the loop in case of minc != 1
        for m in range(0, self.m_max_cmb+1, self.minc):
            self.Em[:,m] = ((self.ell+1)*(self.glm[:, :, m]**2+self.hlm[:, :, m]**2)).sum(axis=1)
        #self.Em = ((self.ell+1)*(self.glm**2+self.hlm**2)).sum(axis=1)

        # Time-averaged energy
        self.ElM = facT * N.trapz(self.El, self.time, axis=0)
        self.EmM = facT * N.trapz(self.Em, self.time, axis=0)

        # Secular variation
        self.ESVl = (self.ell+1)*(self.dglmdt**2+self.dhlmdt**2).sum(axis=2)
        self.ESVlM = facT * N.trapz(self.ESVl, self.time, axis=0)
        self.taul = N.sqrt(self.ElM[1:]/self.ESVlM[1:])

        #if iplot:
        #    self.plot()
        """
