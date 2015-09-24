# -*- coding: utf-8 -*-
from magic import npfile, scanDir, MagicSetup
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
    if m > 0:
        fac *= N.sqrt(2.)
    glm = scale_b*ratio_cmb_surface**(ell+2.)/rcmb**2*fac*alm
    hlm = -scale_b*ratio_cmb_surface**(ell+2.)/rcmb**2*fac*blm

    return glm, hlm


class BcoeffCmb(MagicSetup):
    """
    This class allows to read the :ref:`B_coeff_cmb.TAG <secCoeffFiles>` files.
    It first read the poloidal potential at the CMB and then transform
    it to the Gauss coefficients :math:`g_{\ell m}` and :math:`h_{\ell m}`
    using the getGauss function.

    >>> # Reads the files B_coeff_cmb.testa, B_coeff_cmb.testb
    >>> # and B_coeff_cmb.testc and stack them in one single time series
    >>> cmb = BcoeffCmb(tag='test[a-c]')
    >>> print(cmb.ell, cmb.glm) # print \ell and g_{\ell m}
    >>> print(cmb.glm[:, 1, 0]) # time-series of the axisymmetric dipole
    >>> plot(cmb.time, cmb.dglmdt[:, 1, 0]) # Secular variation of the dipole
    """
    
    def __init__(self, tag, ratio_cmb_surface=1, scale_b=1, iplot=True,
                 precision='Float64'):
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
        """

        logFiles = scanDir('log.*')
        if len(logFiles) != 0:
            MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])
        else:
            str1 = 'Aspect ratio ?\n'
            self.radratio = float(input(str1))

        rcmb = 1./(1.-self.radratio)
        ricb = self.radratio/(1.-self.radratio)

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
        data = N.array(data, dtype=precision)

        self.ell = N.arange(self.l_max_cmb+1)
        self.nstep = data.shape[0]

        # Get time
        self.time = N.zeros(self.nstep, precision)
        self.time = data[:, 0]

        # Rearange and get Gauss coefficients
        self.alm = N.zeros((self.nstep, self.l_max_cmb+1, self.m_max_cmb+1), precision)
        self.blm = N.zeros((self.nstep, self.l_max_cmb+1, self.m_max_cmb+1), precision)
        self.glm = N.zeros((self.nstep, self.l_max_cmb+1, self.m_max_cmb+1), precision)
        self.hlm = N.zeros((self.nstep, self.l_max_cmb+1, self.m_max_cmb+1), precision)

        # Axisymmetric coefficients (m=0)
        self.alm[:, 1:, 0] = data[:, 1:self.l_max_cmb+1]
        self.glm[:, 1:, 0], self.hlm[:, 1:, 0] = getGauss(self.alm[:, 1:, 0], 
                     self.blm[:, 1:, 0], self.ell[1:], 0, scale_b, 
                     ratio_cmb_surface, rcmb)
        # Other coefficients (m!=0)
        k = self.l_max_cmb+1
        for m in range(self.minc, self.l_max_cmb+1, self.minc):
            for l in range(m, self.l_max_cmb+1):
                self.alm[:, l, m] = data[:, k]
                self.blm[:, l, m] = data[:, k+1]
                self.glm[:, l, m], self.hlm[:, l, m] = getGauss(self.alm[:, l, m], 
                                    self.blm[:, l, m], l, m,
                                    scale_b, ratio_cmb_surface, rcmb)
                k += 2

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
        self.Em = N.zeros((self.nstep, self.m_max_cmb+1), precision)
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

        if iplot:
            self.plot()


    def plot(self):
        """
        Display some results when iplot is set to True
        """
        fig = P.figure()
        ax = fig.add_subplot(211)
        ax.semilogy(self.ell[1:], self.ElM[1:], 'b-o')
        if labTex:
            ax.set_xlabel(r'$\ell$')
        else:
            ax.set_xlabel('Degree l')
        ax.set_ylabel('Magnetic energy')
        ax.set_xlim(1., self.l_max_cmb)

        ax1 = fig.add_subplot(212)
        ax1.semilogy(self.ell[0:self.m_max_cmb+1:self.minc], self.EmM[::self.minc], 'b-o')   
        if labTex:
            ax1.set_xlabel(r'$m$')
        else:
            ax1.set_xlabel('Order m')
        ax1.set_ylabel('Magnetic energy')

        fig1 = P.figure()
        ax = fig1.add_subplot(111)
        ax.loglog(self.ell[1:], self.taul, 'b-o')
        if labTex:
            ax.set_xlabel(r'$\ell$')
            ax.set_ylabel(r'$\tau_\ell$')
        else:
            ax.set_xlabel('Degree l')
            ax.set_ylabel('tau  l')
        ax.set_xlim(1, self.l_max_cmb)

        fig2 = P.figure()
        ax = fig2.add_subplot(111)
        ax.plot(self.time, self.glm[:, 1, 0], label='g10')
        ax.plot(self.time, self.glm[:, 2, 0], label='g20')
        ax.plot(self.time, self.glm[:, 3, 0], label='g30')
        ax.set_xlabel('Time')
        ax.set_ylabel('Gauss coefficients')


class Bcoeff_r(MagicSetup):
    """
    This class allows to read the :ref:`B_coeff_r#.TAG <secBcoeffrFile>`
    and :ref:`V_coeff_r#.TAG <secVcoeffrFile>` files.
    It reads the poloidal and toroidal potentials and reconstruct the time
    series (or the energy) contained in any given mode.

    >>> # Reads the files V_coeff_r2.test*
    >>> cr = Bcoeff_r(tag='test*', field='V', r=2)
    >>> print(cr.ell, cr.wlm) # print \ell and w_{\ell m}
    >>> # Time-evolution of the poloidal energy in the (\ell=10, m=10) mode
    >>> plot(cr.time, cr.epolLM[:, 10, 10]) 
    """
    
    def __init__(self, tag, ratio_cmb_surface=1, scale_b=1, iplot=True,
                 field='B', r=1, precision='Float64'):
        """
        :param tag: if you specify a pattern, it tires to read the corresponding files
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

        rcmb = 1./(1.-self.radratio)
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

        print(k, data.shape)
        # ddw in case B is stored
        if field == 'B':
            self.ddwlm = N.zeros((self.nstep, self.l_max_cmb+1, self.m_max_cmb+1), 'Complex64')
            self.ddwlm[:, 1:, 0] = data[:, k:k+self.l_max_cmb]
            k += self.l_max_cmb
            for m in range(self.minc, self.l_max_cmb+1, self.minc):
                for l in range(m, self.l_max_cmb+1):
                    self.ddwlm[:, l, m] = data[:, k]+1j*data[:, k+1]
                    k += 2

        self.epolLM = 0.5*self.ell*(self.ell+1)* (self.ell*(self.ell+1)* \
                                    abs(self.wlm)**2+abs(self.dwlm)**2)
        self.epolAxiL = 0.5*self.ell*(self.ell+1)*(self.ell*(self.ell+1)* \
                                    abs(self.wlm[:,:,0])**2+abs(self.dwlm[:,:,0])**2)
        self.etorLM = 0.5*self.ell*(self.ell+1)*abs(self.zlm)**2
        self.etorAxiL = 0.5*self.ell*(self.ell+1)*abs(self.zlm[:,:,0])**2

        epolTot = self.epolLM.sum(axis=1).sum(axis=1)
        etorTot = self.etorLM.sum(axis=1).sum(axis=1)
        etorAxiTot = self.etorAxiL.sum(axis=1)
        epolAxiTot = self.epolAxiL.sum(axis=1)

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
