# -*- coding: utf-8 -*-
import numpy as np
from .setup import buildSo

if buildSo:
    import magic.legendre as leg


class SpectralTransforms(object):
    """
    This python class is used to compute Legendre and Fourier transforms
    from spectral to physical space. It works in two steps: one first needs
    to initialize the transform

    >>> sh = SpectralTransforms( l_max=256, lm_max=33153, n_theta_max=384)
    >>> print(Tlm[:, 10].shape) # lm_max (Temperature at ir=10)
    >>> T = sh.spec_spat(Tlm) # T[n_phi_max, n_theta_max]
    """

    def __init__(self, l_max=32, minc=1, lm_max=561, n_theta_max=64,
                 verbose=True):
        """
        :param l_max: maximum spherical harmonic degree
        :type l_max: int
        :param minc: azimuthal symmetry
        :type minc: int
        :param lm_max: maximum l,m combination
        :type lm_max: int
        :param n_theta_max: number of grid points in the latitudinal direction
        :type n_theta_max: int
        :param verbose: some info about the SHT layout
        :type verbose: bool
        """
        self._legF90 = leg.legendre
        self._legF90.init(l_max, minc, lm_max, n_theta_max)
        self.l_max = self._legF90.l_max
        self.minc = self._legF90.minc
        self.lm_max = self._legF90.lm_max
        self.n_theta_max = self._legF90.n_theta_max
        self.n_phi_max = self._legF90.n_phi_max
        self.m_max = int((self.l_max/self.minc) * self.minc)

        if verbose:
            print('Spectral transform setup:')
            print('l_max, m_max, minc, lm_max: {}, {}, {}, {}'.format(
                self.l_max, self.m_max, self.minc, self.lm_max))
            print('n_phi_max, n_theta_max: {}, {}'.format(
                self.n_phi_max, self.n_theta_max))

        self.colat = self._legF90.sinth

        self.idx = np.zeros((self.l_max+1, self.m_max+1), 'i')
        self.ell = np.zeros((self.lm_max), 'i')

        self.idx[0:self.l_max+2, 0] = np.arange(self.l_max+1)
        self.ell[0:self.l_max+2] = np.arange(self.l_max+2)
        k = self.l_max+1
        for m in range(self.minc, self.l_max+1, self.minc):
            for l in range(m, self.l_max+1):
                self.idx[l, m] = k
                self.ell[self.idx[l,m]] = l
                k +=1

    def spec_spat(self, *args, **kwargs):
        """
        This subroutine computes a transfrom from spectral to spatial
        for all latitudes. It returns either  one or two 2-D arrays
        (dimension(n_phi_max,n_theta_max)) depending if only the poloidal
        or both the poloidal and the toroidal potentials are given as
        input quantities.

        >>> print(wlmr.shape) # lm_max
        >>> vr = spec_spat_equat(wlmr)
        >>> print(vr.shape) # n_phi, n_theta
        >>> vt, vp = spec_spat_equat(dwdrlmr, zlmr)
        """
        if 'l_axi' in kwargs:
            l_axi = kwargs['l_axi']
            if l_axi:
                n_phi = 1
            else:
                n_phi = self.n_phi_max
        else:
            n_phi = self.n_phi_max

        if len(args) == 1:
            polo = args[0]
            out = self._legF90.specspat_scal(polo, self.n_theta_max, n_phi)
            if n_phi > 1:
                out = np.fft.ifft(out, axis=0)*self.n_phi_max
            out = out.real
            return out
        elif len(args) == 2:
            polo = args[0]
            toro = args[1]
            vt, vp = self._legF90.specspat_vec(polo, toro, self.n_theta_max, n_phi)
            if n_phi > 1:
                vt = np.fft.ifft(vt, axis=0)*self.n_phi_max
                vp = np.fft.ifft(vp, axis=0)*self.n_phi_max
            vt = vt.real
            vp = vp.real

            return vt, vp

    def spec_spat_dtheta(self, polo, l_axi=False):
        """
        This routine computes the theta-derivative and the transform from spectral
        to spatrial spaces. It returns a 2-D array of dimension (n_phi,n_theta)

        >>> p = MagicPotential('V')
        >>> vrlm = p.pol*p.ell*(p.ell+1)/p.radius[ir]**2/p.rho0[ir] # vr at r=ir
        >>> dvrdt = p.sh.spec_spat_dtheta(vrlm) # theta-derivative of vr

        :param polo: the input array(lm_max) in spectral space
        :type polo: numpy.ndarray
        :param l_axi: switch to True, if only the axisymmetric field is needed
        :type l_axi: bool
        :returns: the theta derivative in the physical space (n_phi, n_theta)
        :rtype: numpy.ndarray
        """
        if l_axi:
            n_phi = 1
        else:
            n_phi = self.n_phi_max

        out = self._legF90.specspat_dtheta(polo, self.n_theta_max, n_phi)
        if n_phi > 1:
            out = np.fft.ifft(out, axis=0)*self.n_phi_max
        out = out.real

        return out

    def spec_spat_dphi(self, polo):
        """
        This routine computes the phi-derivative and the transform from spectral
        to spatrial spaces. It returns a 2-D array of dimension (n_phi,n_theta)

        >>> p = MagicPotential('V')
        >>> vrlm = p.pol*p.ell*(p.ell+1)/p.radius[ir]**2/p.rho0[ir] # vr at r=ir
        >>> dvrdp = p.sh.spec_spat_dphi(vrlm) # phi-derivative of vr

        :param polo: the input array(lm_max) in spectral space
        :type polo: numpy.ndarray
        :returns: the phi derivative in the physical space (n_phi, n_theta)
        :rtype: numpy.ndarray
        """
        out = self._legF90.specspat_dphi(polo, self.n_theta_max, self.n_phi_max)
        out = np.fft.ifft(out, axis=0)*self.n_phi_max
        # since dPhilm contains a 1/sin(theta) one has to multiply by sin(theta)
        out = out.real * np.sin(self.colat)

        return out

    def spec_spat_equat(self, *args):
        """
        This subroutine computes a transfrom from spectral to spatial
        at the equator. It returns either one or two 1-D arrays
        (dimension(n_phi_max)) depending if only the poloidal or both the
        poloidal and the toroidal potentials are given as input quantities.

        >>> print(wlmr.shape) # lm_max
        >>> vr = spec_spat_equat(wlmr)
        >>> print(vr.shape) # n_phi
        >>> vt, vp = spec_spat_equat(dwdrlmr, zlmr)
        """
        if len(args) == 1:
            polo = args[0]
            out = np.zeros((self.n_phi_max), np.complex128)
            self._legF90.specspat_equat_scal(polo, out)
            out = np.fft.ifft(out)*self.n_phi_max
            out = out.real

            return out

        elif len(args) == 2:
            polo = args[0]
            toro = args[1]
            vt = np.zeros((self.n_phi_max), np.complex128)
            vp = np.zeros((self.n_phi_max), np.complex128)
            self._legF90.specspat_equat_vec(polo, toro, vt, vp)
            vt = np.fft.ifft(vt)*self.n_phi_max
            vp = np.fft.ifft(vp)*self.n_phi_max
            vt = vt.real
            vp = vp.real

            return vt, vp

    def spat_spec(self, *args):
        """
        This subroutine computes a transfrom from spatial representation
        (n_phi,n_theta) to spectral representation (lm_max). It returns
        one complex 1-D array (dimension(n_phi_max))

        >>> gr = MagicGraph()
        >>> sh = SpectralTransforms(gr.l_max, gr.minc, gr.lm_max, gr.n_theta_max)
        >>> vr = gr.vr[:,:,30] # Radius ir=30
        >>> vrlm = sh.spat_spec(vr) # vrlm is a complex array (lm_max)
        >>> # Caculation of the poloidal potential from vr:
        >>> wlm = np.zeros_like(vrlm)
        >>> wlm[1:] = vrlm[1:]/(sh.ell[1:]*(sh.ell[1:]+1))*gr.radius[30]**2
        >>> # Spheroidal/Toroidal transform
        >>> vtlm, vplm = spec_spat(gr.vtheta, gr.vphi)

        :param input: input array in the physical space (n_phi,n_theta)
        :type input: numpy.ndarray
        :returns: output array in the spectral space (lm_max)
        :rtype: numpy.ndarray
        """
        if len(args) == 1:
            input = args[0]
            out = np.fft.fft(input, axis=0)/self.n_phi_max
            outLM = self._legF90.spatspec(out, self.lm_max)

            return outLM
        elif len(args) == 2:
            in1 = args[0]
            in2 = args[1]
            out1 = np.fft.fft(in1, axis=0)/self.n_phi_max
            out2 = np.fft.fft(in2, axis=0)/self.n_phi_max

            utLM, upLM = self._legF90.spatspec_sphertor(out1, out2, self.lm_max)
            return utLM, upLM


if __name__ == '__main__':
    sh = SpectralTransforms( l_max=256, lm_max=33153, n_theta_max=384)
    print( a.lm_max )
