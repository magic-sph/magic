# -*- coding: utf-8 -*-
from magic import MagicGraph, BLayers
from .spectralTransforms import SpectralTransforms
from .setup import labTex
from .libmagic import symmetrize
import matplotlib.pyplot as plt
import numpy as np
import sys



class ExtraPot:
    """
    This class is used to compute the potential field extrapolation of the magnetic
    field in an arbitrary outer spherical shell domain. It takes as an input
    the magnetic field at the CMB.
    """

    def __init__(self, rcmb, brcmb, minc, ratio_out=2., nrout=32, cutCMB=False,
                 deminc=True):
        """
        :param bcmb: the surface radial field, array of dimension [np, nt]
        :type bcmb: numpy.ndarary
        :param rcmb: the value of the radius at the surface
        :type rcmb: float
        :param minc: azimuthal symmetry
        :type minc: int
        :param ratio_out: the ratio of the outer sphere radius to the surface radius
        :type ratio_out: float
        :param nrout: the number of radial point (linearly spaced) of the
                      extrapolated field in the outer spherical domain
        :type nrout: int
        :param cutCMB: a logical if one wants to remove the first grid point
                       (useful if one then wants to merge the graphic file with
                       the extrapolation)
        :type cutCMB: bool
        :param deminc: a logical to indicate if one wants do get rid of the
                       possible azimuthal symmetry
        :type deminc: bool
        """
        self.rcmb = rcmb
        self.brcmb = brcmb
        self.minc = minc
        self.nrout = nrout
        self.np, self.nt = self.brcmb.shape


        self.rout = np.linspace(self.rcmb, ratio_out*rcmb, self.nrout)
        if cutCMB:
            self.rout = self.rout[1:]
            self.nrout = self.nrout -1

        l_max = (20*self.nt)/30
        m_max = (l_max/self.minc)*self.minc
        lm_max = m_max*(l_max+1)/minc-m_max*(m_max-minc)/(2*minc)+(l_max+1-m_max)

        sh = SpectralTransforms(l_max, self.minc, lm_max, self.nt)
        brlm = sh.spat_spec(self.brcmb)
        bpolcmb = np.zeros_like(brlm)
        bpolcmb[1:] = brlm[1:]/(sh.ell[1:]*(sh.ell[1:]+1))*self.rcmb**2
        btor = np.zeros_like(brlm)

        self.brout = np.zeros((self.np, self.nt, self.nrout), dtype=self.brcmb.dtype)
        self.btout = np.zeros_like(self.brout)
        self.bpout = np.zeros_like(self.brout)


        for k, radius  in enumerate(self.rout):
            radratio = self.rcmb/radius
            # Potential extrapolation
            bpol = bpolcmb * radratio**sh.ell
            # Derivative is analytical:
            dbpoldr = -sh.ell/radius * bpol
            brlm = sh.ell*(sh.ell+1) * bpol / radius**2
            self.brout[..., k] = sh.spec_spat(brlm)
            self.btout[..., k], self.bpout[..., k] = sh.spec_spat(dbpoldr, btor)
            self.btout[..., k] /= radius
            self.bpout[..., k] /= radius

        if deminc:
            self.brout = symmetrize(self.brout, self.minc)
            self.btout = symmetrize(self.btout, self.minc)
            self.bpout = symmetrize(self.bpout, self.minc)

    def avg(self, field='br', levels=12, cm='RdYlBu_r', normed=True,
            vmax=None, vmin=None):
        """
        A small routine to plot the azimuthal averages of the extrapolated
        fields.

        :param field: the quantity you want to plot: 'br' or 'bp'
        :type field: str
        :param levels: the number of contour levels
        :type levels: int
        :param cm: the name of the colormap
        :type cm: str
        :param vmax: maximum value of the contour levels
        :type vmax: float
        :param vmin: minimum value of the contour levels
        :type vmin: float
        :param normed: when set to True, the colormap is centered around zero.
                       Default is True, except for entropy/temperature plots.
        :type normed: bool
        """

        if field == 'br':
            if labTex:
                label = r'$B_r$'
            else:
                label = 'Br'
            data = self.brout
        elif field == 'bp':
            if labTex:
                label = r'$B_\phi$'
            else:
                label = 'Bphi'
            data = self.bpout

        phiavg = data.mean(axis=0)
        th = np.linspace(0, np.pi, self.nt)
        rr, tth = np.meshgrid(self.rout, th)
        xx = rr * np.sin(tth)
        yy = rr * np.cos(tth)

        fig = plt.figure(figsize=(4,8))
        fig.subplots_adjust(top=0.99, bottom=0.01, right=0.99, left=0.01)
        ax = fig.add_subplot(111, frameon=False)
        cmap = plt.get_cmap(cm)
        if vmax is not None and vmin is not None:
            normed = False
            cs = np.linspace(vmin, vmax, levels)
            im = ax.contourf(xx, yy, phiavg, cs, cmap=cmap, extend='both')
        else:
            cs = levels
            im = ax.contourf(xx, yy, phiavg, cs, cmap=cmap)
        ax.plot(self.rcmb*np.sin(th), self.rcmb*np.cos(th), 'k-')
        ax.axis('off')
