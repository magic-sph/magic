# -*- coding: utf-8 -*-
from magic import npfile, MagicSetup, scanDir, chebgrid
import numpy as np
import matplotlib.pyplot as plt


class MagicRSpec(MagicSetup):
    """
    This class allows to read the :ref:`rB[r|p]Spec.TAG files <secrBspecFiles>`.
    Those files contain the time-evolution of the poloidal/toroidal magnetic energy
    for all radii and for spherical harmonic degrees from 1 to 6. This is an unformatted
    fortran file.

    >>> # Read all the `BrSpec.test*` files in the current working directory and
    >>> # stack them.
    >>> rsp = MagicRSpec(tag='test*', field='Br')
    """

    def __init__(self, tag, field='Br', precision=np.float32, avg=False):
        """
        :param tag: if you specify a pattern, it tries to read the corresponding
                    files and stack them.
        :type tag: str
        :param field: nature of the radial spectra. Possible choices are
                      'Bt' or 'Bp'
        :type field: str
        :param precision: single or double precision (default single, i.e.
                          np.float32)
        :type precision: str
        :param avg: when set to True, display time averaged quantities
        :type avg: bool
        """

        logFiles = scanDir('log.*')
        if len(logFiles) != 0:
            MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])
            self.n_r_max = int(self.n_r_max)
            self.n_r_ic_max = int(self.n_r_ic_max)
        else:
            n_r_max = 'n_r_max ?\n'
            self.n_r_max = int(input(str1))
            n_r_ic_max = 'n_r_ic_max ?\n'
            self.n_r_ic_max = int(input(str1))
            str1 = 'Aspect ratio ?\n'
            self.radratio = float(input(str1))

        self.n_r_tot = self.n_r_max+self.n_r_ic_max

        self.ricb = self.radratio/(1.-self.radratio)
        self.rcmb = 1./(1.-self.radratio)
        outerCoreGrid = chebgrid(self.n_r_max-1, self.rcmb, self.ricb)
        n_r_ic_tot = 2*self.n_r_ic_max-1
        innerCoreGrid = chebgrid(n_r_ic_tot-1, self.ricb, -self.ricb)

        self.radius = np.zeros((self.n_r_tot-1), dtype=precision)

        self.radius[:self.n_r_max] = outerCoreGrid
        self.radius[self.n_r_max-1:] = innerCoreGrid[:self.n_r_ic_max]

        pattern = 'r{}'.format(field) +'Spec'
        files = scanDir('{}.{}'.format(pattern, tag))

        # Read the rB[rp]Spec.TAG files (stack them)
        data = []
        for k, file in enumerate(files):
            print('Reading {}'.format(file))
            f = npfile(file, endian='B')

            while 1:
                try:
                    data.append(f.fort_read(precision))
                except TypeError:
                    break
        data = np.array(data, dtype=precision)

        # Time (every two lines only)
        self.time = data[::2, 0]

        # Poloidal/Toroidal energy for all radii for the 6 first spherical harmonic
        # degrees
        self.e_pol = np.zeros((len(self.time), self.n_r_tot-1, 6), dtype=precision)
        self.e_pol_axi = np.zeros_like(self.e_pol)

        self.e_pol[:, :, 0] = data[::2, 1:self.n_r_tot]
        self.e_pol[:, :, 1] = data[::2, self.n_r_tot-1:2*(self.n_r_tot-1)]
        self.e_pol[:, :, 2] = data[::2, 2*(self.n_r_tot-1):3*(self.n_r_tot-1)]
        self.e_pol[:, :, 3] = data[::2, 3*(self.n_r_tot-1):4*(self.n_r_tot-1)]
        self.e_pol[:, :, 4] = data[::2, 4*(self.n_r_tot-1):5*(self.n_r_tot-1)]
        self.e_pol[:, :, 5] = data[::2, 5*(self.n_r_tot-1):6*(self.n_r_tot-1)]

        self.e_pol_axi[:, :, 0] = data[1::2, 1:self.n_r_tot]
        self.e_pol_axi[:, :, 1] = data[1::2, self.n_r_tot-1:2*(self.n_r_tot-1)]
        self.e_pol_axi[:, :, 2] = data[1::2, 2*(self.n_r_tot-1):3*(self.n_r_tot-1)]
        self.e_pol_axi[:, :, 3] = data[1::2, 3*(self.n_r_tot-1):4*(self.n_r_tot-1)]
        self.e_pol_axi[:, :, 4] = data[1::2, 4*(self.n_r_tot-1):5*(self.n_r_tot-1)]
        self.e_pol_axi[:, :, 5] = data[1::2, 5*(self.n_r_tot-1):6*(self.n_r_tot-1)]

        if avg:
            self.plotAvg()

    def plotAvg(self):
        """
        Plotting function for time-averaged profiles
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.radius, self.e_pol[:, :, 0].mean(axis=0), 'b-', label='Dip.')
        ax.plot(self.radius, abs(self.e_pol_axi[:, :, 0].mean(axis=0)), 'b--',
                label='Dip. axi.')

        ax.plot(self.radius, self.e_pol[:, :, 1].mean(axis=0), 'r-', label='Quad.')
        ax.plot(self.radius, abs(self.e_pol_axi[:, :, 1].mean(axis=0)), 'r--',
                label='Quad. axi.')

        ax.plot(self.radius, self.e_pol[:, :, 2].mean(axis=0), 'g-', label='Octu.')
        ax.plot(self.radius, abs(self.e_pol_axi[:, :, 2].mean(axis=0)), 'g--',
                label='Octu. axi.')

        ax.plot(self.radius, self.e_pol[:, :, 3].mean(axis=0), 'g-', label='Hexadeca.')
        ax.plot(self.radius, abs(self.e_pol_axi[:, :, 3].mean(axis=0)), 'g--',
                label='Hexadeca. axi.')

        ax.set_xlabel('Radius')
        ax.set_ylabel('Energy')

        ax.legend(frameon=False, loc='best')



if __name__ == '__main__':
    r = MagicRSpec(tag='test', field='Bp', avg=True)
    plt.show()
