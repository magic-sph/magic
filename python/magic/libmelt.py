# -*- coding: utf-8 -*-
import numpy as np
import os
import re
import matplotlib.pyplot as plt
from .plotlib import radialContour
from .libmagic import scanDir, symmetrize
from .npfile import npfile
from .movie import getNlines
from .log import MagicSetup
from scipy.interpolate import interp1d


class MagicMelt(MagicSetup):
    """
    This python class is used to read the rmelt.TAG files produced when the
    computation of the phase field is switched on.

    >>> rm = MagicMelt()
    >>> # This reads the most recent rmelt.TAG file

    >>> rm = MagicMelt(all=True)
    >>> # This reads all the rmelt.TAG files found in the working directory
    >>> print(rm.time, rm.rmelt)
    """

    def __init__(self, tag=None, datadir='.', endian='B', precision=np.float32,
                 iplot=True, levels=17, cm='viridis', all=False, nstep=1):
        """
        :param tag: extension TAG of the rmelt file. If not specified, the
                    most recent rmelt.TAG file found in the directory will
                    be selected.
        :type tag: str
        :param datadir: directory of the rmelt file (default is . )
        :type datadir: str
        :param precision: single or double precision (default np.float32)
        :type precision: str
        :param endian: endianness of the file ('B' or 'l')
        :type endian: str
        :param iplot: to plot the output, default is True
        :type iplot: bool
        :param levels: the number of contour levels
        :type levels: int
        :param cm: the name of the color map
        :type cm: str
        :param all: when set to True, the complete time series is reconstructed by
                    stacking all the corresponding files from the working directory
                    (default False)
        :type all: bool
        :param nstep: the stepping between two timesteps
        :type nstep: int
        """

        prefix = 'rmelt'
        pattern = os.path.join(datadir, 'log.*')
        logFiles = scanDir(pattern)
        self.endian = endian
        self.precision = precision

        if tag is not None:
            pattern = os.path.join(datadir, '{}.{}'.format(prefix, tag))
            files = scanDir(pattern)

            # Either the log.tag directly exists and the setup is easy to obtain
            if os.path.exists(os.path.join(datadir, 'log.{}'.format(tag))):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.{}'.format(tag))
            # Or the tag is a bit more complicated and we need to find
            # the corresponding log file
            else:
                st = os.path.join(datadir, '{}\.(.*)'.format(prefix))
                mask = re.compile(st)
                if mask.match(files[-1]):
                    ending = mask.search(files[-1]).groups(0)[0]
                    pattern = os.path.join(datadir, 'log.{}'.format(ending))
                    if logFiles.__contains__(pattern):
                        MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml='log.{}'.format(ending))
            for k, file in enumerate(files):
                if k == 0:
                    self.colatitude, self.time, self.rmelt = self.read_one(file)
                else:
                    colat, time, rmelt = self.read_one(file)
                    self.stack_two(colat, time, rmelt)


        # If no tag is specified, the most recent is plotted
        elif not all:
            if len(logFiles) != 0:
                MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])
                name = '{}.{}'.format(prefix, self.tag)
                filename = os.path.join(datadir, name)
            else:
                mot = '{}.*'.format(prefix)
                dat = [(os.stat(i).st_mtime, i) for i in glob.glob(mot)]
                dat.sort()
                filename = dat[-1][1]

            self.colatitude, self.time, self.rmelt = self.read_one(filename)

        # If no tag is specified but all=True, all the directory is plotted
        else:
            if len(logFiles) != 0:
                MagicSetup.__init__(self, quiet=True, nml=logFiles[-1])
            pattern = os.path.join(datadir, '{}.*'.format(prefix))
            files = scanDir(pattern)
            for k, file in enumerate(files):
                if k == 0:
                    self.colatitude, self.time, self.rmelt = self.read_one(file)
                else:
                    colat, time, rmelt = self.read_one(file)
                    self.stack_two(colat, time, rmelt)

        if nstep != 1:
            self.time = self.time[::nstep]
            self.rmelt = self.rmelt[::nstep, :]

        if iplot:
            self.plot(levels, cm)

    def read_one(self, filename):
        """
        This routine is used to read one file

        :param filename: name of the rmelt file
        :type filename: str
        :returns: a list of 3 numpy.ndarrays, the first one is the colatitude,
                  the second one the time and the third one the melt radius
        :rtype: list
        """

        nlines = getNlines(filename, endian=self.endian,
                           precision=self.precision)
        nlines -= 2  # Remove the two lines of the header

        infile = npfile(filename, endian=self.endian)
        n_theta_max = infile.fort_read('i4')[0]
        colatitude = infile.fort_read(self.precision)
        time = np.zeros((nlines), self.precision)
        rmelt = np.zeros((nlines, n_theta_max), self.precision)

        for i in range(nlines):
            data = infile.fort_read(self.precision)
            time[i] = data[0]
            rmelt[i, :] = data[1:]

        infile.close()

        return colatitude, time, rmelt

    def stack_two(self, theta_new, time_new, rmelt_new):
        """
        This routine is used to stack two time series

        :param theta_new: the colatitude of the new file
        :type theta_new: numpy.ndarray
        :param time_new: the time of the new file
        :type time_new: numpy.ndarray
        :param rmelt: the melt radius of the new file
        :type rmelt: numpy.ndarray
        """

        # If the number of latitudinal grid points has changed, one
        # needs to interpolate
        if len(self.colatitude) != len(theta_new):
            it = interp1d(self.colatitude, self.rmelt, axis=-1, 
                          fill_value='extrapolate')
            self.rmelt = it(theta_new)
            self.colatitude = theta_new

        if abs(time_new[0]-self.time[-1]) <= 1e-10:
            self.time = np.concatenate((self.time, time_new[1:]))
            self.rmelt = np.concatenate((self.rmelt, rmelt_new[1:, :]), axis=0)
        else:
            self.time = np.concatenate((self.time, time_new))
            self.rmelt = np.concatenate((self.rmelt, rmelt_new), axis=0)

    def plot(self, levels=17, cm='viridis'):
        """
        Plotting routine for rmelt.TAG files

        :param levels: the number of contour levels
        :type levels: int
        :param cm: the name of the color map
        :type cm: str
        """

        fig = plt.figure()
        ax = fig.add_subplot(111)

        im = ax.contourf(self.time, self.colatitude, self.rmelt.T, levels, 
                         cmap=plt.get_cmap(cm))
        ax.set_xlabel('Time')
        ax.set_ylabel('Colatitude')
        fig.colorbar(im)
        fig.tight_layout()

        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        rmelt_mean = self.rmelt.mean(axis=0)
        rmelt_std = self.rmelt.std(axis=0)
        ax1.fill_between(self.colatitude, rmelt_mean-rmelt_std,
                         rmelt_mean+rmelt_std, alpha=0.3)
        ax1.plot(self.colatitude, rmelt_mean)
        ax1.set_xlabel('Colatitude')
        ax1.set_ylabel('Melt radius')
        ax1.set_xlim(0., np.pi)
        fig1.tight_layout()

def get_rmelt(rad, phase, order=3):
    """
    This subroutine determines the melting radius by doing a spline interpolation
    around phase=0.5

    :param rad: the radius
    :type rad: numpy.ndarray
    :param phase: the phase field as a function of radius
    :type phase: numpy.ndarray
    :param order: order of the interpolation
    :type order: int
    :returns: the melting radius
    :rtype: float
    """

    if order == 1:
        mask = np.where(phase<0.5, 1, 0)
        idx = np.nonzero(mask)[0][0]
        if idx != 1:
            slope = (phase[idx]-phase[idx-1])/(rad[idx]-rad[idx-1])
            intersect = phase[idx]-slope*rad[idx]
            rmelt = (0.5-intersect)/slope
        else:
            rmelt = rad
    else:
        if rad[0] > rad[-1]:
            rr = rad[::-1]
            ph = phase[::-1]
        else:
            rr = rad
            ph = phase

        iph = interp1d(ph, rr)
        rmelt = np.float64(iph(0.5))

    return rmelt

def get_tmelt(rad, temp, rmelt):
    """
    This subroutine determines the melting temperature by doing a spline interpolation
    around the melting radius

    :param rad: the radius
    :type rad: numpy.ndarray
    :param temp: temperature as a function of radius
    :type temp: numpy.ndarray
    :returns: the melting temperature
    :rtype: float
    """

    if rad[0] > rad[-1]:
        rr = rad[::-1]
        tt = temp[::-1]
    else:
        rr = rad
        tt = temp

    ir = interp1d(rr, tt)
    tmelt = np.float64(ir(rmelt))

    return tmelt

def plot_rmelt(gr, levels=17, cm='viridis', order=1):
    """
    :param levels: the number of contour levels
    :type levels: int
    :param cm: the name of the color map
    :type cm: str
    :param order: order of the interpolation
    :type order: int
    :returns: the melting radius as a function of (phi, theta)
    :rtype: numpy.ndarray
    """

    rmelt = np.zeros((gr.n_phi_max, gr.n_theta_max), np.float64)
    for ip in range(gr.n_phi_max):
        for it in range(gr.n_theta_max):
            rmelt[ip, it] = get_rmelt(gr.radius, gr.phase[ip, it, :], order=1)

    rmelt = symmetrize(rmelt, gr.minc)
    radialContour(rmelt, levels=levels, cm=cm, normed=False)

    return rmelt
