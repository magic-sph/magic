# -*- coding: utf-8 -*-
import numpy as np
import os
import re
from .log import MagicSetup
from .libmagic import scanDir
from magic.setup import buildSo
from .npfile import npfile

if buildSo:
    try:
        import magic.greader_single as Gsngl
        import magic.greader_double as Gdble
        readingMode = 'f2py'
    except ImportError:
        readingMode = 'python'
    # print('read with {}'.format(readingMode))
else:
    readingMode = 'python'


def getGraphEndianness(filename):
    """
    This function determines the endianness of the graphic files and
    also tries to detect wheter there record markers or not.

    :param filename: input of the filename
    :type filename: str
    :returns: the endianness of the file ('B'='big_endian' or
              'l'='little_endian') and the presence of record
              markers
    :rtype: list
    """
    try:
        f = npfile(filename, endian='B')
        try:
            st = f.fort_read('S20')
            if len(st) > 1:
                raise TypeError
            endian = 'B'
        except TypeError:
            endian = 'l'
        access = 'rm'
    except ValueError:
        f = open(filename, 'rb')
        # Little endian
        version = np.fromfile(f, 'i4', count=1)[0]
        endian = 'l'
        if abs(version) > 100:
            f.close()
            f = open(filename, 'rb')
            version = np.fromfile(f, '>i4', count=1)[0]
            endian = 'B'
        access = 'st'
    f.close()

    return endian, access


class MagicGraph(MagicSetup):
    """
    This class allows to read the 3-D graphic outputs of the MagIC code
    (:ref:`G_#.TAG <secGraphFile>` and G_ave.TAG) files. Those are
    binary unformatted outputs, there are therefore two ways to load them:

       * If buildLib=True in magic.cfg and the fortran libraries were correctly
         built, then the reader uses a fortran program that is expected to be
         much faster than the pure python routine.
       * If buildLib=False, then a pure python program is used to read the G
         files.

    >>> # Regular G files
    >>> gr = MagicGraph(ivar=1, tag='N0m2a')
    >>> print(gr.vr.shape) # shape of vr
    >>> print(gr.ek) # print ekman number
    >>> print(gr.minc) # azimuthal symmetry
    >>> # Averaged G file with double precision
    >>> gr = MagicGraph(ave=True, tag='N0m2', precision=np.float64)
    """

    def __init__(self, ivar=None, datadir='.', quiet=True,
                 ave=False, tag=None, precision=np.float32):
        """
        :param ave: when set to True, it tries to find an average
                    G file (G_ave.TAG)
        :type ave: bool
        :param ivar: the number of the G file
        :type ivar: int
        :param tag: extension TAG of the G file. If not specified, the
                    most recent G_#.TAG file found in the directory will
                    be selected.
        :type tag: str
        :param quiet: when set to True, makes the output silent
        :type quiet: bool
        :param datadir: directory of the G file (default is . )
        :type datadir: str
        :param precision: single or double precision (default np.float32)
        :type precision: str
        """
        self.precision = precision

        if ave:
            self.name = 'G_ave'
        else:
            self.name = 'G_'

        if tag is not None:
            if ivar is not None:
                file = '{}{}.{}'.format(self.name, ivar, tag)
                filename = os.path.join(datadir, file)
            else:
                files = scanDir('{}*{}'.format(self.name, tag))
                if len(files) != 0:
                    filename = os.path.join(datadir, files[-1])
                else:
                    print('No such tag... try again')
                    return

            if os.path.exists('log.{}'.format(tag)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.{}'.format(tag))
        else:
            if ivar is not None:
                files = scanDir('{}{}*'.format(self.name, ivar))
                filename = os.path.join(datadir, files[-1])
            else:
                files = scanDir('{}*'.format(self.name))
                filename = os.path.join(datadir, files[-1])
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists('log.{}'.format(ending)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.{}'.format(ending))

        if not os.path.exists(filename):
            print('No such file')
            return

        if not quiet:
            print('Reading {}'.format(filename))

        # Get file endianness
        endian, access = getGraphEndianness(filename)

        if readingMode != 'python':

            if self.precision == np.float32:
                G = Gsngl.greader_single
            elif self.precision == np.float64:
                G = Gdble.greader_double

            if access == 'st':
                G.readg_stream(filename, endian=endian)
            elif access == 'rm':
                G.readg(filename, endian=endian)
            self.nr = int(G.nr)
            self.n_r_ic_max = int(G.nric)-1
            self.ntheta = int(G.nt)
            self.npI = int(G.np)
            self.minc = int(G.minc)
            self.time = G.time
            self.ra = G.ra
            self.ek = G.ek
            self.pr = G.pr
            self.raxi = G.raxi
            self.sc = G.sc
            self.prmag = G.prmag
            self.radratio = G.radratio
            self.sigma = G.sigma
            if self.npI == self.ntheta*2:
                self.npI = int(self.npI/self.minc)
            self.nphi = int(self.npI*self.minc+1)
            self.radius = G.radius
            self.colatitude = G.colat
            self.entropy = G.entropy
            self.vr = G.vr
            self.vtheta = G.vt
            self.vphi = G.vp
            self.pre = G.pre
            self.xi = G.xi
            if self.prmag != 0:
                self.Br = G.br
                self.Btheta = G.bt
                self.Bphi = G.bp

            if self.prmag != 0 and self.n_r_ic_max > 1:
                self.radius_ic = G.radius_ic
                self.Br_ic = G.br_ic
                self.Btheta_ic = G.bt_ic
                self.Bphi_ic = G.bp_ic
        else:
            if access == 'rm':
                self.read_record_marker(filename, endian, quiet=quiet)
            elif access == 'st':
                self.read_stream(filename, endian)

    def read_stream(self, filename, endian):
        """
        This function is used to read a Graphic file that has no record marker.

        :param filename: name of the graphic file
        :type filename: str
        :param endian: endianness of the file
        :type endian: str
        """
        if endian == 'B':
            prefix = '>'
        else:
            prefix = ''
        if self.precision == np.float32:
            suffix = 4
        else:
            suffix = 8
        f = open(filename, 'rb')

        # Header
        fmt = '{}i{}'.format(prefix, suffix)
        version = np.fromfile(f, fmt, count=1)[0]
        fmt = '{}S64'.format(prefix)
        runID = np.fromfile(f, fmt, count=1)[0]
        fmt = '{}f{}'.format(prefix, suffix)
        self.time = np.fromfile(f, fmt, count=1)[0]
        self.ra, self.pr, self.raxi, self.sc, self.ek, self.prmag, self.radratio, \
            self.sigma = np.fromfile(f, fmt, count=8)
        fmt = '{}i{}'.format(prefix, suffix)
        self.nr, self.ntheta, self.npI, self.minc, self.n_r_ic_max = \
            np.fromfile(f, fmt, count=5)
        if self.npI == self.ntheta*2:
            self.npI = int(self.npI/self.minc)
        self.nphi = self.npI*self.minc+1
        l_heat, l_chem, l_mag, l_press, l_cond_ic = np.fromfile(f, fmt, count=5)
        fmt = '{}f{}'.format(prefix, suffix)
        self.colatitude = np.fromfile(f, fmt, count=self.ntheta)
        self.radius = np.fromfile(f, fmt, count=self.nr)
        if ( l_mag > 0 and self.n_r_ic_max > 1 ):
            self.radius_ic = np.fromfile(f, fmt, count=self.n_r_ic_max)

        self.vr = np.zeros((self.npI, self.ntheta, self.nr), self.precision)
        self.vtheta = np.zeros_like(self.vr)
        self.vphi = np.zeros_like(self.vr)
        if l_heat > 0:
            self.entropy = np.zeros_like(self.vr)
        if l_chem > 0:
            self.xi = np.zeros_like(self.vr)
        if l_press > 0:
            self.pre = np.zeros_like(self.vr)
        if l_mag > 0:
            self.Br = np.zeros_like(self.vr)
            self.Btheta = np.zeros_like(self.vr)
            self.Bphi = np.zeros_like(self.vr)
            if self.n_r_ic_max > 1:
                self.Br_ic  = np.zeros((self.npI, self.ntheta, self.n_r_ic_max), \
                                       self.precision)
                self.Btheta_ic = np.zeros_like(self.Br_ic)
                self.Bphi_ic = np.zeros_like(self.Br_ic)

        # Outer core
        for i in range(self.nr):
            dat = np.fromfile(f, fmt, count=self.ntheta*self.npI)
            self.vr[:, :, i] = dat.reshape(self.npI, self.ntheta)
            dat = np.fromfile(f, fmt, count=self.ntheta*self.npI)
            self.vtheta[:, :, i] = dat.reshape(self.npI, self.ntheta)
            dat = np.fromfile(f, fmt, count=self.ntheta*self.npI)
            self.vphi[:, :, i] = dat.reshape(self.npI, self.ntheta)
            if l_heat > 0:
                dat = np.fromfile(f, fmt, count=self.ntheta*self.npI)
                self.entropy[:, :, i] = dat.reshape(self.npI, self.ntheta)
            if l_chem > 0:
                dat = np.fromfile(f, fmt, count=self.ntheta*self.npI)
                self.xi[:, :, i] = dat.reshape(self.npI, self.ntheta)
            if l_press > 0:
                dat = np.fromfile(f, fmt, count=self.ntheta*self.npI)
                self.pre[:, :, i] = dat.reshape(self.npI, self.ntheta)
            if l_mag > 0:
                dat = np.fromfile(f, fmt, count=self.ntheta*self.npI)
                self.Br[:, :, i] = dat.reshape(self.npI, self.ntheta)
                dat = np.fromfile(f, fmt, count=self.ntheta*self.npI)
                self.Btheta[:, :, i] = dat.reshape(self.npI, self.ntheta)
                dat = np.fromfile(f, fmt, count=self.ntheta*self.npI)
                self.Bphi[:, :, i] = dat.reshape(self.npI, self.ntheta)

        # Inner core
        if ( l_mag > 0 and self.n_r_ic_max > 1 ):
            for i in range(self.n_r_ic_max):
                dat = np.fromfile(f, fmt, count=self.ntheta*self.npI)
                self.Br_ic[:, :, i] = dat.reshape(self.npI, self.ntheta)
                dat = np.fromfile(f, fmt, count=self.ntheta*self.npI)
                self.Btheta_ic[:, :, i] = dat.reshape(self.npI, self.ntheta)
                dat = np.fromfile(f, fmt, count=self.ntheta*self.npI)
                self.Bphi_ic[:, :, i] = dat.reshape(self.npI, self.ntheta)

        f.close()

    def read_record_marker(self, filename, endian, quiet=True):
        """
        This function is used to read a Graphic file that contains record markers.

        :param filename: name of the graphic file
        :type filename: str
        :param endian: endianness of the file
        :type endian: str
        :param quiet: when set to True, makes the output silent
        :type quiet: bool
        """
        # Read data
        inline = npfile(filename, endian=endian)

        # Read the header
        version = inline.fort_read('|S20')[0]
        version = version.rstrip()
        runID = inline.fort_read('|S64')[0]
        self.time, n_r_max, n_theta_max, self.npI, n_r_ic_max, minc, \
            nThetaBs, self.ra, self.ek, self.pr, self.prmag, \
            self.radratio, self.sigma = inline.fort_read(self.precision)

        self.nr = int(n_r_max)
        self.ntheta = int(n_theta_max)
        self.npI = int(self.npI)
        self.minc = int(minc)
        if self.npI == self.ntheta*2:
            self.npI = int(self.npI/self.minc)
        self.nphi = self.npI*self.minc+1
        self.n_r_ic_max = int(n_r_ic_max) - 1
        nThetaBs = int(nThetaBs)

        if not quiet:
            st = 'Rayleigh = {:.1e}, Ekman = {:.1e}, Prandtl = {:.1e}'.format(
                self.ra, self.ek, self.pr)
            print(st)
            print('nr = {}, nth = {}, nphi = {}'.format(
                  self.nr, self.ntheta, self.npI))

        self.colatitude = inline.fort_read(self.precision)
        self.radius = np.zeros((self.nr), self.precision)

        if self.prmag != 0 and self.n_r_ic_max > 1:
            self.radius_ic = np.zeros((self.n_r_ic_max), self.precision)

        entropy = np.zeros((self.npI, self.ntheta, self.nr),
                           self.precision)
        vr = np.zeros_like(entropy)
        vtheta = np.zeros_like(entropy)
        vphi = np.zeros_like(entropy)
        if self.prmag != 0:
            Br = np.zeros_like(entropy)
            Btheta = np.zeros_like(entropy)
            Bphi = np.zeros_like(entropy)
            if self.sigma != 0:
                Br_ic = np.zeros((self.npI, self.ntheta, self.n_r_ic_max),
                                 self.precision)
                Btheta_ic = np.zeros_like(Br_ic)
                Bphi_ic = np.zeros_like(Br_ic)
        if version == b'Graphout_Version_8' or \
           version == b'Graphout_Version_10' or \
           version == b'Graphout_Version_12':
            pressure = np.zeros_like(entropy)
        if version == b'Graphout_Version_11' or \
           version == b'Graphout_Version_12':
            xi = np.zeros_like(entropy)

        for k in range(self.nr*nThetaBs):
            # radius and Thetas in this block
            ir, rad, ilat1, ilat2 = inline.fort_read(self.precision)
            ir = int(ir)
            ilat1 = int(ilat1) - 1
            ilat2 = int(ilat2) - 1
            self.radius[ir] = rad
            nth_loc = ilat2 - ilat1 + 1

            if version == b'Graphout_Version_9' or \
               version == b'Graphout_Version_10' or \
               version == b'Graphout_Version_11' or \
               version == b'Graphout_Version_12':
                if self.prmag != 0:
                    data = inline.fort_read(self.precision,
                                            shape=(nth_loc, self.npI))
                    entropy[:, ilat1:ilat2+1, ir] = data.T
                    data = inline.fort_read(self.precision,
                                            shape=(nth_loc, self.npI))
                    vr[:, ilat1:ilat2+1, ir] = data.T
                    data = inline.fort_read(self.precision,
                                            shape=(nth_loc, self.npI))
                    vtheta[:, ilat1:ilat2+1, ir] = data.T
                    data = inline.fort_read(self.precision,
                                            shape=(nth_loc, self.npI))
                    vphi[:, ilat1:ilat2+1, ir] = data.T
                    if version == b'Graphout_Version_11' or \
                       version == b'Graphout_Version_12':
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        xi[:, ilat1:ilat2+1, ir] = data.T
                    if version == b'Graphout_Version_10' or \
                       version == b'Graphout_Version_12':
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        pressure[:, ilat1:ilat2+1, ir] = data.T
                    data = inline.fort_read(self.precision,
                                            shape=(nth_loc, self.npI))
                    Br[:, ilat1:ilat2+1, ir] = data.T
                    data = inline.fort_read(self.precision,
                                            shape=(nth_loc, self.npI))
                    Btheta[:, ilat1:ilat2+1, ir] = data.T
                    data = inline.fort_read(self.precision,
                                            shape=(nth_loc, self.npI))
                    Bphi[:, ilat1:ilat2+1, ir] = data.T

                else:
                    # vectorize
                    data = inline.fort_read(self.precision,
                                            shape=(nth_loc, self.npI))
                    entropy[:, ilat1:ilat2+1, ir] = data.T
                    data = inline.fort_read(self.precision,
                                            shape=(nth_loc, self.npI))
                    vr[:, ilat1:ilat2+1, ir] = data.T
                    data = inline.fort_read(self.precision,
                                            shape=(nth_loc, self.npI))
                    vtheta[:, ilat1:ilat2+1, ir] = data.T
                    data = inline.fort_read(self.precision,
                                            shape=(nth_loc, self.npI))
                    vphi[:, ilat1:ilat2+1, ir] = data.T
                    if version == b'Graphout_Version_11' or \
                       version == b'Graphout_Version_12':
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        xi[:, ilat1:ilat2+1, ir] = data.T
                    if version == b'Graphout_Version_10' or \
                       version == b'Graphout_Version_12':
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        pressure[:, ilat1:ilat2+1, ir] = data.T
            else:
                if self.prmag != 0:
                    # To vectorize, one must also read the end-of-line
                    # symbol that accounts for 2 additionnal floats
                    # For the last line of the chunk, we should not
                    # read the last symbol
                    if version == b'Graphout_Version_8':
                        data = inline.fort_read(
                            self.precision,
                            shape=(8*(self.npI+2)*nth_loc-2))
                        # Add 2 zeros for the last lines
                        data = np.append(data, 0.)
                        data = np.append(data, 0.)
                        data = data.reshape((8, nth_loc, self.npI+2))
                    else:
                        data = inline.fort_read(self.precision,
                                                shape=(7*(self.npI+2)*nth_loc-2))
                        # Add 2 zeros for the last lines
                        data = np.append(data, 0.)
                        data = np.append(data, 0.)
                        data = data.reshape((7, nth_loc, self.npI+2))
                    data = data[:, :, :-2:]
                    entropy[:, ilat1:ilat2+1, ir] = data[0, ...].T
                    vr[:, ilat1:ilat2+1, ir] = data[1, ...].T
                    vtheta[:, ilat1:ilat2+1, ir] = data[2, ...].T
                    vphi[:, ilat1:ilat2+1, ir] = data[3, ...].T
                    if version == b'Graphout_Version_8':
                        pressure[:, ilat1:ilat2+1, ir] = data[4, ...].T
                        Br[:, ilat1:ilat2+1, ir] = data[5, ...].T
                        Btheta[:, ilat1:ilat2+1, ir] = data[6, ...].T
                        Bphi[:, ilat1:ilat2+1, ir] = data[7, ...].T
                    else:
                        Br[:, ilat1:ilat2+1, ir] = data[4, ...].T
                        Btheta[:, ilat1:ilat2+1, ir] = data[5, ...].T
                        Bphi[:, ilat1:ilat2+1, ir] = data[6, ...].T
                else:
                    if version == b'Graphout_Version_8':
                        data = inline.fort_read(
                            self.precision,
                            shape=(5*(self.npI+2)*nth_loc-2))
                        # Add 2 zeros for the last lines
                        data = np.append(data, 0.)
                        data = np.append(data, 0.)
                        data = data.reshape((5, nth_loc, self.npI+2))
                    else:
                        data = inline.fort_read(
                            self.precision,
                            shape=(4*(self.npI+2)*nth_loc-2))
                        # Add 2 zeros for the last lines
                        data = np.append(data, 0.)
                        data = np.append(data, 0.)
                        data = data.reshape((4, nth_loc, self.npI+2))
                    data = data[:, :, :-2:]
                    entropy[:, ilat1:ilat2+1, ir] = data[0, ...].T
                    vr[:, ilat1:ilat2+1, ir] = data[1, ...].T
                    vtheta[:, ilat1:ilat2+1, ir] = data[2, ...].T
                    vphi[:, ilat1:ilat2+1, ir] = data[3, ...].T
                    if version == b'Graphout_Version_8':
                        pressure[:, ilat1:ilat2+1, ir] = data[4, ...].T

        if self.prmag != 0 and self.n_r_ic_max > 1:
            for k in range(self.n_r_ic_max-1):
                # radius and Thetas in this block
                ir, rad, ilat1, ilat2 = inline.fort_read(self.precision)
                ir = ir-self.nr+1
                ir = int(ir)
                ilat1 = int(ilat1) - 1
                ilat2 = int(ilat2) - 1
                self.radius_ic[ir] = rad
                nth_loc = ilat2 - ilat1 + 1

                if version == b'Graphout_Version_9' or \
                   version == b'Graphout_Version_10' or \
                   version == b'Graphout_Version_11' or \
                   version == b'Graphout_Version_12':

                    data = inline.fort_read(self.precision,
                                            shape=(nth_loc, self.npI))
                    Br_ic[:, ilat1:ilat2+1, ir] = data.T
                    data = inline.fort_read(self.precision,
                                            shape=(nth_loc, self.npI))
                    Btheta_ic[:, ilat1:ilat2+1, ir] = data.T
                    data = inline.fort_read(self.precision,
                                            shape=(nth_loc, self.npI))
                    Bphi_ic[:, ilat1:ilat2+1, ir] = data.T

                else:
                    data = inline.fort_read(self.precision,
                                            shape=(3*(self.npI+2)*nth_loc-2))
                    data = np.append(data, 0.)
                    data = np.append(data, 0.)
                    data = data.reshape((3, nth_loc, self.npI+2))
                    data = data[:, :, :-2:]
                    Br_ic[:, ilat1:ilat2+1, ir] = data[0, ...].T
                    Btheta_ic[:, ilat1:ilat2+1, ir] = data[1, ...].T
                    Bphi_ic[:, ilat1:ilat2+1, ir] = data[2, ...].T

            Br_ic[:, :, 0] = Br[:, :, -1]
            Bt_ic[:, :, 0] = Bt[:, :, -1]
            Bp_ic[:, :, 0] = Bp[:, :, -1]

        inline.close()
        # Sorting of data (strange hemispherical way that the magic
        # code uses
        entropy = self.rearangeLat(entropy)
        vr = self.rearangeLat(vr)
        vtheta = self.rearangeLat(vtheta)
        vphi = self.rearangeLat(vphi)
        if self.prmag != 0:
            Br = self.rearangeLat(Br)
            Btheta = self.rearangeLat(Btheta)
            Bphi = self.rearangeLat(Bphi)
            if self.n_r_ic_max > 1:
                Br_ic = self.rearangeLat(Br_ic)
                Btheta_ic = self.rearangeLat(Btheta_ic)
                Bphi_ic = self.rearangeLat(Bphi_ic)
        if version == b'Graphout_Version_11' or \
           version == b'Graphout_Version_12':
            xi = self.rearangeLat(xi)
        if version == b'Graphout_Version_8' or \
           version == b'Graphout_Version_10' or \
           version == b'Graphout_Version_12':
            pressure = self.rearangeLat(pressure)

        # Normalise r
        self.radius = self.radius/(1.-self.radratio)
        if self.prmag != 0 and self.n_r_ic_max > 1:
            self.radius_ic = self.radius_ic/(1.-self.radratio)

        # Full solution by repeating minc times the structure
        self.nphi = self.minc * self.npI + 1
        self.entropy = entropy
        self.vr = vr
        self.vtheta = vtheta
        self.vphi = vphi
        if self.prmag != 0:
            self.Br = Br
            self.Btheta = Btheta
            self.Bphi = Bphi
            if self.n_r_ic_max > 1:
                self.Br_ic = Br_ic
                self.Btheta_ic = Btheta_ic
                self.Bphi_ic = Bphi_ic

        if version == b'Graphout_Version_11' or \
           version == b'Graphout_Version_12':
            self.xi = xi
        if version == b'Graphout_Version_8' or \
           version == b'Graphout_Version_10' or \
           version == b'Graphout_Version_12':
            self.pre = pressure

    def rearangeLat(self, field):
        """
        This function is used to unfold the colatitudes

        :param field: input array with MagIC ordering of colatitudes (i.e.
                      successively Northern Hemisphere and Southern
                      Hemisphere)
        :type field: numpy.ndarray
        :return: an array with the regular ordering of the colatitudes
        :rtype: numpy.ndarray
        """
        even = field[:, ::2, :]
        odd = field[:, 1::2, :]

        return np.concatenate((even, odd[:, ::-1, :]), axis=1)
