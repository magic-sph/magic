import numpy as np
import os
import scipy.interpolate as sint
from magic.libmagic import chebgrid, fd_grid, scanDir

def get_truncation(n_theta_max, nalias, minc):
    """
    This routine determines l_max, m_max and lm_max from the values
    of n_theta_max, minc and nalias.

    :param n_theta_max: number of points along the colatitude
    :type n_theta_max: int
    :param nalias: dealiasing paramete (20 is fully dealiased)
    :type nalias: int
    :param minc: azimuthal symmetry
    :type minc: int
    :returns: returns a list of three integers: l_max, m_max and lm_max
    :rtype: list
    """
    lmax = nalias*n_theta_max // 30
    mmax = (lmax//minc) * minc
    lm_max = mmax*(lmax+1)//minc - \
             mmax*(mmax-minc)//(2*minc)+(lmax+1-mmax)

    return lmax, mmax, lm_max

def get_map(lm_max, lmax, mmax, minc):
    """
    This routine determines the look-up tables to convert the indices
    (l, m) to the single index lm.

    :param lm_max: total number of lm combinations.
    :type lm_max: int
    :param lmax: maximum spherical harmonic degree
    :type lmax: int
    :param mmax: maximum spherical harmonic order
    :type mmax: int
    :param minc: azimuthal symmetry
    :type minc: int
    :returns: returns a list of three look-up tables: idx, lm2l, lm2m
    :rtype: list
    """
    idx = np.zeros((lmax+1, mmax+1), np.int8)
    lm2l = np.zeros(lm_max, np.int8)
    lm2m = np.zeros(lm_max, np.int8)
    idx[0:lmax+2, 0] = np.arange(lmax+1)
    lm2l[0:lmax+1] = np.arange(lmax+1)
    lm2m[0:lmax+2] = 0
    k = lmax+1
    for m in range(minc, lmax+1, minc):
        for l in range(m, lmax+1):
            idx[l, m] = k
            lm2l[k] = l
            lm2m[k] = m
            k +=1

    return idx, lm2l, lm2m

def interp_one_field(field, rold, rnew, rfac=None):
    """
    This routine interpolates a complex input field from an old radial grid
    to a new one. 

    :param field: the field to be interpolated
    :type field: numpy.ndarray
    :param rold: the old radial grid points
    :type rold: numpy.ndarray
    :param rnew: the new radial grid points
    :type rnew: numpy.ndarray
    :param rfac: a rescaling function that depends on the radius
    :type rfac: numpy.ndarray
    :returns: the field interpolated on the new radial grid
    :rtype: numpy.ndarray
    """
    nrold = field.shape[0]
    lm_max = field.shape[1]
    nrnew = rnew.shape[0]

    if rfac is None:
        rfac = np.ones_like(rold)

    if rold[1] > rold[0]:
        old_ordering = 'ascending'
    else:
        old_ordering = 'descending'
        rold = rold[::-1]
        field = field[::-1, :]
        rfac = rfac[::-1]
    if rnew[1] > rnew[0]:
        new_ordering = 'ascending'
    else:
        new_ordering = 'descending'
        rnew = rnew[::-1]

    field_new = np.zeros((nrnew, lm_max), np.complex128)
    for lm in range(lm_max):
        tckp = sint.splrep(rold, rfac*field[:, lm].real, k=5)
        finter_real = sint.splev(rnew, tckp)
        tckp = sint.splrep(rold, rfac*field[:, lm].imag, k=5)
        finter_imag = sint.splev(rnew, tckp)
        field_new[:, lm]  = (finter_real+1j*finter_imag)

    if old_ordering == 'descending':
        rold = rold[::-1]
        field = field[::-1, :]
        rfac = rfac[::-1]
    if new_ordering == 'descending':
        field_new = field_new[::-1, :]

    return field_new

def Graph2Rst(gr, filename='checkpoint_ave'):
    """
    This function allows to transform an input Graphic file into a checkpoint
    file format that can be read by MagIC to restart a simulation.

    >>> # Load a Graphic File
    >>> gr = MagicGraph()
    >>> # Produce the file checkpoint_ave.from_G
    >>> Graph2Rst(gr, filename='checkpoint_ave.from_G')

    :param gr: the input graphic file one wants to convert into a restart file
    :type gr: magic.MagicGraph
    :param filename: name of the checkpoint file
    :type filename: str
    """
    chk = MagicCheckpoint(l_read=False)
    chk.graph2rst(gr, filename)


class MagicCheckpoint:
    """
    This class allows to manipulate checkpoint files produced by MagIC. It
    can read it as

    >>> chk = MagicCheckpoint(filename='checkpoint_end.test')
    >>> print(chk.wpol.shape, chk.l_max)

    This class can also be used to intepolate from FD to Cheb or the opposite
    >>> chk.cheb2fd(96)
    >>> chk.write('checkpoint_fd.test')

    One can also transform a Graphic file into a checkpoint
    >>> gr = MagicGraph()
    >>> chk = MagicCheckpoint(l_read=False)
    >>> chk.graph2rst(gr)

    Finally one can convert checkpoints from XSHELLS
    >>> chk = MagicCheckpoint(l_read=False)
    >>> chk.xshells2magic('st0', 161, rscheme='cheb', cond_state='deltaT')
    """

    def __init__(self, l_read=True, filename=None):
        """
        :param l_read: a boolean to decide whether one reads a checkpoint or not
        :type l_read: bool
        :param filename: name of the checkpoint file to be read
        :type filename: str
        """
        if l_read:
            if filename is None:
                chks = scanDir('checkpoint*')
                filename = chks[-1]
            self.read(filename)

    def read(self, filename):
        """
        This routine is used to read a checkpoint file. 

        :param filename: name of the checkpoint file
        :type filename: str
        """
        file = open(filename, 'rb')

        self.version = np.fromfile(file, dtype=np.int32, count=1)[0]
        self.time = np.fromfile(file, dtype=np.float64, count=1)[0]
        
        # Time scheme
        self.tscheme_family = file.read(10).decode()
        nexp, nimp, nold = np.fromfile(file, dtype=np.int32, count=3)
        if self.tscheme_family.startswith('MULTISTEP'):
            self.dt = np.fromfile(file, dtype=np.float64, count=nexp)
        else:
            self.dt = np.fromfile(file, dtype=np.float64, count=1)[0]
        n_time_step = np.fromfile(file, dtype=np.int32, count=1)[0]
        self.ra, self.pr, self.raxi, self.sc, self.prmag, self.ek, \
                 self.radratio, self.sigma_ratio = \
                 np.fromfile(file, dtype=np.float64, count=8)

        # Truncation
        self.n_r_max, self.n_theta_max, self.n_phi_tot, self.minc,\
                      self.nalias, self.n_r_ic_max = \
                      np.fromfile(file, dtype=np.int32, count=6)
        self.l_max, self.m_max, self.lm_max = get_truncation(self.n_theta_max,
                                                             self.nalias, self.minc)
        # Define maps
        self.idx, self.lm2l, self.lm2m = get_map(self.lm_max, self.l_max,
                                                 self.m_max, self.minc)

        # Radial scheme
        self.rscheme_version = file.read(72).decode()
        if self.rscheme_version.startswith('cheb'):
            self.n_cheb_max, self.map = np.fromfile(file, dtype=np.int32, count=2)
            self.alph1, self.alph2 = np.fromfile(file, dtype=np.float64, count=2)
        else:
            order, obound = np.fromfile(file, dtype=np.int32, count=2)
            self.fd_stretch, self.fd_ratio = \
                                      np.fromfile(file, dtype=np.float64, count=2)

        # Radial grid
        self.radius = np.fromfile(file, dtype=np.float64, count=self.n_r_max)

        # Torques
        if self.tscheme_family.startswith('MULTISTEP'):
            domega_ic = np.fromfile(file, dtype=np.float64, count=(nexp+nimp+nold-3))
            domega_ma = np.fromfile(file, dtype=np.float64, count=(nexp+nimp+nold-3))
            lotorque_ic = np.fromfile(file, dtype=np.float64, count=(nexp+nimp+nold-3))
            lotorque_ma = np.fromfile(file, dtype=np.float64, count=(nexp+nimp+nold-3))

        om = np.fromfile(file, dtype=np.float64, count=12)
        self.omega_ic = om[0]
        self.omega_ma = om[6]

        # Logicals
        self.l_heat, self.l_chem, self.l_mag, self.l_press, self.l_cond_ic = \
            np.fromfile(file, dtype=np.int32, count=5)

        # Fields
        self.wpol = np.fromfile(file, dtype=np.complex128,
                                count=self.n_r_max*self.lm_max)
        self.wpol = self.wpol.reshape((self.n_r_max, self.lm_max))
        if self.tscheme_family.startswith('MULTISTEP'):
            tmp = np.fromfile(file, dtype=np.complex128,
                              count=self.lm_max*self.n_r_max*(nexp+nimp+nold-3))

        self.ztor = np.fromfile(file, dtype=np.complex128,
                                count=self.n_r_max*self.lm_max)
        self.ztor = self.ztor.reshape((self.n_r_max, self.lm_max))
        if self.tscheme_family.startswith('MULTISTEP'):
            tmp = np.fromfile(file, dtype=np.complex128,
                              count=self.lm_max*self.n_r_max*(nexp+nimp+nold-3))

        if self.l_press:
            self.pre = np.fromfile(file, dtype=np.complex128,
                                   count=self.n_r_max*self.lm_max)
            self.pre = self.pre.reshape((self.n_r_max, self.lm_max))
            if self.tscheme_family.startswith('MULTISTEP'):
                tmp = np.fromfile(file, dtype=np.complex128,
                                  count=self.lm_max*self.n_r_max*(nexp+nimp+nold-3))

        if self.l_heat:
            self.entropy = np.fromfile(file, dtype=np.complex128,
                                       count=self.n_r_max*self.lm_max)
            self.entropy = self.entropy.reshape((self.n_r_max, self.lm_max))
            if self.tscheme_family.startswith('MULTISTEP'):
                tmp = np.fromfile(file, dtype=np.complex128,
                                  count=self.lm_max*self.n_r_max*(nexp+nimp+nold-3))
        if self.l_chem:
            self.xi = np.fromfile(file, dtype=np.complex128,
                                  count=self.n_r_max*self.lm_max)
            self.xi = self.xi.reshape((self.n_r_max, self.lm_max))
            if self.tscheme_family.startswith('MULTISTEP'):
                tmp = np.fromfile(file, dtype=np.complex128,
                                  count=self.lm_max*self.n_r_max*(nexp+nimp+nold-3))
        if self.l_mag:
            self.bpol = np.fromfile(file, dtype=np.complex128,
                                    count=self.n_r_max*self.lm_max)
            self.bpol = self.bpol.reshape((self.n_r_max, self.lm_max))
            if self.tscheme_family.startswith('MULTISTEP'):
                tmp = np.fromfile(file, dtype=np.complex128,
                                  count=self.lm_max*self.n_r_max*(nexp+nimp+nold-3))

            self.btor = np.fromfile(file, dtype=np.complex128,
                                    count=self.n_r_max*self.lm_max)
            self.btor = self.btor.reshape((self.n_r_max, self.lm_max))
            if self.tscheme_family.startswith('MULTISTEP'):
                tmp = np.fromfile(file, dtype=np.complex128,
                                  count=self.lm_max*self.n_r_max*(nexp+nimp+nold-3))

        if self.l_cond_ic:
            self.radius_ic = chebgrid(2*self.n_r_ic_max-2, self.radius[-1],
                                      -self.radius[-1])
            self.radius_ic = self.radius_ic[:self.n_r_ic_max]
            self.radius_ic[-1] = 0.
            self.bpol_ic = np.fromfile(file, dtype=np.complex128,
                                       count=self.lm_max*self.n_r_ic_max)
            self.bpol_ic = self.bpol_ic.reshape((self.n_r_ic_max, self.lm_max))
            if self.tscheme_family.startswith('MULTISTEP'):
                tmp = np.fromfile(file, dtype=np.complex128,
                                  count=self.lm_max*self.n_r_ic_max*(nexp+nimp+nold-3))

            self.btor_ic = np.fromfile(file, dtype=np.complex128,
                                       count=self.lm_max*self.n_r_ic_max)
            self.btor_ic = self.btor_ic.reshape((self.n_r_ic_max, self.lm_max))
            if self.tscheme_family.startswith('MULTISTEP'):
                tmp = np.fromfile(file, dtype=np.complex128,
                                  count=self.lm_max*self.n_r_ic_max*(nexp+nimp+nold-3))
        file.close()

    def write(self, filename):
        """
        This routine is used to store a checkpoint file. It only stores the
        state vector not the past quantities required to restart a multistep
        scheme.

        :param filename: name of the checkpoint file
        :type filename: str
        """
        file = open(filename, 'wb')

        # Header
        version = np.array([2], np.int32)
        version.tofile(file)
        time = np.array([self.time], np.float64)
        time.tofile(file)

        # Time scheme
        tscheme = 'DIRK      '.encode()
        file.write(tscheme)
        par = np.array([1, 1, 1], np.int32)
        par.tofile(file)
        if hasattr(self, 'dt'):
            dt = np.array([self.dt], np.float64)
        else:
            dt = np.array([1.0e-6], np.float64)
        dt.tofile(file)
        par = np.array([1], np.int32)
        par.tofile(file)

        # Control parameters
        if hasattr(self, 'ra') and hasattr(self, 'sc') and hasattr(self, 'prmag'):
            x = np.array([self.ra, self.pr, self.raxi, self.sc, self.prmag,
                          self.ek, self.radratio, self.sigma_ratio], np.float64)
        else:
            x = np.array([1e5, 1.0,  0.0, 1.0, 5.0, 1.0e-3, self.radratio, 1.0],
                         np.float64)
        x.tofile(file)

        # Truncation
        x = np.array([self.n_r_max, self.n_theta_max,  self.n_phi_tot, self.minc,
                      self.nalias, self.n_r_ic_max], np.int32)
        x.tofile(file)

        # Radial scheme
        file.write(self.rscheme_version.encode())
        if self.rscheme_version.startswith('cheb'):
            x = np.array([self.n_cheb_max, self.map], np.int32)
            x.tofile(file)
            x = np.array([self.alph1, self.alph2], np.float64)
            x.tofile(file)
        else:
            x = np.array([2, 2], np.int32)
            x.tofile(file)
            x = np.array([self.fd_stretch, self.fd_ratio], np.float64)
            x.tofile(file)

        # Radial grid
        self.radius.tofile(file)

        # torques
        dumm = np.zeros(12, np.float64)
        dumm[0] = self.omega_ic
        dumm[6] = self.omega_ma
        dumm.tofile(file)

        # Logicals
        flags = np.array([self.l_heat, self.l_chem, self.l_mag, False,
                          self.l_cond_ic], np.int32)
        flags.tofile(file)

        # Fields
        self.wpol.tofile(file)
        self.ztor.tofile(file)
        if self.l_heat:
            self.entropy.tofile(file)
        if self.l_chem:
            self.xi.tofile(file)
        if self.l_mag:
            self.bpol.tofile(file)
            self.btor.tofile(file)
        if self.l_cond_ic:
            self.bpol_ic.tofile(file)
            self.btor_ic.tofile(file)

        file.close()

    def cheb2fd(self, n_r_max, fd_stretch=0.3, fd_ratio=0.1):
        """
        This routine is used to convert a checkpoint that has a Gauss-Lobatto
        grid into a finite-difference grid.

        :param n_r_max: number of radial grid points of the finite difference grid
        :type n_r_max: int
        :param fd_stretch: stretching of the radial grid
        :type fd_stretch: float
        :param fd_ratio: ratio of smallest to largest grid spacing
        :type fd_ratio: float
        """
        self.fd_ratio = 0.1
        self.fd_stretch = 0.3
        self.rscheme_version = 'fd'+'{:>70s}'.format('')
        self.fd_stretch = fd_stretch
        self.fd_ratio = fd_ratio
        rnew = fd_grid(n_r_max, self.radius[0], self.radius[-1],
                       self.fd_stretch, self.fd_ratio)

        tmp = interp_one_field(self.wpol, self.radius, rnew)
        self.wpol = tmp
        tmp = interp_one_field(self.ztor, self.radius, rnew)
        self.ztor = tmp
        if self.l_heat:
            tmp = interp_one_field(self.entropy, self.radius, rnew)
            self.entropy = tmp
        if self.l_chem:
            tmp = interp_one_field(self.xi, self.radius, rnew)
            self.xi = tmp
        if self.l_mag:
            tmp = interp_one_field(self.bpol, self.radius, rnew)
            self.bpol = tmp
            tmp = interp_one_field(self.btor, self.radius, rnew)
            self.btor = tmp

        self.radius = rnew
        self.n_r_max = n_r_max

    def fd2cheb(self, n_r_max):
        """
        This routine is used to convert a checkpoint that has finite differences
        in radius into a Gauss-Lobatto grid.

        :param n_r_max: number of radial grid points of the Gauss-Lobatto grid
        :type n_r_max: int
        """
        self.n_cheb_max = n_r_max-2
        self.map = 0
        self.alph1 = 1
        self.alph2 = 0

        self.rscheme_version = 'cheb'+'{:>68s}'.format('')
        rnew = chebgrid(n_r_max-1, self.radius[0], self.radius[-1])

        tmp = interp_one_field(self.wpol, self.radius, rnew)
        self.wpol = tmp
        tmp = interp_one_field(self.ztor, self.radius, rnew)
        self.ztor = tmp
        if self.l_heat:
            tmp = interp_one_field(self.entropy, self.radius, rnew)
            self.entropy = tmp
        if self.l_chem:
            tmp = interp_one_field(self.xi, self.radius, rnew)
            self.xi = tmp
        if self.l_mag:
            tmp = interp_one_field(self.bpol, self.radius, rnew)
            self.bpol = tmp
            tmp = interp_one_field(self.btor, self.radius, rnew)
            self.btor = tmp

        self.radius = rnew
        self.n_r_max = n_r_max

    def xshells2magic(self, xsh_trailing, n_r_max, rscheme='cheb',
                      cond_state='deltaT', scale_b=1.,
                      filename='checkpoint_end.from_xhells'):
        """
        This routine is used to convert XSHELLS field[U,B,T].xsh_trailing files
        into a MagIC checkpoint file.

        >>> chk = MagicCheckPoint()
        >>> # Convert field[U,T,B].st1ns_hr2 into a MagIC checkpoint file
        >>> chk.xshells2magic('st1ns_hr2', 512, rscheme='fd', cond_state='mixed',
                              scale_b=4.472136e-4)

        :param xsh_trailing: trailing of the field[U,B,T].xsh_trailing files
        :type xsh_trailing: str
        :param n_r_max: number of radial grid points to be used
        :type n_r_max: int
        :param rscheme: the type of radial scheme ('cheb' or 'fd')
        :type rscheme: str
        :param cond_state: the type of conducting state:
                              - 'deltaT': fixed temperature contrast
                              - 'mixed': hybrid forcing (STEP1-2 like)
        :type cond_state: str
        :param scale_b: a rescaling factor for the magnetic field
        :type scale_b: float
        """
        import pyxshells as pyx

        # xshells grid
        f = pyx.load_field('fieldU.{}'.format(xsh_trailing))
        if os.path.exists('fieldT.{}'.format(xsh_trailing)):
            self.l_heat = True
        if os.path.exists('fieldB.{}'.format(xsh_trailing)):
            self.l_mag = True
        self.l_press = False
        self.l_chem = False
        # Right now don't know where it is stored
        self.l_cond_ic = False

        rr_xsh = f.grid.r
        nr_xsh = len(rr_xsh)
        ro = rr_xsh[-1]
        ri = rr_xsh[0]
        self.radratio = ri/ro

        self.n_r_max = n_r_max
        if rscheme.startswith('cheb'):
            self.radius = chebgrid(self.n_r_max-1, ro, ri)
            self.n_cheb_max = self.n_r_max-2
            self.map = 0
            self.alph1 = 1.
            self.alph2 = 0.
            self.rscheme_version = 'cheb'+'{:>68s}'.format('')
        else:
            self.radius = fd_grid(n_r_max, ro, ri)
            self.fd_stretch = 0.3
            self.fd_ratio = 0.1
            self.rscheme_version = 'fd'+'{:>70s}'.format('')

        self.l_max = f.lmax
        self.m_max = f.mmax
        self.minc = f.mres
        self.lm_max = f.pol_full().shape[-1]

        # When nalias=60 ntheta=lmax: trick to have lmax in MagIC's header
        self.nalias = 60
        self.n_theta_max = self.l_max 
        self.n_phi_tot = self.l_max
        self.n_r_ic_max = 1

        self.time = f.time
        # Dummy rotation rates: don't know where to get them from xSHELLs
        self.omega_ic = 0.
        self.omega_ma = 0.

        self.wpol = interp_one_field(f.pol_full(), rr_xsh, self.radius,
                                     rfac=rr_xsh)
        self.ztor = interp_one_field(f.tor_full(), rr_xsh, self.radius,
                                     rfac=rr_xsh)

        if self.l_heat:
            f = pyx.load_field('fieldT.{}'.format(xsh_trailing))
            field_xsh = np.zeros((nr_xsh, self.lm_max), np.complex128)
            field_xsh = f.data[1:-1, 0, :]
            self.entropy = interp_one_field(field_xsh, rr_xsh, self.radius)

            if cond_state == 'deltaT':
                temp0 = -ri**2/(ri**2+ro**2)
                tcond = ro*ri/(ro-ri)/self.radius+temp0-ri/(ro-ri)
            elif cond_state == 'mixed':
                fi = 0.75
                ci = (2.*fi-1.)/(ro**3-ri**3)
                co = (fi*ro**3-(1.-fi)*ri**3)/(ro**3-ri**3)
                tcond = ci*self.radius**2/2.+co/self.radius
                tcondo = ci*ro**2/2.+co/ro
                tcond = tcond-tcondo

            self.entropy[:, 0] += np.sqrt(4.*np.pi) * tcond

        if self.l_mag:
            f = pyx.load_field('fieldB.{}'.format(xsh_trailing))
            field_xsh = scale_b * f.pol_full()
            self.bpol = interp_one_field(field_xsh, rr_xsh, self.radius,
                                         rfac=rr_xsh)
            field_xsh = scale_b * f.tor_full()
            self.btor = interp_one_field(field_xsh, rr_xsh, self.radius,
                                         rfac=rr_xsh)

        self.write(filename)

    def graph2rst(self, gr, filename='checkpoint_ave.from_chk'):
        """
        :param gr: the input graphic file one wants to convert into a restart
                   file
        :type gr: magic.MagicGraph
        :param filename: name of the checkpoint file
        :type filename: str
        """

        from magic import SpectralTransforms, thetaderavg, phideravg, MagicRadial

        if hasattr(gr, 'tag'):
            tag = gr.tag

            if os.path.exists('anel.{}'.format(tag)):
                r = MagicRadial(field='anel', iplot=False)
                rho0 = r.rho0
            else:
                rho0 = np.ones_like(gr.radius)
        else:
            rho0 = np.ones_like(gr.radius)

        self.n_r_max = gr.n_r_max
        self.n_theta_max = gr.n_theta_max
        self.n_phi_tot = gr.n_phi_tot
        self.minc = gr.minc
        self.n_r_ic_max = gr.nr_ic+2
        self.nalias = 20

        # Spectral truncation
        self.l_max, self.m_max, self.lm_max = get_truncation(self.n_theta_max,
                                                             self.nalias, self.minc)
        # Define maps
        self.idx, self.lm2l, self.lm2m = get_map(self.lm_max, self.l_max,
                                                 self.m_max, self.minc)
        self.radius = gr.radius.astype(np.float64)
        ri = self.radius[-1]
        ro = self.radius[0]
        self.radratio = ri/ro
        if gr.radial_scheme == 'CHEB':
            self.rscheme_version = 'cheb'+'{:>68s}'.format('')
            self.n_cheb_max = self.n_r_max-2
            if gr.l_newmap == 'F':
                self.map = 0
            else:
                self.map = 1
            self.alph1 = gr.alph1
            self.alph2 = gr.alph2
        else:
            self.rscheme_version = 'fd'+'{:>70s}'.format('')
            self.fd_stretch = gr.fd_stretch
            self.fd_ratio = gr.fd_ratio

        # Flags
        if gr.mode in [2, 3, 7, 8, 9, 10] or gr.ra == 0.:
            self.l_heat = False
        else:
            self.l_heat = True
        if gr.raxi > 0. or gr.raxi < 0.:
            self.l_chem = True
        else:
            self.l_chem = False
        if gr.mode in [0, 2, 3, 6, 8, 9]:
            self.l_mag = True
        else:
            self.l_mag = False
        if gr.sigma_ratio == 0.:
            self.l_cond_ic = False
        else:
            self.l_cond_ic = True
        self.l_press = False

        if self.l_cond_ic:
            self.radius_ic = np.zeros(self.n_r_ic_max, np.float64)
            self.radius_ic[0] = ri
            self.radius_ic[1:] = gr.radius_ic

        self.time = gr.time.astype(np.float64)

        # Rotation rates: dummy
        self.omega_ic = 0.
        self.omega_ma = 0.

        sh = SpectralTransforms(l_max=self.l_max, lm_max=self.lm_max, minc=self.minc,
                                n_theta_max=self.n_theta_max)

        # Calculate and store the poloidal potential using vr
        self.wpol = np.zeros((self.n_r_max, self.lm_max), dtype=np.complex128)
        for i in range(self.n_r_max):
            vr = sh.spat_spec(gr.vr[:, :, i])
            self.wpol[i, 1:] = vr[1:]/(sh.ell[1:]*(sh.ell[1:]+1)) * \
                               self.radius[i]**2 * rho0[i]

        # Calculate the toroidal potential using wr
        self.ztor = np.zeros_like(self.wpol)

        th3D = np.zeros_like(gr.vr)
        rr3D = np.zeros_like(th3D)
        for i in range(self.n_theta_max):
            th3D[:, i, :] = gr.colatitude[i]
        for i in range(self.n_r_max):
            rr3D[:, :, i] = self.radius[i]
        s3D = rr3D*np.sin(th3D)
        omr = 1./s3D*(thetaderavg(np.sin(th3D)*gr.vphi, order=4) -
                      phideravg(gr.vtheta, minc=self.minc))
 
        for i in range(self.n_r_max):
            om = sh.spat_spec(omr[:, :, i])
            self.ztor[i, 1:] = om[1:]/(sh.ell[1:]*(sh.ell[1:]+1)) * \
                               self.radius[i]**2 * rho0[i]

        # Calculate the entropy
        if self.l_heat:
            self.entropy = np.zeros_like(self.wpol)
            for i in range(self.n_r_max):
                p = sh.spat_spec(gr.entropy[:, :, i])
                self.entropy[i, :] = p[:]

        # Calculate the chemical composition
        if self.l_chem:
            self.xi = np.zeros_like(self.wpol)
            for i in range(self.n_r_max):
                p = sh.spat_spec(gr.xi[:, :, i])
                self.xi[i, :] = p[:]

        # Calculate the magnetic field
        if self.l_mag:
            self.bpol = np.zeros_like(self.wpol)
            for i in range(self.n_r_max):
                Br = sh.spat_spec(gr.Br[:, :, i])
                self.bpol[i, 1:] = Br[1:]/(sh.ell[1:]*(sh.ell[1:]+1)) * \
                                   self.radius[i]**2

            self.btor = np.zeros_like(self.ztor)
            jr = 1./s3D*(thetaderavg(np.sin(th3D)*gr.Bphi, order=4) -
                         phideravg(gr.Btheta, minc=self.minc))

            for i in range(self.n_r_max):
                om = sh.spat_spec(jr[:, :, i])
                self.btor[i, 1:] = om[1:]/(sh.ell[1:]*(sh.ell[1:]+1)) * \
                                   self.radius[i]**2

        if self.l_mag and self.l_cond_ic:
            self.bpol_ic = np.zeros((self.n_r_ic_max, self.lm_max), np.complex128)
            for i in range(self.n_r_ic_max):
                rdep = np.ones(sh.ell.shape, dtype=np.float64)
                if i == 0:  # ICB radius
                    vr = sh.spat_spec(gr.Br[:, :, -1])
                    rr = self.radius[-1]
                    rdep[:] = 1.
                else:
                    vr = sh.spat_spec(gr.Br_ic[:, :, i-1])
                    rr = self.radius_ic[i-1]
                    rdep[1:] = (self.radius_ic[i-1]/ri)**(sh.ell[1:]+1)
                self.bpol_ic[i, 1:] = vr[1:]/(sh.ell[1:]*(sh.ell[1:]+1))*rr**2
                # Not stable: 
                #if self.radius_ic[i] >= 0.01:
                #    mask = ( self.lm2l <= 2 ) * ( self.lm2m <= 2)
                #    self.bpol_ic[i, mask] /= rdep[mask]

            # Calculate the toroidal potential using jr
            self.btor_ic = np.zeros_like(self.bpol_ic)

            th3D = np.zeros_like(gr.Br_ic)
            rr3D = np.zeros_like(th3D)
            for i in range(self.n_theta_max):
                th3D[:, i, :] = gr.colatitude[i]
            for i in range(self.n_r_ic_max-1):
                rr3D[:, :, i] = self.radius_ic[i]
            rr3D[:, :, -1] = 1e-4
            s3D = rr3D*np.sin(th3D)
            jr_ic = np.zeros_like(th3D)
            jr_ic = 1./s3D*(thetaderavg(np.sin(th3D)*gr.Bphi_ic, order=4) -
                            phideravg(gr.Btheta_ic, minc=self.minc))

            for i in range(self.n_r_ic_max):
                rdep = np.ones(sh.ell.shape, dtype=np.float64)
                if i == 0:  # ICB radius
                    om = sh.spat_spec(jr[:, :, -1])
                    rr = self.radius[-1]
                    rdep[:] = 1.
                else:
                    om = sh.spat_spec(jr_ic[:, :, i-1])
                    rr = self.radius_ic[i-1]
                    rdep[1:] = (self.radius_ic[i-1]/ri)**(sh.ell[1:]+1)
                self.btor_ic[i, 1:] = om[1:]/(sh.ell[1:]*(sh.ell[1:]+1))*rr**2
                # Not stable
                #if self.radius_ic[i] >= 0.1:
                #    mask = ( self.lm2l <= 5 ) * ( self.lm2m <= 5)
                #    self.btor_ic[i, mask] /= rdep[mask]

        self.write(filename)



if __name__ == '__main__':
    from magic import MagicGraph
    chk = MagicCheckpoint(l_read=False)
    #chk.fd2cheb(33)
    #chk.write('checkpoint_cheb.tmp')
    gr = MagicGraph()
    chk.graph2rst(gr)
    #chk.cheb2fd(96)
    #chk.write('checkpoint_alpha.tmp')
