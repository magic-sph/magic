# -*- coding: utf-8 -*-
import scipy.interpolate as S
import numpy as np
import glob, os, re, sys
from .npfile import *


def selectField(obj, field, labTex=True, ic=False):
    """
    This function selects for you which field you want to display. It actually
    allows to avoid possible variables miss-spelling: i.e. 'Bphi'='bp'='Bp'='bphi'

    :param obj: a graphic output file
    :type obj: :py:class:`magic.MagicGraph`
    :param field: the name of the field one wants to select
    :type field: str
    :param labTex: when set to True, format the labels using LaTeX fonts
    :type labTex: bool
    :returns: a tuple that contains the selected physical field and its label
    :rtype: (numpy.ndarray, str)
    """
    if field in ('Bp', 'bp', 'bphi', 'Bphi'):
        data = obj.Bphi
        if labTex:
            label = r'$B_{\phi}$'
        else:
            label = 'Bphi'
        if ic:
            data_ic = obj.Bphi_ic
        else:
            data_ic = None
    elif field in ('Bt', 'bt', 'btheta', 'Btheta'):
        data = obj.Btheta
        if labTex:
            label = r'$B_{\theta}$'
        else:
            label = 'Btheta'
        if ic:
            data_ic = obj.Btheta_ic
        else:
            data_ic = None
    elif field in ('Br', 'br'):
        data = obj.Br
        if labTex:
            label = r'$B_r$'
        else:
            label = 'Br'
        if ic:
            data_ic = obj.Br_ic
        else:
            data_ic = None
    elif field in ('pressure', 'pre', 'Pre', 'Pressure', 'press', 'Press'):
        data = obj.pre
        if labTex:
            label = r'$p$'
        else:
            label = 'p'
        data_ic = None
    elif field in ('composition', 'xi', 'Xi', 'Comp', 'comp', 'chem', 'Chem'):
        data = obj.xi
        if labTex:
            label = r'$\xi$'
        else:
            label = 'xi'
        data_ic = None
    elif field in ('Vr', 'vr', 'Ur', 'ur'):
        data = obj.vr
        if labTex:
            label = r'$v_r$'
        else:
            label = 'vr'
        data_ic = None
    elif field in ('Vtheta', 'vtheta', 'Utheta', 'utheta', 'vt', 'Vt',
                   'Ut', 'ut'):
        data = obj.vtheta
        if labTex:
            label = r'$v_{\theta}$'
        else:
            label = 'vtheta'
        data_ic = None
    elif field in ('Vphi', 'vphi', 'Uphi', 'uphi', 'up', 'Up', 'Vp', 'vp'):
        data = obj.vphi
        if labTex:
            label = r'$v_{\phi}$'
        else:
            label = 'vphi'
        data_ic = None
    elif field in ('entropy', 's', 'S'):
        data = obj.entropy
        label = 'Entropy'
        data_ic = None
    elif field in ('u2'):
        data = obj.vphi**2+obj.vr**2+obj.vtheta**2
        if labTex:
            label = r'$u^2$'
        else:
            label = 'u2'
        data_ic = None
    elif field in ('nrj'):
        temp0, rho0, beta = anelprof(obj.radius, obj.strat, obj.polind)
        data = 1./2.*rho0*(obj.vphi**2+obj.vr**2+obj.vtheta**2)
        if labTex:
            label = r'$E_{\hbox{kin}}$'
        else:
            label = 'Ekin'
        data_ic = None
    elif field in ('b2', 'B2'):
        data = obj.Bphi**2+obj.Br**2+obj.Btheta**2
        if labTex:
            label = r'$B^2$'
        else:
            label = 'B2'
        if ic:
            data_ic = obj.Bphi_ic**2+obj.Br_ic**2+obj.Btheta_ic**2
        else:
            data_ic = None
    elif field in ('vrconv', 'vrc'):
        data = obj.vr-obj.vr.mean(axis=0)
        if labTex:
            label = r'$v_{r}$ conv'
        else:
            label = 'vr conv'
        data_ic = None
    elif field in ('vtconv', 'vtc'):
        data = obj.vtheta-obj.vtheta.mean(axis=0)
        if labTex:
            label = r'$v_{\theta}$ conv'
        else:
            label = 'vt conv'
        data_ic = None
    elif field in ('vpconv', 'vpc'):
        data = obj.vphi-obj.vphi.mean(axis=0)
        if labTex:
            label = r'$v_{\phi}$ conv'
        else:
            label = 'vp conv'
        data_ic = None
    elif field in ('bpfluct'):
        data = obj.Bphi-obj.Bphi.mean(axis=0)
        if labTex:
            label = r"$B_{\phi}'$"
        else:
            label = "Bp'"
        if ic:
            data_ic = obj.Bphi_ic-obj.Bphi_ic.mean(axis=0)
        else:
            data_ic = None
    elif field in ('brfluct'):
        data = obj.Br-obj.Br.mean(axis=0)
        if labTex:
            label = r"$B_r'$"
        else:
            label = "Br'"
        if ic:
            data_ic = obj.Br_ic-obj.Br_ic.mean(axis=0)
        else:
            data_ic = None
    elif field in ('entropyfluct'):
        data = obj.entropy-obj.entropy.mean(axis=0)
        if labTex:
            label = r"$s'$"
        else:
            label = "s'"
        data_ic = None
    elif field in ('xifluct'):
        data = obj.xi-obj.xi.mean(axis=0)
        if labTex:
            label = r"$\xi'$"
        else:
            label = "xi'"
        data_ic = None
    elif field in ('prefluct'):
        data = obj.pre-obj.pre.mean(axis=0)
        if labTex:
            label = r"$p'$"
        else:
            label = "p'"
        data_ic = None
    elif field in ('vrea'):
        data = np.zeros_like(obj.vr)
        for i in range(obj.ntheta):
            data[:, i, :] = (obj.vr[:, i, :]-obj.vr[:, -i-1, :])/2.
        if labTex:
            label = r'$v_{r}$ ea'
        else:
            label = 'vr ea'
        data_ic = None
    elif field in ('vra'):
        data = np.zeros_like(obj.vr)
        for i in range(obj.ntheta):
            data[:, i, :] = (obj.vr[:, i, :]+obj.vr[:, -i-1, :])/2.
        if labTex:
            label = r'$v_{r}$ es'
        else:
            label = 'vr es'
        data_ic = None
    elif field in ('vpea'):
        data = np.zeros_like(obj.vr)
        for i in range(obj.ntheta):
            data[:, i, :] = (obj.vphi[:, i, :]-obj.vphi[:, -i-1, :])/2.
        if labTex:
            label = r'$v_{\phi}$ ea'
        else:
            label = r'vp ea'
        data_ic = None
    elif field in ('vpa'):
        data = np.zeros_like(obj.vr)
        for i in range(obj.ntheta):
            data[:, i, :] = (obj.vphi[:, i, :]+obj.vphi[:, -i-1, :])/2.
        if labTex:
            label = r'$v_{\phi}$ es'
        else:
            label = r'vp es'
        data_ic = None
    elif field in ('tea'):
        data = np.zeros_like(obj.vr)
        for i in range(obj.ntheta):
            data[:, i, :] = (obj.entropy[:, i, :]-obj.entropy[:, -i-1, :])/2.
        if labTex:
            label = r'$s$ ea'
        else:
            label = r's ea'
        data_ic = None

    return data, data_ic, label

def avgField(time, field, tstart=None, std=False):
    """
    This subroutine computes the time-average (and the std) of a time series

    >>> ts = MagicTs(field='misc', iplot=False, all=True)
    >>> nuavg = avgField(ts.time, ts.topnuss, 0.35)
    >>> print(nuavg)

    :param time: time
    :type time: numpy.ndarray
    :param field: the time series of a given field
    :type field: numpy.ndarray
    :param tstart: the starting time of the averaging
    :type tstart: float
    :param std: when set to True, the standard deviation is also calculated
    :type std: bool
    :returns: the time-averaged quantity
    :rtype: float
    """
    if tstart is not None:
        mask = np.where(abs(time-tstart) == min(abs(time-tstart)), 1, 0)
        ind = np.nonzero(mask)[0][0]
    else: # the whole input array is taken!
        ind = 0
    fac = 1./(time[-1]-time[ind])
    avgField = fac*np.trapz(field[ind:], time[ind:])

    if std:
        stdField = np.sqrt(fac*np.trapz((field[ind:]-avgField)**2, time[ind:]))
        return avgField, stdField
    else:
        return avgField

def writeVpEq(par, tstart):
    """
    This function computes the time-averaged surface zonal flow (and Rolc) and
    format the output

    >>> # Reads all the par.* files from the current directory
    >>> par = MagicTs(field='par', iplot=False, all=True)
    >>> # Time-average
    >>> st = writeVpEq(par, tstart=2.1)
    >>> print(st)

    :param par: a :py:class:`MagicTs <magic.MagicTs>` object containing the par file
    :type par: :py:class:`magic.MagicTs`
    :param tstart: the starting time of the averaging
    :type tstart: float
    :returns: a formatted string
    :rtype: str
    """
    mask = np.where(abs(par.time-tstart) == min(abs(par.time-tstart)), 1, 0)
    ind = np.nonzero(mask)[0][0]
    fac = 1./(par.time.max()-par.time[ind])
    avgReEq = fac*np.trapz(par.reEquat[ind:], par.time[ind:])
    roEq = avgReEq*par.ek*(1.-par.radratio)
    avgRolC = fac*np.trapz(par.rolc[ind:], par.time[ind:])
    st = '{:10.3e}{:5.2f}{:6.2f}{:11.3e}{:11.3e}{:11.3e}'.format(par.ek,
        par.strat, par.pr, par.ra, roEq, avgRolC)

    return st

def progressbar(it, prefix="", size=60):
    """
    Fancy progress-bar for loops

    .. code-block:: python

           for i in progressbar(range(1000000)):
               x = i

    :type it: iterator
    :param prefix: prefix string before progress bar
    :type prefix: str
    :param size: width of the progress bar (in points of xterm width)
    :type size: int
    :type size: int
    """
    count = len(it)
    def _show(_i):
        x = int(size*_i/count)
        sys.stdout.write("{}[{}{}] {}/{}\r".format(prefix, "#"*x, "."*(size-x),
                                                   _i, count))
        sys.stdout.flush()

    _show(0)
    for i, item in enumerate(it):
        yield item
        _show(i+1)
    sys.stdout.write("\n")
    sys.stdout.flush()

def scanDir(pattern, tfix=None):
    """
    This function sorts the files which match a given input pattern from the oldest
    to the most recent one (in the current working directory)

    >>> dat = scanDir('log.*')
    >>> print(log)

    :param pattern: a classical regexp pattern
    :type pattern: str
    :param tfix: in case you want to add only the files that are more recent than
                 a certain date, use tfix (computer 1970 format!!)
    :type tfix: float
    :returns: a list of files that match the input pattern
    :rtype: list
    """
    dat = [(os.stat(i).st_mtime, i) for i in glob.glob(pattern)]
    dat.sort()
    if tfix is not None:
        out = []
        for i in dat:
            if i[0] > tfix:
                out.append(i[1])
    else:
        out = [i[1] for i in dat]
    return out

def symmetrize(data, ms, reversed=False):
    """
    Symmetrise an array which is defined only with an azimuthal symmetry minc=ms

    :param data: the input array
    :type data: numpy.ndarray
    :param ms: the azimuthal symmetry
    :type ms: int
    :param reversed: set to True, in case the array is reversed (i.e. n_phi is the last column)
    :type reversed: bool
    :returns: an output array of dimension (data.shape[0]*ms+1)
    :rtype: numpy.ndarray
    """
    if reversed:
        nphi = data.shape[-1]*ms+1
        size = [nphi]
        size.insert(0,data.shape[-2])
        if len(data.shape) == 3:
            size.insert(0,data.shape[-3])
        out = np.zeros(size, dtype=data.dtype)
        for i in range(ms):
            out[..., i*data.shape[-1]:(i+1)*data.shape[-1]] = data
        out[..., -1] = out[..., 0]
    else:
        nphi = data.shape[0]*ms +1
        size = [nphi]
        if len(data.shape) >= 2:
            size.append(data.shape[1])
        if len(data.shape) == 3:
            size.append(data.shape[2])
        out = np.zeros(size, dtype=data.dtype)
        for i in range(ms):
            out[i*data.shape[0]:(i+1)*data.shape[0], ...] = data
        out[-1, ...] = out[0, ...]
    return out

def fast_read(file, skiplines=0, binary=False, precision=np.float64):
    """
    This function reads an input ascii table
    (can read both formatted or unformatted fortran)

    >>> # Read 'e_kin.test', skip the first 10 lines
    >>> data = fast_read('e_kin.test', skiplines=10)

    :param file: name of the input file
    :type file: str
    :param skiplines: number of header lines to be skept during reading
    :type skiplines: int
    :param binary: when set to True, try to read an unformatted binray Fortran
                   file (default is False)
    :type binary: bool
    :param precision: single (np.float32) or double precision (np.float64)
    :type precision: str
    :returns: an array[nlines, ncols] that contains the data of the ascii file
    :rtype: numpy.ndarray
    """
    if not binary:
        f = open(file, 'r')
        X = []
        for k, line in enumerate(f.readlines()):
            st = line.replace('D', 'E')
            if k >= skiplines:
                X.append(st.split())
        X = np.array(X, dtype=precision)
        f.close()
    else:
        f = npfile(file, endian='B')
        X = []
        while 1:
            try:
                X.append(f.fort_read(precision))
            except TypeError:
                break
        X = np.array(X, dtype=precision)
        f.close()
    return X



def anelprof(radius, strat, polind, g0=0., g1=0., g2=1.):
    """
    This functions calculates the reference temperature and density profiles
    of an anelastic model.

    >>> rad = chebgrid(65, 1.5, 2.5)
    >>> temp, rho, beta = anelprof(rad, strat=5., polind=2.)

    :param radius: the radial gridpoints
    :type radius: numpy.ndarray
    :param polind: the polytropic index
    :type polind: float
    :param strat: the number of the density scale heights between the inner
                  and the outer boundary
    :type strat: float
    :param g0: gravity profile: g=g0
    :type g0: float
    :param g1: gravity profile: g=g1*r/r_o
    :type g1: float
    :param g2: gravity profile: g=g2*(r_o/r)**2
    :type g2: float
    :returns: a tuple that contains the temperature profile, the density profile
              and the log-derivative of the density profile versus radius
    :rtype: (numpy.ndarray, numpy.ndarray, numpy.ndarray)
    """
    if radius[-1] < radius[0]:
        ro = radius[0]
        radratio = radius[-1]/radius[0]
    else:
        ro = radius[-1]
        radratio = radius[0]/radius[-1]
    grav = g0 + g1 * radius/ro + g2 * (ro/radius)**2
    ofr=( np.exp(strat/polind)-1. )/ ( g0+0.5*g1*(1.+radratio) +g2/radratio )

    temp0 = -ofr*( g0*radius +0.5*g1*radius**2/ro-g2*ro**2/radius ) + \
            1.+ ofr*ro*(g0+0.5*g1-g2)
    rho0 = temp0**polind
    beta = -ofr*grav*polind/temp0
    return temp0, rho0, beta

def fd_grid(nr, a, b, fd_stretching=0.3, fd_ratio=0.1):
    """
    This function defines a stretched grid between a and b

    >>> r_icb = 0.5 ; r_cmb = 1.5; n_r_max=64
    >>> rr = fd_grid(n_r_max, r_cmb, r_icb)

    :param nr: number of radial grid points
    :type nr: int
    :param a: upper boundary of the grid
    :type a: float
    :param b: lower boundary of the grid
    :type b: float
    :param fd_stretching: fraction of points in the bulk
    :type fd_stretching: float
    :param fd_ratio: ratio of minimum to maximum spacing
    :type fd_ratio: float
    :returns: the radial grid
    :returns: the radial grid
    :rtype: numpy.ndarray
    """

    ratio1 = fd_stretching
    ratio2 = fd_ratio

    if abs(a-b-1.0) > 1.0e-12:
        sys.exit('Not implemented yet')

    rr = np.zeros(nr, dtype=np.float64)
    rr[0] = a

    if ratio2 == 1.0: # Regular grid
        dr_before = (a-b)/(float(nr)-1.)
        dr_after = dr_before
        for i in range(1, nr):
            rr[i]=rr[i-1]-dr_before
    
    else:
        n_boundary_points = int( float(nr-1)/(2.*(1+ratio1)) )
        ratio1 = float(nr-1)/float(2.*n_boundary_points) -1.

        n_bulk_points = nr-1-2*n_boundary_points

        dr_after = np.exp(np.log(ratio2)/float(n_boundary_points))
        dr_before = 1.
        for i in range(n_boundary_points):
            dr_before = dr_before*dr_after
        dr_before = 1./(float(n_bulk_points)+ \
                    2.*dr_after*((1-dr_before)/(1.-dr_after)))
    
        for i in range(n_boundary_points):
            dr_before = dr_before*dr_after

        for i in range(1, n_boundary_points+1):
            rr[i] = rr[i-1]-dr_before
            dr_before = dr_before/dr_after

        for i in range(n_bulk_points):
            rr[n_boundary_points+1+i] = rr[n_boundary_points+i]-dr_before

        for i in range(n_boundary_points):
            dr_before = dr_before*dr_after
            rr[n_boundary_points+1+n_bulk_points+i] = \
                rr[n_boundary_points+n_bulk_points+i]-dr_before

    return rr

def chebgrid(nr, a, b):
    """
    This function defines a Gauss-Lobatto grid from a to b.

    >>> r_icb = 0.5 ; r_cmb = 1.5; n_r_max=65
    >>> rr = chebgrid(n_r_max, r_icb, r_cmb)

    :param nr: number of radial grid points plus one (Nr+1)
    :type nr: int
    :param a: lower limit of the Gauss-Lobatto grid
    :type a: float
    :param b: upper limit of the Gauss-Lobatto grid
    :type b: float
    :returns: the Gauss-Lobatto grid
    :rtype: numpy.ndarray
    """
    rst = (a+b)/(b-a)
    rr = 0.5*(rst+np.cos(np.pi*(1.-np.arange(nr+1.)/nr)))*(b-a)
    return rr

def matder(nr, z1, z2):
    """
    This function calculates the derivative in Chebyshev space.

    >>> r_icb = 0.5 ; r_cmb = 1.5; n_r_max=65
    >>> d1 = matder(n_r_max, r_icb, r_cmb)
    >>> # Chebyshev grid and data
    >>> rr = chebgrid(n_r_max, r_icb, r_cmb)
    >>> f = sin(rr)
    >>> # Radial derivative
    >>> df = dot(d1, f)

    :param nr: number of radial grid points
    :type nr: int
    :param z1: lower limit of the Gauss-Lobatto grid
    :type z1: float
    :param z2: upper limit of the Gauss-Lobatto grid
    :type z2: float
    :returns: a matrix of dimension (nr,nr) to calculate the derivatives
    :rtype: numpy.ndarray
    """
    nrp = nr+1
    w1 = np.zeros((nrp, nrp), dtype=np.float64)
    zl = z2-z1
    for i in range(nrp):
        for j in range(nrp):
            w1[i, j] = spdel(i, j, nr, zl)

    return w1


def intcheb(f, nr, z1, z2):
    """
    This function integrates an input function f defined on the Gauss-Lobatto grid.

    >>> print(intcheb(f, 65, 0.5, 1.5))

    :param f: an input array
    :type: numpy.ndarray
    :param nr: number of radial grid points
    :type nr: int
    :param z1: lower limit of the Gauss-Lobatto grid
    :type z1: float
    :param z2: upper limit of the Gauss-Lobatto grid
    :type z2: float
    :returns: the integrated quantity
    :rtype: float
    """
    func = lambda i, j: 2.*np.cos(np.pi*i*j/nr)/nr
    w1 = np.fromfunction(func, (nr+1, nr+1))

    w1[:, 0] = w1[:, 0]/2.
    w1[:, nr] = w1[:, nr]/2.
    w1[0, :] = w1[0, :]/2.
    w1[nr, :] = w1[nr, :]/2.

    w2 = np.dot(w1, f)
    int = 0.
    for i in range(0, nr+1, 2):
        int = int-(z2-z1)/(i**2-1)*w2[i]

    return int


def spdel(kr, jr, nr, zl):
    if kr != nr :
        fac = 1.
        k = kr
        j = jr
    else:
        fac = -1.
        k = 0.
        j = nr-jr

    spdel = fac*dnum(k, j, nr)/den(k, j, nr)
    return -spdel*(2./zl)

def dnum(k, j, nr):
    if k == 0:
        if (j == 0 or j == nr):
            dnum = 0.5
            a = nr % 2
            if a == 1:
                dnum = -dnum
            if j == 0:
                dnum = 1./3.*float(nr*nr)+1./6.
            return dnum

        dnum = 0.5*(float(nr)+0.5)*((float(nr)+0.5)+(1./np.tan(np.pi*float(j) \
               /float(2.*nr)))**2)+1./8.-0.25/(np.sin(np.pi*float(j)/ \
               float(2*nr))**2) - 0.5*float(nr*nr)
        return dnum

    dnum = ff(k+j, nr)+ff(k-j, nr)
    return dnum

def ff(i, nr):
    if i == 0:
        return 0
    ff = float(nr)*0.5/np.tan(np.pi*float(i)/float(2.*nr))

    a = i % 2
    if a == 0:
        ff = -ff
    return ff

def den(k, j, nr):
    if k == 0:
        den = 0.5*float(nr)
        a = j % 2
        if a == 1:
            den = -den
        if (j == 0 or j == nr):
            den = 1.
        return den

    den = float(nr)*np.sin(np.pi*float(k)/float(nr))
    if (j == 0 or j == nr):
        den = 2.*den
    return den

def timeder(time,y):
    """
    time derivative of an input array

    computed with central differences (numpy.gradient)

    >>> ts = MagicTs(field='e_kin')
    >>> dt_ekinpol = timeder(ts,field='ekin_pol')

    """
    out = np.gradient(y, time, edge_order=1)

    return out

def secondtimeder(time,y):
    """
    second time derivative of an input array

    computed with central differences (numpy.gradient)

    >>> ts = MagicTs(field='e_kin')
    >>> dt_ekinpol = secondtimeder(ts,field='ekin_pol')

    """
    tmp = np.gradient(y, time, edge_order=1)
    out = np.gradient(tmp, time, edge_order=1)
    return out

def phideravg(data, minc=1, order=4):
    """
    phi-derivative of an input array

    >>> gr = MagicGraph()
    >>> dvphidp = phideravg(gr.vphi, minc=gr.minc)

    :param data: input array
    :type data: numpy.ndarray
    :param minc: azimuthal symmetry
    :type minc: int
    :param order: order of the finite-difference scheme (possible values are 2 or 4)
    :type order: int
    :returns: the phi-derivative of the input array
    :rtype: numpy.ndarray
    """
    nphi = data.shape[0]
    dphi = 2.*np.pi/minc/(nphi-1.)
    if order == 2:
        der = (np.roll(data, -1,  axis=0)-np.roll(data, 1, axis=0))/(2.*dphi)
        der[0, ...] = (data[1, ...]-data[-2, ...])/(2.*dphi)
        der[-1, ...] = der[0, ...]
    elif order == 4:
        der = (   -np.roll(data,-2,axis=0) \
               +8.*np.roll(data,-1,axis=0) \
               -8.*np.roll(data, 1,axis=0)  \
                  +np.roll(data, 2,axis=0)   )/(12.*dphi)
        der[1, ...] = (-data[3, ...]+8.*data[2, ...]-\
                       8.*data[0, ...] +data[-2, ...])/(12.*dphi)
        der[-2, ...] = (-data[0, ...]+8.*data[-1, ...]-\
                       8.*data[-3, ...]+data[-4, ...])/(12.*dphi)
        der[0, ...] = (-data[2, ...]+8.*data[1, ...]-\
                       8.*data[-2, ...] +data[-3, ...])/(12.*dphi)
        der[-1, ...] = der[0, ...]
    return der

def rderavg(data, eta=0.35, spectral=True, exclude=False):
    """
    Radial derivative of an input array

    >>> gr = MagiGraph()
    >>> dvrdr = rderavg(gr.vr, eta=gr.radratio)

    :param data: input array
    :type data: numpy.ndarray
    :param eta: aspect ratio of the spherical shell
    :type eta: float
    :param spectral: when set to True use Chebyshev derivatives, otherwise use
                     finite differences (default is True)
    :type spectral: bool
    :param exclude: when set to True, exclude the first and last radial grid points
                    and replace them by a spline extrapolation (default is False)
    :type exclude: bool
    :returns: the radial derivative of the input array
    :rtype: numpy.ndarray
    """
    r1 = 1./(1.-eta)
    r2 = eta/(1.-eta)
    nr = data.shape[-1]
    grid = chebgrid(nr-1, r1, r2)
    if exclude:
        g = grid[::-1]
        gnew = np.linspace(r2, r1, 1000)
        if len(data.shape) == 2:
            for i in range(data.shape[0]):
                val = data[i, ::-1]
                tckp = S.splrep(g[1:-1], val[1:-1])
                fnew = S.splev(gnew, tckp)
                data[i, 0] = fnew[-1]
                data[i, -1] = fnew[0]
        else:
            for j in range(data.shape[0]):
                for i in range(data.shape[1]):
                    val = data[j, i, ::-1]
                    tckp = S.splrep(g[1:-1], val[1:-1])
                    fnew = S.splev(gnew, tckp)
                    data[j, i, 0] = fnew[-1]
                    data[j, i, -1] = fnew[0]
    if spectral:
        d1 = matder(nr-1, r1, r2)
        if len(data.shape) == 2:
            der = np.tensordot(data, d1, axes=[1, 1])
        else:
            der = np.tensordot(data, d1, axes=[2, 1])
    else:
        denom = np.roll(grid, -1) - np.roll(grid, 1)
        denom[0] = grid[1]-grid[0]
        denom[-1] = grid[-1]-grid[-2]
        der = (np.roll(data, -1,  axis=-1)-np.roll(data, 1, axis=-1))/denom
        der[..., 0] = (data[..., 1]-data[..., 0])/(grid[1]-grid[0])
        der[..., -1] = (data[..., -1]-data[..., -2])/(grid[-1]-grid[-2])
    return der

def thetaderavg(data, order=4):
    """
    Theta-derivative of an input array (finite differences)

    >>> gr = MagiGraph()
    >>> dvtdt = thetaderavg(gr.vtheta)

    :param data: input array
    :type data: numpy.ndarray
    :param order: order of the finite-difference scheme (possible values are 2 or 4)
    :type order: int
    :returns: the theta-derivative of the input array
    :rtype: numpy.ndarray
    """
    if len(data.shape) == 3: # 3-D
        ntheta = data.shape[1]
        dtheta = np.pi/(ntheta-1.)
        if order == 2:
            der = (np.roll(data, -1,  axis=1)-np.roll(data, 1, axis=1))/(2.*dtheta)
            der[:, 0, :] = (data[:, 1, :]-data[:, 0, :])/dtheta
            der[:, -1, :] = (data[:, -1, :]-data[:, -2, :])/dtheta
        elif order == 4:
            der = (   -np.roll(data,-2,axis=1) \
                   +8.*np.roll(data,-1,axis=1) \
                   -8.*np.roll(data, 1,axis=1)  \
                      +np.roll(data, 2,axis=1)   )/(12.*dtheta)
            der[:, 1, :] = (data[:, 2, :]-data[:, 0, :])/(2.*dtheta)
            der[:, -2, :] = (data[:, -1, :]-data[:, -3, :])/(2.*dtheta)
            der[:, 0, :] = (data[:, 1, :]-data[:, 0, :])/dtheta
            der[:, -1, :] = (data[:, -1, :]-data[:, -2, :])/dtheta

    elif len(data.shape) == 2: #2-D
        ntheta = data.shape[0]
        dtheta = np.pi/(ntheta-1.)
        if order == 2:
            der = (np.roll(data, -1,  axis=0)-np.roll(data, 1, axis=0))/(2.*dtheta)
            der[0, :] = (data[1, :]-data[0, :])/dtheta
            der[-1, :] = (data[-1, :]-data[-2, :])/dtheta
        elif order == 4:
            der = (-np.roll(data,-2,axis=0)+8.*np.roll(data,-1,axis=0)-\
                  8.*np.roll(data,1,axis=0)+np.roll(data,2,axis=0))/(12.*dtheta)
            der[1, :] = (data[2, :]-data[0, :])/(2.*dtheta)
            der[-2, :] = (data[-1, :]-data[-3, :])/(2.*dtheta)
            der[0, :] = (data[1, :]-data[0, :])/dtheta
            der[-1, :] = (data[-1, :]-data[-2, :])/dtheta

    return der


def zderavg(data, eta=0.35, spectral=True, colat=None, exclude=False):
    """
    z derivative of an input array

    >>> gr = MagiGraph()
    >>> dvrdz = zderavg(gr.vr, eta=gr.radratio, colat=gr.colatitude)

    :param data: input array
    :type data: numpy.ndarray
    :param eta: aspect ratio of the spherical shell
    :type eta: float
    :param spectral: when set to True use Chebyshev derivatives, otherwise use
                     finite differences (default is True)
    :type spectral: bool
    :param exclude: when set to True, exclude the first and last radial grid points
                    and replace them by a spline extrapolation (default is False)
    :type exclude: bool
    :param colat: colatitudes (when not specified a regular grid is assumed)
    :type colat: numpy.ndarray
    :returns: the z derivative of the input array
    :rtype: numpy.ndarray
    """
    if len(data.shape) == 3:  # 3-D
        ntheta = data.shape[1]
    elif len(data.shape) == 2:  # 2-D
        ntheta = data.shape[0]
    nr = data.shape[-1]
    r1 = 1./(1.-eta)
    r2 = eta/(1.-eta)
    if colat is not None:
        th = colat
    else:
        th = np.linspace(0., np.pi, ntheta)
    rr = chebgrid(nr-1, r1, r2)

    if len(data.shape) == 3:  # 3-D
        thmD = np.zeros_like(data)
        for i in range(ntheta):
            thmD[:,i,:] = th[i]
    elif len(data.shape) == 2:  # 2-D
        thmD = np.zeros((ntheta, nr), np.float64)
        for i in range(ntheta):
            thmD[i, :] = th[i]

    dtheta = thetaderavg(data)
    dr = rderavg(data, eta, spectral, exclude)
    dz = np.cos(thmD)*dr - np.sin(thmD)/rr*dtheta
    return dz

def sderavg(data, eta=0.35, spectral=True, colat=None, exclude=False):
    """
    s derivative of an input array

    >>> gr = MagiGraph()
    >>> dvpds = sderavg(gr.vphi, eta=gr.radratio, colat=gr.colatitude)

    :param data: input array
    :type data: numpy.ndarray
    :param eta: aspect ratio of the spherical shell
    :type eta: float
    :param spectral: when set to True use Chebyshev derivatives, otherwise use
                     finite differences (default is True)
    :type spectral: bool
    :param exclude: when set to True, exclude the first and last radial grid points
                    and replace them by a spline extrapolation (default is False)
    :type exclude: bool
    :param colat: colatitudes (when not specified a regular grid is assumed)
    :type colat: numpy.ndarray
    :returns: the s derivative of the input array
    :rtype: numpy.ndarray
    """
    ntheta = data.shape[0]
    nr = data.shape[-1]
    r1 = 1./(1.-eta)
    r2 = eta/(1.-eta)
    if colat is not None:
        th = colat
    else:
        th = np.linspace(0., np.pi, ntheta)
    rr = chebgrid(nr-1, r1, r2)
    rr2D, th2D = np.meshgrid(rr,th)
    dtheta = thetaderavg(data)
    dr = rderavg(data, eta, spectral, exclude)
    ds = np.sin(th2D)*dr + np.cos(th2D)/rr2D*dtheta
    return ds


def cylSder(radius, data, order=4):
    """
    This function computes the s derivative of an input array defined on
    a regularly-spaced cylindrical grid.

    >>> s = linspace(0., 1., 129 ; dat = cos(s)
    >>> ddatds = cylSder(s, dat)

    :param radius: cylindrical radius
    :type radius: numpy.ndarray
    :param data: input data
    :type data: numpy.ndarray
    :param order: order of the finite-difference scheme (possible values are 2 or 4)
    :type order: int
    :returns: s derivative
    :rtype: numpy.ndarray
    """
    ns = data.shape[-1]
    ds = (radius.max()-radius.min())/(ns-1.)
    if order == 2:
        der = (np.roll(data, -1,  axis=-1)-np.roll(data, 1, axis=-1))/(2.*ds)
        der[..., 0] = (data[..., 1]-data[..., 0])/ds
        der[..., -1] = (data[..., -1]-data[..., -2])/ds
    elif order == 4:
        der = (   -np.roll(data,-2,axis=-1) \
               +8.*np.roll(data,-1,axis=-1) \
               -8.*np.roll(data, 1,axis=-1) \
                  +np.roll(data, 2,axis=-1)   )/(12.*ds)
        der[..., 1] = (data[..., 2]-data[..., 0])/(2.*ds)
        der[..., -2] = (data[..., -1]-data[..., -3])/(2.*ds)
        der[..., 0] = (data[..., 1]-data[..., 0])/ds
        der[..., -1] = (data[..., -1]-data[..., -2])/ds

    return der

def cylZder(z, data):
    """
    This function computes the z derivative of an input array defined on
    a regularly-spaced cylindrical grid.

    >>> z = linspace(-1., 1., 129 ; dat = cos(z)
    >>> ddatdz = cylZder(z, dat)

    :param z: height of the cylinder
    :type z: numpy.ndarray
    :param data: input data
    :type data: numpy.ndarray
    :returns: z derivative
    :rtype: numpy.ndarray
    """
    nz = data.shape[1]
    dz = (z.max()-z.min())/(nz-1.)
    der = (np.roll(data, -1,  axis=1)-np.roll(data, 1, axis=1))/(2.*dz)
    der[:, 0, :] = (data[:, 1, :]-data[:, 0, :])/dz
    der[:, -1, :] = (data[:, -1, :]-data[:, -2, :])/dz
    return der

def getCpuTime(file):
    """
    This function calculates the CPU time from one given log file

    :param file: the log file you want to analyze
    :type file: file
    :returns: the total CPU time
    :rtype: float
    """
    threads_old = re.compile(r'[\s]*\![\s]*nThreads\:[\s]*(.*)')
    threads_new = re.compile(r'[\s]*\![\s]*Number of OMP threads\:[\s]*(.*)')
    ranks = re.compile(r'[\s]*\![\s\w]*ranks[\s\w]*\:[\s]*(.*)')
    runTime = re.compile(r'[\s\!\w]*time:[\s]*([0-9]*)d[\s\:]*([0-9]*)h[\s\:]*([0-9]*)m[\s\:]*([0-9]*)s[\s\:]*([0-9]*)ms.*')
    runTime_new = re.compile(r' \! Total run time:[\s]*([0-9]*)[\s]*h[\s]*([0-9]*)[\s]*m[\s]*([0-9]*)[\s]*s[\s]*([0-9]*)[\s]*ms[\s]*')

    f = open(file, 'r')
    tab = f.readlines()
    nThreads = 1 # In case a pure MPI version is used
    nRanks = 1 # In case the old OpenMP version is used
    realTime = 0.
    for line in tab:
        if threads_old.match(line):
            nThreads = int(threads_old.search(line).groups()[0])
        elif threads_new.match(line):
            nThreads = int(threads_new.search(line).groups()[0])
        elif ranks.match(line):
            nRanks = int(ranks.search(line).groups()[0])
        elif runTime.match(line):
            days = int(runTime.search(line).groups()[0])
            hours = int(runTime.search(line).groups()[1])
            min = int(runTime.search(line).groups()[2])
            sec = int(runTime.search(line).groups()[3])
            ms = int(runTime.search(line).groups()[4])
            realTime = 24*days+hours+1./60*min+1./3600*sec+1./3.6e6*ms
        elif runTime_new.match(line):
            hours = int(runTime_new.search(line).groups()[0])
            min = int(runTime_new.search(line).groups()[1])
            sec = int(runTime_new.search(line).groups()[2])
            ms = int(runTime_new.search(line).groups()[3])
            realTime = hours+1./60*min+1./3600*sec+1./3.6e6*ms
    f.close()
    cpuTime = nThreads*nRanks*realTime

    return cpuTime

def ReadBinaryTimeseries(infile,
                         ncols,
                         datatype='f8',
                         endianness='>'):
    """
    This function reads binary timeseries. It is then faster than
    the fast_read function.

    :param infile: the file to read
    :type infile: string
    :param ncols: number of columns of the file
    :type ncols: int
    :param datatype: 'f8' = 64-bit floating-point number
                     'f4' = 32-bit floating-point number
    :type datatype: string
    :param endianness: '>' = big-endian ; '<' = small-endian
    :type endianness: string
    :returns: an array[nlines, ncols] that contains
              the data of the binary file
    :rtype: numpy.ndarray
    """
    DUMM = endianness+'i4'
    FTYP = endianness+datatype

    size = os.path.getsize(infile)
    #nline = size/(2*4 + ncols*4) # line = 2*i4 + ncols*f4
    typeG = np.dtype([('dum1',DUMM,1),
                      ('line',FTYP,ncols),
                      ('dum2',DUMM,1)])

    with open(infile,'rb') as f:
        data = np.fromfile(f,dtype=typeG,count=size)['line']

    return data

def getTotalRunTime():
    """
    This function calculates the total CPU time of one run directory

    :returns: the total RUN time
    :rtype: float
    """
    logFiles = glob.glob('log.*')
    totCpuTime = 0
    for file in logFiles:
        totCpuTime += getCpuTime(file)

    return totCpuTime

def prime_factors(n):
    """
    This function returns all prime factors of a number

    :type n: int
    :returns: all prime factors
    :rtype: list
    """
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)

    return factors
