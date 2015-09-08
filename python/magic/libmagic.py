# -*- coding: utf-8 -*-
import scipy.interpolate as S
import numpy as N
import glob, os, re, sys
from .npfile import *

__author__  = "$Author$"
__date__   = "$Date$"
__version__ = "$Revision$"


def selectField(obj, field, labTex=True):
    """
    This function selects for you which field you want to display.
    """
    if field in ('Bp', 'bp', 'bphi', 'Bphi'):
        data = obj.Bphi
        if labTex:
            label = r'$B_{\phi}$'
        else:
            label = 'Bphi'
    elif field in ('Bt', 'bt', 'btheta', 'Btheta'):
        data = obj.Btheta
        if labTex:
            label = r'$B_{\theta}$'
        else:
            label = 'Btheta'
    elif field in ('Br', 'br'):
        data = obj.Br
        if labTex:
            label = r'$B_r$'
        else:
            label = 'Br'
    elif field in ('Vr', 'vr', 'Ur', 'ur'):
        data = obj.vr
        if labTex:
            label = r'$v_r$'
        else:
            label = 'vr'
    elif field in ('Vtheta', 'vtheta', 'Utheta', 'utheta', 'vt', 'Vt',
                   'Ut', 'ut'):
        data = obj.vtheta
        if labTex:
            label = r'$v_{\theta}$'
        else:
            label = 'vtheta'
    elif field in ('Vphi', 'vphi', 'Uphi', 'uphi', 'up', 'Up', 'Vp', 'vp'):
        data = obj.vphi
        if labTex:
            label = r'$v_{\phi}$'
        else:
            label = 'vphi'
    elif field in ('entropy', 's', 'S'):
        data = obj.entropy
        label = 'Entropy'
    elif field in ('u2'):
        data = obj.vphi**2+obj.vr**2+obj.vtheta**2
        if labTex:
            label = r'$u^2$'
        else:
            label = 'u2'
    elif field in ('nrj'):
        temp0, rho0, beta = anelprof(obj.radius, obj.strat, obj.polind)
        data = 1./2.*rho0*(obj.vphi**2+obj.vr**2+obj.vtheta**2)
        if labTex:
            label = r'$E_{\hbox{kin}}$'
        else:
            label = 'Ekin'
    elif field in ('b2', 'B2'):
        data = (obj.Bphi**2+obj.Br**2+obj.Btheta**2)
        if labTex:
            label = r'$B^2$'
        else:
            label = 'B2'
    elif field in ('vrconv', 'vrc'):
        data = obj.vr-obj.vr.mean(axis=0)
        if labTex:
            label = r'$v_{r}$ conv'
        else:
            label = 'vr conv'
    elif field in ('vtconv', 'vtc'):
        data = obj.vtheta-obj.vtheta.mean(axis=0)
        if labTex:
            label = r'$v_{\theta}$ conv'
        else:
            label = 'vt conv'
    elif field in ('vpconv', 'vpc'):
        data = obj.vphi-obj.vphi.mean(axis=0)
        if labTex:
            label = r'$v_{\phi}$ conv'
        else:
            label = 'vp conv'
    elif field in ('bpfluct'):
        data = obj.Bphi-obj.Bphi.mean(axis=0)
        if labTex:
            label = r"$B_{\phi}'$"
        else:
            label = "Bp'"
    elif field in ('brfluct'):
        data = obj.Br-obj.Br.mean(axis=0)
        if labTex:
            label = r"$B_r'$"
        else:
            label = "Br'"
    elif field in ('entropyfluct'):
        data = obj.entropy-obj.entropy.mean(axis=0)
        if labTex:
            label = r"$s'$"
        else:
            label = "s'"
    elif field in ('vrea'):
        data = N.zeros_like(obj.vr)
        for i in range(obj.ntheta):
            data[:, i, :] = (obj.vr[:, i, :]-obj.vr[:, -i-1, :])/2.
        if labTex:
            label = r'$v_{r}$ ea'
        else:
            label = 'vr ea'
    elif field in ('vra'):
        data = N.zeros_like(obj.vr)
        for i in range(obj.ntheta):
            data[:, i, :] = (obj.vr[:, i, :]+obj.vr[:, -i-1, :])/2.
        if labTex:
            label = r'$v_{r}$ es'
        else:
            label = 'vr es'
    elif field in ('vpea'):
        data = N.zeros_like(obj.vr)
        for i in range(obj.ntheta):
            data[:, i, :] = (obj.vphi[:, i, :]-obj.vphi[:, -i-1, :])/2.
        if labTex:
            label = r'$v_{\phi}$ ea'
        else:
            label = r'vp ea'
    elif field in ('vpa'):
        data = N.zeros_like(obj.vr)
        for i in range(obj.ntheta):
            data[:, i, :] = (obj.vphi[:, i, :]+obj.vphi[:, -i-1, :])/2.
        if labTex:
            label = r'$v_{\phi}$ es'
        else:
            label = r'vp es'
    elif field in ('tea'):
        data = N.zeros_like(obj.vr)
        for i in range(obj.ntheta):
            data[:, i, :] = (obj.entropy[:, i, :]-obj.entropy[:, -i-1, :])/2.
        if labTex:
            label = r'$s$ ea'
        else:
            label = r's ea'

    return data, label

def avgField(time, field, tstart):
    """
    This subroutine computes the time-average of a time series

    :param time: time
    :param field: the time series of a given field
    :param tstart: the starting time of the averaging
    """
    mask = N.where(abs(time-tstart) == min(abs(time-tstart)), 1, 0)
    ind = N.nonzero(mask)[0][0]
    fac = 1./(time.max()-time[ind])
    avgField = fac*N.trapz(field[ind:], time[ind:])
    return avgField

def writeVpEq(par, tstart):
    """
    This function computes the time-averaged surface zonal flow (and Rolc)

    :param par: a MagicTs object containing the par file
    :param tstart: the starting time of the averaging
    """
    mask = N.where(abs(par.time-tstart) == min(abs(par.time-tstart)), 1, 0)
    ind = N.nonzero(mask)[0][0]
    fac = 1./(par.time.max()-par.time[ind])
    avgReEq = fac*N.trapz(par.reEquat[ind:], par.time[ind:])
    roEq = avgReEq*par.ek*(1.-par.radratio)
    avgRolC = fac*N.trapz(par.rolc[ind:], par.time[ind:])
    st = '%10.3e%5.2f%6.2f%11.3e%11.3e%11.3e' % (par.ek, par.strat, par.pr, 
                                                 par.ra, roEq, avgRolC)
    return st

def progressbar(it, prefix = "", size = 60):
    count = len(it)
    def _show(_i):
        x = int(size*_i/count)
        sys.stdout.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), _i, count))
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
    :param tfix: in case you want to add only the files that are more recent than
                 a certain date, use tfix
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

def hammer2cart(ttheta, pphi, colat=False):
    if not colat: # for lat and phi \in [-pi, pi]
        xx = 2.*N.sqrt(2.) * N.cos(ttheta)*N.sin(pphi/2.)\
             /N.sqrt(1.+N.cos(ttheta)*N.cos(pphi/2.))
        yy = N.sqrt(2.) * N.sin(ttheta)\
             /N.sqrt(1.+N.cos(ttheta)*N.cos(pphi/2.))
    else:  # for colat and phi \in [0, 2pi]
        xx = -2.*N.sqrt(2.) * N.sin(ttheta)*N.cos(pphi/2.)\
             /N.sqrt(1.+N.sin(ttheta)*N.sin(pphi/2.))
        yy = N.sqrt(2.) * N.cos(ttheta)\
             /N.sqrt(1.+N.sin(ttheta)*N.sin(pphi/2.))
    return xx, yy

def cut(dat, vmax=None, vmin=None):
    """
    Replace data by vmax if data > vmax
    or by vmin if data < vmin
    """
    if vmax is not None:
        mask = N.where(dat>=vmax, 1, 0)
        dat = dat*(mask == 0) + vmax*(mask == 1)
        normed = False
    if vmin is not None:
        mask = N.where(dat<=vmin, 1, 0)
        dat = dat*(mask == 0) + vmin*(mask == 1)
        normed = False
    return dat

def symmetrize(data, ms):
    """
    Symmetrise an array which is defined only on minc=ms

    :param data: the input array
    :param ms: the initial symmetry
    """
    np = data.shape[0]*ms +1
    size = [np]
    size.append(data.shape[1])
    if len(data.shape) == 3:
        size.append(data.shape[2])
    out = N.zeros(size, dtype=data.dtype)
    for i in range(ms):
        out[i*data.shape[0]:(i+1)*data.shape[0], ...] = data
    out[-1, ...] = out[0, ...]
    return out

def fast_read(file, skiplines=0, binary=False, precision='Float64'):
    """
    This function reads an input ascii table 
    (can read both formatted or unformatted fortran)

    :param file: name of the input file
    :param skiplines: number of header lines to be skept
    :param binary: is it a formatted or unformatted file?
    :param precision: single ('Float32') or double precision ('Float64')
    """
    if not binary:
        f = open(file, 'r')
        X = []
        for k, line in enumerate(f.readlines()):
            st = line.replace('D', 'E')
            if k >= skiplines:
                X.append(st.split())
        X = N.array(X, dtype=precision)
        f.close()
    else:
        f = npfile(file, endian='B')
        X = []
        while 1:
            try:
                X.append(f.fort_read(precision))
            except TypeError:
                break
        X = N.array(X, dtype=precision)
        f.close()
    return X


def varmax(tag):
    return len(glob.glob('G_[0-9]*.%s' % tag))


def anelprof(radius, strat, polind, g0=0., g1=0., g2=1.):
    """
    This functions calculates the reference temperature and density profiles
    of an anelastic model.

    :param radius: the radial gridpoints
    :param polind: the polytropic index
    :param strat: the number of the density scale heights between the inner
                  and the outer boundary
    :param g0: gravity profile: g=g0
    :param g1: gravity profile: g=g1*r/r_o
    :param g2: gravity profile: g=g2*(r_o/r)**2
    """
    if radius[-1] < radius[0]:
        ro = radius[0]
        radratio = radius[-1]/radius[0]
    else:
        ro = radius[-1]
        radratio = radius[0]/radius[-1]
    grav = g0 + g1 * radius/ro + g2 * (ro/radius)**2
    ofr=( N.exp(strat/polind)-1. )/ ( g0+0.5*g1*(1.+radratio) +g2/radratio )

    temp0 = -ofr*( g0*radius +0.5*g1*radius**2/ro-g2*ro**2/radius ) + \
            1.+ ofr*ro*(g0+0.5*g1-g2)
    rho0 = temp0**polind
    beta = -ofr*grav*polind/temp0
    return temp0, rho0, beta

def chebgrid(nr, a, b):
    """
    This function defines a Gauss-Lobatto grid from a to b.

    :param nr: number of radial grid points
    :param a: lower limit of the Gauss-Lobatto grid
    :param b: upper limit of the Gauss-Lobatto grid
    """
    rst = (a+b)/(b-a)
    rr = 0.5*(rst+N.cos(N.pi*(1.-N.arange(nr+1.)/nr)))*(b-a)
    return rr

def matder(nr, z1, z2):
    """ 
    This function calculates the derivative in Chebyshev space.

    >>> r_icb = 0.5 ; r_cmb = 1.5
    >>> d1 = matder(65, r_icb, r_cmb)

    :param nr: number of radial grid points
    :param z1: lower limit of the Gauss-Lobatto grid
    :param z2: upper limit of the Gauss-Lobatto grid
    """
    nrp = nr+1
    w1 = N.zeros((nrp, nrp), dtype='Float64')
    zl = z2-z1
    for i in range(nrp):
        for j in range(nrp):
            w1[i, j] = spdel(i, j, nr, zl)
     
    return w1


def intcheb(f, nr, z1, z2):
    """ 
    This function integrates an input function f defined on the Gauss-Lobatto grid.

    :param f: an input array
    :param nr: number of radial grid points
    :param z1: lower limit of the Gauss-Lobatto grid
    :param z2: upper limit of the Gauss-Lobatto grid
    """
    func = lambda i, j: 2.*N.cos(N.pi*i*j/nr)/nr
    w1 = N.fromfunction(func, (nr+1, nr+1))
     
    w1[:, 0] = w1[:, 0]/2.
    w1[:, nr] = w1[:, nr]/2.
    w1[0, :] = w1[0, :]/2.
    w1[nr, :] = w1[nr, :]/2.
     
    w2 = N.dot(w1, f)
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
     
        dnum = 0.5*(float(nr)+0.5)*((float(nr)+0.5)+(1./N.tan(N.pi*float(j) \
               /float(2.*nr)))**2)+1./8.-0.25/(N.sin(N.pi*float(j)/ \
               float(2*nr))**2) - 0.5*float(nr*nr)
        return dnum
     
    dnum = ff(k+j, nr)+ff(k-j, nr)
    return dnum

def ff(i, nr):
    if i == 0:
        return 0
    ff = float(nr)*0.5/N.tan(N.pi*float(i)/float(2.*nr))
     
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
     
    den = float(nr)*N.sin(N.pi*float(k)/float(nr))
    if (j == 0 or j == nr):
        den = 2.*den
    return den


def phideravg(data, minc=1):
    """
    phi-derivative of an array

    :param data: input array
    :param minc: azimuthal symmetry
    """
    nphi = data.shape[0]
    dphi = 2.*N.pi/minc/(nphi-1.)
    der = (N.roll(data, -1,  axis=0)-N.roll(data, 1, axis=0))/(2.*dphi)
    der[0, ...] = (data[1, ...]-data[-2, ...])/(2.*dphi)
    der[-1, ...] = der[0, ...]
    return der

def rderavg(data, eta=0.35, spectral=True, exclude=False):
    r1 = 1./(1.-eta)
    r2 = eta/(1.-eta)
    nr = data.shape[-1]
    grid = chebgrid(nr-1, r1, r2)
    if exclude:
        g = grid[::-1]
        gnew = N.linspace(r2, r1, 1000)
        for i in range(data.shape[0]):
            val = data[i, ::-1]
            tckp = S.splrep(g[1:-1], val[1:-1])
            fnew = S.splev(gnew, tckp)
            data[i, 0] = fnew[-1]
            data[i, -1] = fnew[0]
    if spectral:
        d1 = matder(nr-1, r1, r2)
        if len(data.shape) == 2:
            der = N.tensordot(data, d1, axes=[1, 1])
        else:
            der = N.tensordot(data, d1, axes=[2, 1])
    else:
        denom = N.roll(grid, -1) - N.roll(grid, 1)
        der = (N.roll(data, -1,  axis=-1)-N.roll(data, 1, axis=-1))/denom
        der[:, 0] = (data[:, 1]-data[:, 0])/(grid[1]-grid[0])
        der[:, -1] = (data[:, -1]-data[:, -2])/(grid[-1]-grid[-2])
    return der

def thetaderavg(data, order=4):
    if len(data.shape) == 3: # 3-D
        ntheta = data.shape[1]
        dtheta = N.pi/(ntheta-1.)
        if order == 2:
            der = (N.roll(data, -1,  axis=1)-N.roll(data, 1, axis=1))/(2.*dtheta)
            der[:, 0, :] = (data[:, 1, :]-data[:, 0, :])/dtheta
            der[:, -1, :] = (data[:, -1, :]-data[:, -2, :])/dtheta
        elif order == 4:
            der = (   -N.roll(data,-2,axis=1) \
                   +8.*N.roll(data,-1,axis=1) \
                   -8.*N.roll(data, 1,axis=1)  \
                      +N.roll(data, 2,axis=1)   )/(12.*dtheta)
            der[:, 1, :] = (data[:, 2, :]-data[:, 0, :])/(2.*dtheta)
            der[:, -2, :] = (data[:, -1, :]-data[:, -3, :])/(2.*dtheta)
            der[:, 0, :] = (data[:, 1, :]-data[:, 0, :])/dtheta
            der[:, -1, :] = (data[:, -1, :]-data[:, -2, :])/dtheta

    elif len(data.shape) == 2: #2-D
        ntheta = data.shape[0]
        dtheta = N.pi/(ntheta-1.)
        if order == 2:
            der = (N.roll(data, -1,  axis=0)-N.roll(data, 1, axis=0))/(2.*dtheta)
            der[0, :] = (data[1, :]-data[0, :])/dtheta
            der[-1, :] = (data[-1, :]-data[-2, :])/dtheta
        elif order == 4:
            der = (-N.roll(data,-2,axis=0)+8.*N.roll(data,-1,axis=0)-\
                  8.*N.roll(data,1,axis=0)+N.roll(data,2,axis=0))/(12.*dtheta)
            der[1, :] = (data[2, :]-data[0, :])/(2.*dtheta)
            der[-2, :] = (data[-1, :]-data[-3, :])/(2.*dtheta)
            der[0, :] = (data[1, :]-data[0, :])/dtheta
            der[-1, :] = (data[-1, :]-data[-2, :])/dtheta

    return der


def zderavg(data, eta=0.35, spectral=True, colat=None, exclude=False):
    ntheta = data.shape[0]
    nr = data.shape[-1]
    r1 = 1./(1.-eta)
    r2 = eta/(1.-eta)
    if colat is not None:
        th = colat
    else:
        th = N.linspace(0., N.pi, ntheta)
    rr = chebgrid(nr-1, r1, r2)
    rr2D, th2D = N.meshgrid(rr,th)
    dtheta = thetaderavg(data)
    dr = rderavg(data, eta, spectral, exclude)
    dz = N.cos(th2D)*dr - N.sin(th2D)/rr2D*dtheta
    return dz

def sderavg(data, eta=0.35, spectral=True, colat=None, exclude=False):
    ntheta = data.shape[0]
    nr = data.shape[-1]
    r1 = 1./(1.-eta)
    r2 = eta/(1.-eta)
    if colat is not None:
        th = colat
    else:
        th = N.linspace(0., N.pi, ntheta)
    rr = chebgrid(nr-1, r1, r2)
    rr2D, th2D = N.meshgrid(rr,th)
    dtheta = thetaderavg(data)
    dr = rderavg(data, eta, spectral, exclude)
    ds = N.sin(th2D)*dr + N.cos(th2D)/rr2D*dtheta
    return ds

def sder3D(data, eta=0.35, spectral=True, colat=None, exclude=False):
    ntheta = data.shape[1]
    nr = data.shape[-1]
    r1 = 1./(1.-eta)
    r2 = eta/(1.-eta)
    if colat is not None:
        th = colat
    else:
        th = N.linspace(0., N.pi, ntheta)
    rr = chebgrid(nr-1, r1, r2)
    th3D = N.zeros_like(data)
    rr3D = N.zeros_like(data)
    for i in range(ntheta):
        th3D[:, i, :] = th[i]
    for i in range(nr):
        rr3D[:, :, i] = rr[i]
    dtheta = thetaderavg(data)
    dr = rderavg(data, eta, spectral, exclude)
    ds = N.sin(th3D)*dr + N.cos(th3D)/rr3D*dtheta
    return ds

def cylSder(radius, data):
    ns = data.shape[-1]
    ds = radius.max()/(ns-1.)
    der = (N.roll(data, -1,  axis=-1)-N.roll(data, 1, axis=-1))/(2.*ds)
    der[..., 0] = (data[..., 1]-data[..., 0])/ds
    der[..., -1] = (data[..., -1]-data[..., -2])/ds
    return der

def cylZder(z, data):
    nz = data.shape[1]
    dz = (z.max()-z.min())/(nz-1.)
    der = (N.roll(data, -1,  axis=1)-N.roll(data, 1, axis=1))/(2.*dz)
    der[:, 0, :] = (data[:, 1, :]-data[:, 0, :])/dz
    der[:, -1, :] = (data[:, -1, :]-data[:, -2, :])/dz
    return der

def getCpuTime(file):
    """
    This function calculates the CPU time from one given log file

    :param file: the log file you want to analyze
    """
    threads = re.compile(r'[\s]*\![\s]*nThreads\:[\s]*(.*)')
    ranks = re.compile(r'[\s]*\![\s\w]*ranks[\s\w]*\:[\s]*(.*)')
    runTime = re.compile(r'[\s\!\w]*time:[\s]*([0-9]*)d[\s\:]*([0-9]*)h[\s\:]*([0-9]*)m[\s\:]*([0-9]*)s[\s\:]*([0-9]*)ms.*')
    f = open(file, 'r')
    tab = f.readlines()
    nThreads = 1 # In case a pure MPI version is used
    nRanks = 1 # In case the old OpenMP version is used
    realTime = 0.
    for line in tab:
        if threads.match(line):
            nThreads = int(threads.search(line).groups()[0])
        elif ranks.match(line):
            nRanks = int(ranks.search(line).groups()[0])
        elif runTime.match(line):
            days = int(runTime.search(line).groups()[0])
            hours = int(runTime.search(line).groups()[1])
            min = int(runTime.search(line).groups()[2])
            sec = int(runTime.search(line).groups()[3])
            ms = int(runTime.search(line).groups()[4])
            realTime = 24*days+hours+1./60*min+1./3600*sec+1./3.6e6*ms
    f.close()
    cpuTime = nThreads*nRanks*realTime
    return cpuTime

def getTotalRunTime():
    """
    This function calculates the toal CPU time of one run directory
    """
    logFiles = glob.glob('log.*')
    totCpuTime = 0
    for file in logFiles:
        totCpuTime += getCpuTime(file)
    return totCpuTime
