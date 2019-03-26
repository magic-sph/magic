# -*- coding: utf-8 -*-
import subprocess as sp
import os
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
try:
    import configparser as CoPa
except ImportError:
    import ConfigParser as CoPa

pythonVersion = sys.version_info.major

if pythonVersion == 3:
    pythonSuffix = '%i%i' % (pythonVersion, sys.version_info.minor)

if 'MAGIC_HOME' in os.environ:
    path = os.environ['MAGIC_HOME']
    magicdir = path+'/python/magic'

    parser = CoPa.RawConfigParser()
    parser.read(magicdir + '/magic.cfg')
    backend = parser.get('plots', 'backend')
    labTex = parser.getboolean('plots', 'labTex')
    try:
        defaultCm = parser.get('plots', 'defaultCm')
        defaultLevels = parser.getint('plots', 'defaultLevels')
    except CoPa.NoOptionError:
        defaultCm = 'seismic'
        defaultLevels = 65
    # Do you want to build so libraries -> need f2py + ifort or gfortran (> 4.1)
    buildSo = parser.getboolean('libraries', 'buildLib')
    fcompiler = parser.get('libraries', 'fcompiler')
    ccompiler = parser.get('libraries', 'ccompiler')
    f2pycmd = parser.get('libraries', 'f2pyexec')
else: # Default if the PATH is messed up
    backend = 'GTKAgg'
    labTex = False
    defaultCm = 'seismic'
    defaultLevels = 65
    buildSo = False

#
# Plots setup
#
plt.switch_backend(backend)
if int(mpl.__version__[0]) < 2:
    plt.rc('xtick.major', size=7, width=1)
    plt.rc('xtick.minor', size=3.5, width=1)
    plt.rc('ytick.major', size=7, width=1)
    plt.rc('ytick.minor', size=3.5, width=1)
    plt.rc('axes.formatter', limits=(-5,5))
    if mpl.__version__ >= '1.5':
        from cycler import cycler
        colors = ['#30a2da', '#6d904f', '#fc4f30', '#e5ae38', '#7a68a6','#ffb5b8',
                  '#8b8b8b', '#988ed5']
        plt.rc('axes', prop_cycle=(cycler('color', colors)))
    else:
        plt.rc('axes', color_cycle=('30a2da', '6d904f', 'fc4f30', 'e5ae38', '7a68a6',
                                    'ffb5b8', '8b8b8b', '988ed5'))
    plt.rc('lines', linewidth=1.5)
    plt.rc('figure.subplot', right=0.97, top=0.96, hspace=0.24)

    if labTex:
        plt.rc('xtick', labelsize=14)
        plt.rc('ytick', labelsize=14)
        plt.rc('legend', fontsize=14)
        plt.rc('axes', labelsize=20)
        plt.rc('text', usetex=True)
        plt.rc('font', size=18, **{'family':'serif','serif':['Computer Modern']})
    else:
        plt.rc('xtick', labelsize=12)
        plt.rc('ytick', labelsize=12)
        plt.rc('legend', fontsize=12)
        plt.rc('axes', labelsize=16)
        plt.rc('font', size=14)

#
# Library setup
#

if buildSo:
    startdir = os.getcwd()
    os.chdir(magicdir)

    if fcompiler.startswith('intel'):
        f90options = '-O3 -xHost'
    elif fcompiler == 'gnu' or fcompiler == 'g95':
        f90options = '-O3 -march=native'
    else:
        f90options = '-O3'

    # For reading G files
    t1 = os.stat('greader_single%i.so' % pythonVersion).st_mtime
    t2 = os.stat('fortranLib/readG_single.f90').st_mtime
    if not os.path.exists('greader_single%i.so' % pythonVersion) or t2 > t1:
        os.chdir('fortranLib')
        print("Please wait: building greader_single...")
        sp.call(['%s' % f2pycmd,
                 '--fcompiler=%s' % fcompiler,
                 '--compiler=%s' % ccompiler,
                 '--opt=%s' % f90options,
                 '-c', '-m',
                 'greader_single%i' % pythonVersion,
                 'readG_single.f90'],  stderr=sp.PIPE, stdout=sp.PIPE)
        if pythonVersion == 3:
            cmd = "mv greader_single3.cpython-%sm* %s/greader_single3.so" % \
                  (pythonSuffix, magicdir)
            sp.call(cmd, shell=True)
        elif pythonVersion == 2:
            sp.call(['mv', 'greader_single2.so', '%s' % magicdir])
        os.chdir(magicdir)

    t1 = os.stat('greader_double%i.so' % pythonVersion).st_mtime
    t2 = os.stat('fortranLib/readG_double.f90').st_mtime
    if not os.path.exists('greader_double%i.so' % pythonVersion) or t2 > t1:
        os.chdir('fortranLib')
        print("Please wait: building greader_double...")
        sp.call(['%s' % f2pycmd,
                 '--fcompiler=%s' % fcompiler,
                 '--compiler=%s' % ccompiler,
                 '--opt=%s' % f90options,
                 '-c', '-m',
                 'greader_double%i' % pythonVersion,
                 'readG_double.f90'],  stderr=sp.PIPE, stdout=sp.PIPE)
        if pythonVersion == 3:
            cmd = "mv greader_double3.cpython-%sm* %s/greader_double3.so" % \
                  (pythonSuffix, magicdir)
            sp.call(cmd, shell=True)
        elif pythonVersion == 2:
            sp.call(['mv', 'greader_double2.so', '%s' % magicdir])
        os.chdir(magicdir)

    # For the reader of the potential files
    t1 = os.stat('lmrreader_single%i.so' % pythonVersion).st_mtime
    t2 = os.stat('fortranLib/readPot_single.f90').st_mtime
    if not os.path.exists('lmrreader_single%i.so' % pythonVersion) or t2 > t1:
        os.chdir('fortranLib')
        print("Please wait: building lmrreader_single...")
        sp.call(['%s' % f2pycmd,
                 '--fcompiler=%s' % fcompiler,
                 '--compiler=%s' % ccompiler,
                 '--opt=%s' % f90options,
                 '-c', '-m',
                 'lmrreader_single%i' % pythonVersion,
                 'readPot_single.f90'],  stderr=sp.PIPE, stdout=sp.PIPE)
        if pythonVersion == 3:
            cmd = "mv lmrreader_single3.cpython-%sm* %s/lmrreader_single3.so" % \
                  (pythonSuffix, magicdir)
            sp.call(cmd, shell=True)
        elif pythonVersion == 2:
            sp.call(['mv', 'lmrreader_single2.so', '%s' % magicdir])
        os.chdir(magicdir)

    # For the Legendre transforms
    t1 = os.stat('legendre%i.so' % pythonVersion).st_mtime
    t2 = os.stat('fortranLib/legendre.f90').st_mtime
    if not os.path.exists('legendre%i.so' % pythonVersion) or t2 > t1:
        os.chdir('fortranLib')
        print("Please wait: building Legendre transforms...")
        sp.call(['%s' % f2pycmd,
                 '--fcompiler=%s' % fcompiler,
                 '--compiler=%s' % ccompiler,
                 '--opt=%s' % f90options,
                 '-c', '-m',
                 'legendre%i' % pythonVersion,
                 'legendre.f90'],  stderr=sp.PIPE, stdout=sp.PIPE)
        if pythonVersion == 3:
            cmd = "mv legendre3.cpython-%sm* %s/legendre3.so" % \
                  (pythonSuffix, magicdir)
            sp.call(cmd, shell=True)
        elif pythonVersion == 2:
            sp.call(['mv', 'legendre2.so', '%s' % magicdir])
        os.chdir(magicdir)

    # For the vtk file format convertion
    t1 = os.stat('vtklib%i.so' % pythonVersion).st_mtime
    t2 = os.stat('fortranLib/vtkLib.f90').st_mtime
    if not os.path.exists('vtklib%i.so' % pythonVersion) or t2 > t1:
        os.chdir('fortranLib')
        print("Please wait: building vtklib...")
        sp.call(['%s' % f2pycmd,
                 '--fcompiler=%s' % fcompiler,
                 '--compiler=%s' % ccompiler,
                 '--opt=%s' % f90options,
                 '-c', '-m',
                 'vtklib%i' % pythonVersion,
                 'vtkLib.f90'],  stderr=sp.PIPE, stdout=sp.PIPE)
        if pythonVersion == 3:
            cmd = "mv vtklib3.cpython-%sm* %s/vtklib3.so" % \
                  (pythonSuffix, magicdir)
            sp.call(cmd, shell=True)
        elif pythonVersion == 2:
            sp.call(['mv', 'vtklib2.so', '%s' % magicdir])
        os.chdir(magicdir)

    os.chdir(startdir)
