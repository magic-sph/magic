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
plt.rc('figure.subplot', right=0.97, top=0.97, hspace=0.24)

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

    # For reading G files
    if not os.path.exists('greader_single%i.so' % pythonVersion):
        os.chdir('fortranLib')
        print("Please wait: building greader_single...")
        sp.call(['%s' % f2pycmd,
                 '--fcompiler=%s' % fcompiler,
                 '--compiler=%s' % ccompiler,
                 '-c', '-m', '--opt=-O3', 
                 'greader_single%i' % pythonVersion,
                 'readG_single.f90'],  stderr=sp.PIPE, stdout=sp.PIPE)
        if pythonVersion == 3:
            cmd = "mv greader_single3.cpython-%sm* %s/greader_single3.so" % \
                  (pythonSuffix, magicdir)
            sp.call(cmd, shell=True)
        elif pythonVersion == 2:
            sp.call(['mv', 'greader_single2.so', '%s' % magicdir])
        os.chdir(magicdir)

    if not os.path.exists('greader_double%i.so' % pythonVersion):
        os.chdir('fortranLib')
        print("Please wait: building greader_double...")
        sp.call(['%s' % f2pycmd,
                 '--fcompiler=%s' % fcompiler,
                 '--compiler=%s' % ccompiler,
                 '-c', '-m', '--opt=-O3', 
                 'greader_double%i' % pythonVersion,
                 'readG_double.f90'],  stderr=sp.PIPE, stdout=sp.PIPE)
        if pythonVersion == 3:
            cmd = "mv greader_double3.cpython-%sm* %s/greader_double3.so" % \
                  (pythonSuffix, magicdir)
            sp.call(cmd, shell=True)
        elif pythonVersion == 2:
            sp.call(['mv', 'greader_double2.so', '%s' % magicdir])
        os.chdir(magicdir)

    # For the potential field extrapolation
    if not os.path.exists('potential%i.so' % pythonVersion):
        os.chdir('fortranLib')
        print("Please wait: building potential extrapolation...")
        sp.call(['%s' % f2pycmd,
                 '--fcompiler=%s' % fcompiler,
                 '--compiler=%s' % ccompiler,
                 '-c', '-m', '--opt=-O3', 
                 'potential%i' % pythonVersion,
                 'spec.f90'],  stderr=sp.PIPE, stdout=sp.PIPE)
        if pythonVersion == 3:
            cmd = "mv potential3.cpython-%sm* %s/potential3.so" % \
                  (pythonSuffix, magicdir)
            sp.call(cmd, shell=True)
        elif pythonVersion == 2:
            sp.call(['mv', 'potential2.so', '%s' % magicdir])
        os.chdir(magicdir)

    # For the vtk file format convertion
    if not os.path.exists('vtklib%i.so' % pythonVersion):
        os.chdir('fortranLib')
        print("Please wait: building vtklib...")
        sp.call(['%s' % f2pycmd,
                 '--fcompiler=%s' % fcompiler,
                 '--compiler=%s' % ccompiler,
                 '-c', '-m', '--opt=-O3', 
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

