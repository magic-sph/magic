# -*- coding: utf-8 -*-
import subprocess as sp
import os
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from .libmagic import scanDir
try:
    import configparser as CoPa
except ImportError:
    import ConfigParser as CoPa

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
#plt.switch_backend(backend)
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
        omp_options = "--f90flags='-qopenmp'"
        omp_link = '-liomp5'
    elif fcompiler == 'gnu' or fcompiler == 'g95' or fcompiler == 'gnu95':
        f90options = '-O3 -march=native'
        omp_options = "--f90flags='-fopenmp'"
        omp_link = '-lgomp'
    else:
        f90options = '-O3'
        omp_options = ''
        omp_link = ''

    def buildLib(fileName,libName):
        t2 = os.stat('fortranLib/' + fileName).st_mtime
        sos = scanDir(libName + '.*')
        if len(sos) >= 1:
            t1 = os.stat(sos[-1]).st_mtime
        else: # in case the file does not exist t2 is set to t1
            t1 = t2
        if len(sos) < 1  or t2 > t1:

            if (sys.version_info.major == 3 and sys.version_info.minor < 12):
                print("Please wait: building %s using distutils..." %libName)
                return_code = sp.call(['{}'.format(f2pycmd),
                        '--fcompiler={}'.format(fcompiler),
                        '--opt={}'.format(f90options),
                        '-c', '-m',
                        '%s' %libName,
                        'fortranLib/%s' %fileName],  stderr=sp.PIPE, stdout=sp.PIPE)
            else:
                print("Please wait: building %s using meson..." %libName)
                my_env = os.environ.copy()
                my_env["FFLAGS"]=f90options
                return_code = sp.call(['{}'.format(f2pycmd),
                                       '-c', '-m',
                                       '%s' %libName,
                                       'fortranLib/%s' %fileName,
                                       "--backend","meson"],
                                       env=my_env,
                                       stderr=sp.PIPE, stdout=sp.PIPE)

            return return_code
        else:
            return 0

    fortranFiles=['readG_single.f90','readG_double.f90',
                  'readPot_single.f90','legendre.f90',
                  'vtkLib.f90','cyl.f90']

    sharedLibFiles=['greader_single','greader_double',
                    'lmrreader_single','legendre',
                    'vtklib','cylavg']

    for fileName,libName in zip(fortranFiles,sharedLibFiles):
        return_code = buildLib(fileName,libName)
        if return_code == 0:
            print("%s built successfully" %libName)
        else:
            print("Error in building %s" %libName)

    os.chdir(startdir)
