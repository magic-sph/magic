# -*- coding: utf-8 -*-
import os
import sys
import matplotlib.pyplot as P
try:
    import configparser as CoPa
except ImportError:
    import ConfigParser as CoPa

__author__  = "$Author$"
__date__   = "$Date$"
__version__ = "$Revision$"


pythonVersion = sys.version_info.major

if 'MAGIC_HOME' in os.environ:
    path = os.environ['MAGIC_HOME']
    magicdir = path+'/python/magic'

    parser = CoPa.RawConfigParser()
    parser.read(magicdir + '/magic.cfg')
    backend = parser.get('plots', 'backend')
    labTex = parser.getboolean('plots', 'labTex')
    # Do you want to build so libraries -> need f2py + ifort or gfortran (> 4.1)
    buildSo = parser.getboolean('libraries', 'buildLib')
    compiler = parser.get('libraries', 'compiler')
    f2pycmd = parser.get('libraries', 'f2pyexec')
else: # Default if the PATH is messed up
    backend = 'GTKAgg'
    labTex = False
    buildSo = False

#
# Plots setup
#

P.switch_backend(backend)
P.rc('xtick.major', size=7)
P.rc('xtick.minor', size=3.5)
P.rc('ytick.major', size=7)
P.rc('ytick.minor', size=3.5)
P.rc('axes.formatter', limits=(-5,5))
#P.rc('axes', color_cycle=('30a2da', 'fc4f30', 'e5ae38', '6d904f', '8b8b8b'))

if labTex:
    P.rc('figure.subplot', right=0.95, top=0.95, hspace=0.24)
    P.rc('xtick', labelsize=14)
    P.rc('ytick', labelsize=14)
    P.rc('legend', fontsize=14)
    P.rc('axes', labelsize=20)
    P.rc('text', usetex=True)
    P.rc('font', size=18, **{'family':'serif','serif':['Computer Modern']})
else:
    P.rc('figure.subplot', right=0.95, top=0.95, hspace=0.24)
    P.rc('xtick', labelsize=12)
    P.rc('ytick', labelsize=12)
    P.rc('legend', fontsize=12)
    P.rc('axes', labelsize=16)
    P.rc('font', size=14)

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
        cmd = '%s --fcompiler=%s -c -m --opt=-O3 greader_single%i readG_single.f90 &> /dev/null' % (f2pycmd, compiler, pythonVersion)
        os.system(cmd)
        if pythonVersion == 3:
            cmd = 'mv greader_single3.cpython-33m.so %s/greader_single3.so' % magicdir
        elif pythonVersion == 2:
            cmd = 'mv greader_single2.so %s' % magicdir
        os.system(cmd)
        os.chdir(magicdir)

    if not os.path.exists('greader_double%i.so' % pythonVersion):
        os.chdir('fortranLib')
        print("Please wait: building greader_double...")
        cmd = '%s --fcompiler=%s -c -m --opt=-O3 greader_double%i readG_double.f90 &> /dev/null' % (f2pycmd, compiler, pythonVersion)
        os.system(cmd)
        if pythonVersion == 3:
            cmd = 'mv greader_double3.cpython-33m.so %s/greader_double3.so' % magicdir
        elif pythonVersion == 2:
            cmd = 'mv greader_double2.so %s' % magicdir
        os.system(cmd)
        os.chdir(magicdir)

    # For the potential field extrapolation
    if not os.path.exists('potential%i.so' % pythonVersion):
        os.chdir('fortranLib')
        print("Please wait: building potential extrapolation...")
        cmd = '%s --fcompiler=%s -c -m --opt=-O3 potential%i spec.f90 &> /dev/null' % (f2pycmd, compiler, pythonVersion)
        os.system(cmd)
        if pythonVersion == 3:
            cmd = 'mv potential3.cpython-33m.so %s/potential3.so' % magicdir
        elif pythonVersion == 2:
            cmd = 'mv potential2.so %s' % magicdir
        os.system(cmd)
        os.chdir(magicdir)

    # For the vtk file format convertion
    if not os.path.exists('vtklib%i.so' % pythonVersion):
        os.chdir('fortranLib')
        print("Please wait: building vtklib...")
        cmd = '%s --fcompiler=%s -c -m --opt=-O3 vtklib%i vtkLib.f90 &> /dev/null' % (f2pycmd, compiler, pythonVersion)
        os.system(cmd)
        if pythonVersion == 3:
            cmd = 'mv vtklib3.cpython-33m.so %s/vtklib3.so' % magicdir
        elif pythonVersion == 2:
            cmd = 'mv vtklib2.so %s' % magicdir
        os.system(cmd)
        os.chdir(magicdir)

    os.chdir(startdir)

