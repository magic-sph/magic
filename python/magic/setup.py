# -*- coding: utf-8 -*-
import os
import pylab as P
import ConfigParser as CoPa

__author__  = "$Author$"
__date__   = "$Date$"
__version__ = "$Revision$"


if os.environ.has_key('MAGIC_HOME'):
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
    if not os.path.exists('greader.so'):
        os.chdir('fortranLib')
        print "Please wait: building greader..."
        cmd = '%s --fcompiler=%s -c -m --opt=-O3 greader readG.f90 &> /dev/null' % (f2pycmd, compiler)
        os.system(cmd)
        cmd = 'mv greader.so %s' % magicdir
        os.system(cmd)
        os.chdir(magicdir)

    # For the potential field extrapolation
    if not os.path.exists('potential.so'):
        os.chdir('fortranLib')
        print "Please wait: building potential extrapolation..."
        cmd = '%s --fcompiler=%s -c -m --opt=-O3 potential spec.f90 &> /dev/null' % (f2pycmd, compiler)
        os.system(cmd)
        cmd = 'mv potential.so %s' % magicdir
        os.system(cmd)
        os.chdir(magicdir)

    # For the vtk file format convertion
    if not os.path.exists('vtklib.so'):
        os.chdir('fortranLib')
        print "Please wait: building vtklib..."
        cmd = '%s --fcompiler=%s -c -m --opt=-O3 vtklib vtkLib.f90 &> /dev/null' % (f2pycmd, compiler)
        os.system(cmd)
        cmd = 'mv vtklib.so %s' % magicdir
        os.system(cmd)
        os.chdir(magicdir)

    os.chdir(startdir)

