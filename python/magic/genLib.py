#!/usr/bin/python2
import os,sys
"""
This routine is used to generate the library files with f2py
"""

# f2py2 -c --help-fcompiler
#compiler = 'gnu95'
compiler = 'intelem'

f2pycmd = ''
path = os.environ['PATH'].split(':')
# This can be f2py of f2py2 depending on your version
for dir in path:
    if os.path.exists(dir+'/f2py2'):
        f2pycmd = 'f2py2'
    elif os.path.exists(dir+'/f2py'):
        f2pycmd = 'f2py'

warning="No f2py command found...\n Check your numpy install or directly specify 'f2pycmd'"
if f2pycmd == '':
    sys.exit(warning)

startdir = os.getcwd()
os.chdir('fortranLib')

# For reading G files
cmd = '%s --fcompiler=%s -c -m --opt=-O3 greader readG.f90' % (f2pycmd, compiler)
os.system(cmd)
cmd = 'mv greader.so %s' % startdir
os.system(cmd)

# For the potential field extrapolation
cmd = '%s --fcompiler=%s -c -m --opt=-O3 potential spec.f90' % (f2pycmd, compiler)
os.system(cmd)
cmd = 'mv potential.so %s' % startdir
os.system(cmd)

# For the vtk file format convertion
cmd = '%s --fcompiler=%s -c -m --opt=-O3 vtklib vtkLib.f90' % (f2pycmd, compiler)

os.system(cmd)
cmd = 'mv vtklib.so %s' % startdir
os.system(cmd)
