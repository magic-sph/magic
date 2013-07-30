import os
from setup import *
from series import *
from radial import *
from graph import *
from surf import *
from spectrum import *
from movie import *
from cyl import *
from compSims import *
from butterfly import *
from bLayers import *
#if os.path.exists('vtklib.so'):
try:
    from potExtra import *
    from movie3D import *
except ImportError,e:
    #print "Potential extrapolation not available"
    #print "Please try to compile the library"
    pass
#    if os.path.exists('potential.so'):
