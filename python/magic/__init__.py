from .setup import buildSo
from .log import *
from .series import MagicTs, AvgField
from .radial import MagicRadial
from .graph import *
from .surf import *
from .spectrum import MagicSpectrum, MagicSpectrum2D
from .movie import *
from .TOreaders import *
from .cyl import *
from .compSims import *
from .butterfly import *
from .bLayers import *
from .checker import *
from .thHeat import *
from .coeff import *
from .radialSpectra import *
from .checkpoint import *
if buildSo:
    try:
        from .potExtra import *
        from .graph2vtk import *
        #from .graph2rst import *
        from .movie3D import *
        from .spectralTransforms import *
        from .potential import *
    except ImportError as e:
        pass
