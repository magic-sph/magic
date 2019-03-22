from .setup import buildSo
from .log import *
from .series import *
from .radial import *
from .graph import *
from .surf import *
from .spectrum import *
from .movie import *
from .TOreaders import *
from .cyl import *
from .compSims import *
from .butterfly import *
from .bLayers import *
from .thHeat import *
from .coeff import *
from .radialSpectra import *
if buildSo:
    try:
        from .potExtra import *
        from .graph2vtk import *
        from .graph2rst import *
        from .movie3D import *
        from .spectralTransforms import *
        from .potential import *
    except ImportError as e:
        pass
