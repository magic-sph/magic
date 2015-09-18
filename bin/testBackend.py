#!/bin/env python
#
# This script tries to find your matplotlib backend 
# and check whether you can make use of the LaTeX fonts
#

try:
    import matplotlib as mpl
    import matplotlib.pyplot as P
    import sys


    fontNames = [f.name for f in mpl.font_manager.fontManager.afmlist]
    cmIsFound = False
    for font in fontNames:
        if 'Computer' in font:
            print(font)
            cmIsFound = True

    if cmIsFound:
        print('LaTex fonts are available')

    try:
        P.switch_backend('GTKAgg')
        sys.exit('GTKAgg has been selected')
    except ImportError:
        pass

    try:
        P.switch_backend('Qt4Agg')
        sys.exit('Qt4Agg has been selected')
    except ImportError:
        pass

    try:
        P.switch_backend('TkAgg')
        sys.exit('TkAgg has been selected')
    except ImportError:
        pass

    sys.exit('No backend found')


except ImportError:
    import sys
    sys.exit('Bye bye')
