#!/usr/bin/env python
#
# This script tries to find your matplotlib backend 
# and check whether you can make use of the LaTeX fonts
#
import sys

def getMyMatplotlibEnv():
    try:
        import numpy as np
        np.seterr('raise')
        import matplotlib.pyplot as plt
        import warnings

        backend = ''
        warnings.filterwarnings('error')

        if sys.platform == 'darwin':  #OSX
            try:
                plt.switch_backend('macosx')
                return 'macosx'
            except ImportError:
                return ''
        else:
            # GTKAgg is the default one
            try:
                try:
                    plt.switch_backend('GTKAgg')
                    return 'GTKAgg'
                except Warning:
                    pass
            except ImportError:
                pass

            # Qt4
            try:
                try:
                    plt.switch_backend('Qt4Agg')
                    return 'Qt4Agg'
                except Warning:
                    pass
            except ImportError:
                pass

            # TkAgg is the fallback value
            try:
                try:
                    plt.switch_backend('TkAgg')
                    return 'TkAgg'
                except Warning:
                    pass
            except ImportError:
                pass

            return ''

    except ImportError:
        return ''

print(getMyMatplotlibEnv())
