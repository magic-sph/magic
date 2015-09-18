#!/usr/bin/env python
#
# This script tries to find your matplotlib backend 
# and check whether you can make use of the LaTeX fonts
#
import sys

def getMyMatplotlibEnv():
    try:
        import numpy as N
        N.seterr('raise')
        import matplotlib.pyplot as P
        import warnings

        backend = ''
        warnings.filterwarnings('error')

        try:
            try:
                P.switch_backend('GTKAgg')
                return 'GTKAgg'
            except Warning:
                pass
        except ImportError:
            pass

        try:
            try:
                P.switch_backend('Qt5Agg')
                return 'Qt5Agg'
            except Warning:
                pass
        except ImportError:
            pass

        try:
            try:
                P.switch_backend('Qt4Agg')
                return 'Qt4Agg'
            except Warning:
                pass
        except ImportError:
            pass

        try:
            try:
                P.switch_backend('TkAgg')
                return 'TkAgg'
            except Warning:
                pass
        except ImportError:
            pass

        return ''

    except ImportError:
        return ''

print(getMyMatplotlibEnv())
