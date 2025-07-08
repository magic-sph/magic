#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def get_numpy_version():
    
    import numpy as np
    major,minor,patch = np.__version__.split('.')
    version = float(major + '.' + minor)

    return version

if __name__=="__main__":
    print(get_numpy_version())

