from __future__ import print_function
import unittest
import numpy as np
import glob
import os
import time
import shutil
import subprocess as sp
import argparse
import sys


def readData(file):
    return np.loadtxt(file)

def verifNum():
    datRef = readData('reference.out')
    datTmp = readData('e_kin.test')
    np.testing.assert_allclose(datRef, datTmp, rtol=1e-8, atol=1e-20)

if __name__ == '__main__':
    verifNum()
    sys.exit()
