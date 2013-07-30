import numpy as N
import greader as G
import os
os.environ['F_UFMTENDIAN'] = 'big'

filename = 'G_1.N0m2a'
G.greader.readg(filename)
print G.greader.vp
