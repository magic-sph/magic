#!/usr/bin/env python
import subprocess as sp
from magic import TOMovie

to = TOMovie(file='TO_mov.start', iplot=False)

out = 'tmp'
file = open(out, 'w')
st = '%.4f %.4f %.4f %.4f %.4f %.4f %.4f' % ( to.asVphi[0, 13, 3], 
     to.rey[1, 21, 22], to.adv[1, 52, 11], to.visc[0, 12, 25], 
     to.lorentz[0, 73, 30], to.coriolis[1, 33, 3], to.dtVp[1, 88, 7] )

# Write output for TO files
file.write(st+'\n')
file.close()

# Cat e_kin.test + misc
with open('e_kin.test', 'w') as outFile:
    sp.call(['cat', 'geos.start', out], stdout=outFile)

sp.call(['rm', 'tmp'])
