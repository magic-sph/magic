#!/usr/bin/env python
import subprocess as sp
from magic import MagicCoeffCmb, MagicCoeffR

out = 'e_kin.test'
file = open(out, 'w')

# Classical CMB file
cmb = MagicCoeffCmb(tag='start',  iplot=False)
st = '%.5e %.5e %.5e %.5e %.5e' % ( cmb.glm[0, cmb.idx[1, 0]], 
     cmb.glm[2, cmb.idx[1, 0]], cmb.hlm[-1, cmb.idx[5, 0]], cmb.glm[2, cmb.idx[5, 4]], 
     cmb.hlm[3, cmb.idx[5, 4]] )
file.write(st+'\n')

# Secular variation CMB file
cmb = MagicCoeffCmb(tag='start',  iplot=False, sv=True)
st = '%.5e %.5e %.5e %.5e %.5e' % ( cmb.glm[0, cmb.idx[1, 0]], 
     cmb.glm[2, cmb.idx[1, 0]], cmb.hlm[-1, cmb.idx[5, 0]], cmb.glm[2, cmb.idx[5, 4]], 
     cmb.hlm[3, cmb.idx[5, 4]] )
file.write(st+'\n')

# Time-averaged CMB file
cmb = MagicCoeffCmb(tag='start',  iplot=False, ave=True)
st = '%.5e %.5e %.5e %.5e' % ( cmb.glm[0, cmb.idx[1, 0]], 
     cmb.hlm[0, cmb.idx[5, 0]], cmb.glm[0, cmb.idx[5, 4]], cmb.hlm[0, cmb.idx[5, 4]] )
file.write(st+'\n')


# Coeff at depth
coeff = MagicCoeffR(tag='start',  field='V', r=1, iplot=False)
st = '%.5e %.5e %.5e %.5e' % ( coeff.wlm[0, coeff.idx[4, 4]].real, 
     coeff.wlm[0, coeff.idx[4, 4]].imag, coeff.zlm[0, coeff.idx[1, 0]].real, coeff.wlm[0, coeff.idx[2, 0]].imag )
file.write(st+'\n')

coeff = MagicCoeffR(tag='start',  field='V', r=2, iplot=False)
st = '%.5e %.5e %.5e %.5e' % ( coeff.wlm[0, coeff.idx[4, 4]].real, 
     coeff.wlm[0, coeff.idx[4, 4]].imag, coeff.zlm[0, coeff.idx[1, 0]].real, coeff.wlm[0, coeff.idx[2, 0]].imag )
file.write(st+'\n')

coeff = MagicCoeffR(tag='start',  field='B', r=1, iplot=False)
st = '%.5e %.5e %.5e %.5e' % ( coeff.wlm[0, coeff.idx[1, 0]].real, 
     coeff.wlm[0, coeff.idx[1, 0]].imag, coeff.zlm[0, coeff.idx[5, 0]].real, coeff.wlm[0, coeff.idx[5, 4]].imag )
file.write(st+'\n')

coeff = MagicCoeffR(tag='start',  field='B', r=2, iplot=False)
st = '%.5e %.5e %.5e %.5e' % ( coeff.wlm[0, coeff.idx[1, 0]].real, 
     coeff.wlm[0, coeff.idx[1, 0]].imag, coeff.zlm[0, coeff.idx[4, 0]].real, coeff.wlm[0, coeff.idx[5, 4]].imag )
file.write(st+'\n')

file.close()
