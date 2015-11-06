#!/usr/bin/env python
from magic import MagicGraph, Movie

gr = MagicGraph(ivar=1)

out = 'e_kin.test'
file = open(out, 'w')
st = '%.4f %.4f %.4f %.4f %.4f %.4f %.4f' % ( gr.entropy[10, 13, 3], 
     gr.Br[33, 0, 2], gr.Btheta[3, 11, 11], gr.Bphi[34, 12, 25], 
     gr.vr[11, 15, 2], gr.vtheta[33, 33, 3], gr.vphi[32, 33, 7] )

# Write output for graphic files
file.write(st+'\n')

av = Movie(file='AV_mov.start', iplot=False)
ahf = Movie(file='AHF_mov.start', iplot=False)
brcmb = Movie(file='Br_CMB_mov.start', iplot=False)
vtr = Movie(file='Vt_R=C1_mov.start', iplot=False)
vreq = Movie(file='Vr_EQU_mov.start', iplot=False)
teq = Movie(file='T_EQU_mov.start', iplot=False)
vortz = Movie(file='VorZ_EQU_mov.start', iplot=False)
hel = Movie(file='HE_R=C2_mov.start', iplot=False)
st = '%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f' % (av.data[1, 121, 12], 
     ahf.data[0, 99, 33], brcmb.data[1, 47, 128], vtr.data[0, 33, 87],
     vreq.data[0, 28, 31], teq.data[1, 29, 29], vortz.data[1, 1, 1],
     hel.data[0, 3, 4])

# Write output for movie files
file.write(st)

file.close()
