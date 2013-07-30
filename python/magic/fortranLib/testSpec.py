from magic import MagicGraph
from magic.libmagic import symmetrize
import pylab as P
import numpy as N
from potential import *


def proj(theta, phi):
    x = 2.*N.sqrt(2.)*N.sin(theta)*N.sin(phi/2.)/N.sqrt(1.+N.sin(theta)*N.cos(phi/2.))
    y = N.sqrt(2.)*N.cos(theta)/N.sqrt(1.+N.sin(theta)*N.cos(phi/2.))
    return x, y

g = MagicGraph()
rcmb = g.radius[0]
rsurf= rcmb
radratio=rcmb/rsurf
brCMB = g.Br[..., 0]
btCMB = g.Btheta[..., 0]
bpCMB = g.Bphi[..., 0]

phi = N.linspace(-N.pi, N.pi, g.nphi)
theta = g.colatitude
ttheta, pphi = N.meshgrid(theta, phi)

x, y = proj(ttheta, pphi)
anlc = N.fft.fft(brCMB, axis=0)/(4.*N.pi*g.npI)
brm, btm, bpm = extrapolate(anlc, radratio, g.minc)

brsurf = N.fft.ifft(brm, axis=0)*g.npI
brsurf = brsurf.real
btsurf = N.fft.ifft(btm, axis=0)*g.npI
btsurf = btsurf.real
bpsurf = N.fft.ifft(bpm, axis=0)*g.npI
bpsurf = bpsurf.real

print 'BrCMB: %.8e %.8e' % (brCMB.max(), brCMB.min())
print 'reconstruct BrCMB: %.8e %.8e' % (brsurf.max(), brsurf.min())
print 'BtCMB: %.8e %.8e' % (btCMB.max(), btCMB.min())
print 'reconstruct BtCMB: %.8e %.8e' % (btsurf.max(), btsurf.min())
print 'BpCMB: %.8e %.8e' % (bpCMB.max(), bpCMB.min())
print 'reconstruct BpCMB: %.8e %.8e' % (bpsurf.max(), bpsurf.min())

fig = P.figure(figsize=(12,4))
P.subplots_adjust(right=0.99, left=0.01, top=0.99, bottom=0.01, wspace=0.02,
                  hspace=0.02)
ax = fig.add_subplot(231, frameon=False)
im = ax.contourf(x,y,symmetrize(brCMB, g.minc), 16, cmap=P.get_cmap('RdYlBu_r'))
P.axis('off')
ax = fig.add_subplot(232, frameon=False)
im = ax.contourf(x,y,symmetrize(btCMB, g.minc), 16, cmap=P.get_cmap('RdYlBu_r'))
P.axis('off')
ax = fig.add_subplot(233, frameon=False)
im = ax.contourf(x,y,symmetrize(bpCMB, g.minc), 16, cmap=P.get_cmap('RdYlBu_r'))
P.axis('off')
ax = fig.add_subplot(234, frameon=False)
im = ax.contourf(x,y,symmetrize(brsurf, g.minc), 16, cmap=P.get_cmap('RdYlBu_r'))
P.axis('off')
ax = fig.add_subplot(235, frameon=False)
im = ax.contourf(x,y,symmetrize(btsurf, g.minc), 16, cmap=P.get_cmap('RdYlBu_r'))
P.axis('off')
ax = fig.add_subplot(236, frameon=False)
im = ax.contourf(x,y,symmetrize(bpsurf, g.minc), 16, cmap=P.get_cmap('RdYlBu_r'))
P.axis('off')
P.show()
