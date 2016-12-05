import numpy as np
import celestmech as cm
import matplotlib.pyplot as plt

yr = 365.25
P = 79.91*yr

ms = 1.
m = np.array([3.7,2.07])
# aAB = cm.per2semi(ms,m[0],P)

a = np.array([0.100027283,1.04164747])
e = np.array([0.08,0.269])
inc = np.array([17.65,2.54])
argp = np.array([0.0,270])
longa = np.array([0,180])
meana = np.array([99.9,66.6])

xBa, vBa = cm.osc2x(ms,m[0],a[0],e[0],inc[0],argp[0],longa[0],meana[0])
xB,vB,xA,vA = cm.astro2bary(ms,m[0],xBa,vBa)
xCoM = xBa*m[0]/(ms+m[0])
vCoM = vBa*m[0]/(ms+m[0])

#m2 = m[1]*ms/(ms+m[0])
xC, vC = cm.osc2x(ms+m[0],m[1],a[1],e[1],inc[1],argp[1],longa[1],meana[1])
xCa = xC+xCoM
vCa = vC+vCoM

xa = np.hstack([xBa,xCa])
va = np.hstack([vBa,vCa])

xbary,vbary,xc,vc = cm.astro2bary(ms,m,xa,va)

import pdb; pdb.set_trace()

cm.print_coords('tmp.hnb','trial.hnb','hnbbvtest',m,xa,va*365.25,ms=ms)
plt.subplot(2,2,1)
plt.plot(vbary[0,0],vbary[1,0],'b.')
plt.plot(vbary[0,1],vbary[1,1],'r.')
plt.plot(vc[0],vc[1],'ks')
plt.plot(0,0,'k+')

plt.subplot(2,2,2)
plt.plot(vbary[0,0],vbary[2,0],'b.')
plt.plot(vbary[0,1],vbary[2,1],'r.')
plt.plot(vc[0],vc[2],'ks')
plt.plot(0,0,'k+')

plt.subplot(2,2,3)
plt.plot(vbary[1,0],vbary[2,0],'b.')
plt.plot(vbary[1,1],vbary[2,1],'r.')
plt.plot(vc[1],vc[2],'ks')
plt.plot(0,0,'k+')


plt.show()