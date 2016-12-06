import numpy as np
import vplot as vpl
import matplotlib.pyplot as plt
import celestmech as cm
import sys

vpldir = 'bv15check'
hnbdir = 'hnbbvtest'

out = vpl.GetOutput(vpldir)
longp1 = cm.fixangles(out.B.ArgP+out.B.LongA)
longp2 = cm.fixangles(out.C.ArgP+out.C.LongA)


t1, a1, e1, i1, ap1, la1, ma1 = np.loadtxt(hnbdir+'/1aei.dat',unpack=True)
t2, a2, e2, i2, ap2, la2, ma2 = np.loadtxt(hnbdir+'/2aei.dat',unpack=True)
lp1 = cm.fixangles(ap1+la1)
lp2 = cm.fixangles(ap2+la2)


plt.subplot(4,2,1)
plt.plot(t1,i1,'.',color='0.5')
plt.plot(out.B.Time,out.B.Inc,'b-')

plt.subplot(4,2,2)
plt.plot(t2,i2,'.',color='0.5')
plt.plot(out.C.Time,out.C.Inc,'b-')


plt.subplot(4,2,3)
plt.plot(t1,la1,'.',color='0.5')
plt.plot(out.B.Time,out.B.LongA,'b.')


plt.subplot(4,2,4)
plt.plot(t2,la2,'.',color='0.5')
plt.plot(out.C.Time,out.C.LongA,'b.')


plt.subplot(4,2,5)
plt.plot(t1,lp1,'.',color='0.5')
plt.plot(out.B.Time,longp1,'b.')


plt.subplot(4,2,6)
plt.plot(t2,lp2,'.',color='0.5')
plt.plot(out.C.Time,longp2,'b.')


plt.subplot(4,2,7)
plt.plot(t1,e1,'.',color='0.5')
plt.plot(out.B.Time,out.B.Eccentricity,'b-')


plt.subplot(4,2,8)
plt.plot(t2,e2,'.',color='0.5')
plt.plot(out.C.Time,out.C.Eccentricity,'b-')


plt.show()