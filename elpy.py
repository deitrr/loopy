import numpy as np
import celestmech as cm
import matplotlib.pyplot as plt

mass = np.array([1.0,3.7,2.07])
t0, x0, y0, z0, u0, v0, w0 = np.loadtxt('b00.out',unpack=True)
t1, x1, y1, z1, u1, v1, w1 = np.loadtxt('b01.out',unpack=True)
t2, x2, y2, z2, u2, v2, w2 = np.loadtxt('b02.out',unpack=True)

r1 = np.vstack([x1-x0,y1-y0,z1-z0])
vv1 = np.vstack([u1-u0,v1-v0,w1-w0])
m1 = np.zeros(len(t0))+mass[1]
a1,e1,i1,ap1,la1,ma1 = cm.x2osc(mass[0],m1,r1,vv1)
lp1 = cm.fixangles(ap1+la1)

xAB = (mass[0]*x0 + mass[1]*x1)/(mass[0]+mass[1])
yAB = (mass[0]*y0 + mass[1]*y1)/(mass[0]+mass[1])
zAB = (mass[0]*z0 + mass[1]*z1)/(mass[0]+mass[1])
uAB = (mass[0]*u0 + mass[1]*u1)/(mass[0]+mass[1])
vAB = (mass[0]*v0 + mass[1]*v1)/(mass[0]+mass[1])
wAB = (mass[0]*w0 + mass[1]*w1)/(mass[0]+mass[1])

r2 = np.vstack([x2-xAB,y2-yAB,z2-zAB])
vv2 = np.vstack([u2-uAB,v2-vAB,w2-wAB])
m2 = np.zeros(len(t0))+mass[2]
a2,e2,i2,ap2,la2,ma2 = cm.x2osc(mass[0]+mass[1],m2,r2,vv2)
lp2 = cm.fixangles(ap2+la2)

plt.subplot(4,2,1)
plt.plot(t1,i1,'k.')

plt.subplot(4,2,2)
plt.plot(t2,i2,'k.')

plt.subplot(4,2,3)
plt.plot(t1,la1,'k.')

plt.subplot(4,2,4)
plt.plot(t2,la2,'k.')

plt.subplot(4,2,5)
plt.plot(t1,lp1,'k.')

plt.subplot(4,2,6)
plt.plot(t2,lp2,'k.')

plt.subplot(4,2,7)
plt.plot(t1,e1,'k.')

plt.subplot(4,2,8)
plt.plot(t2,e2,'k.')

plt.show()