import numpy as np
import os
import celestmech as cm
import matplotlib.pyplot as plt

k = 0.01720209895
def rderiv(rvars,vvars,mass,ibody,nbodies):
  return vvars[ibody,0],vvars[ibody,1],vvars[ibody,2]
  
def vderiv(rvars,vvars,mass,ibody,nbodies):
  accel = np.zeros(3)
  for jbody in range(nbodies):
    if jbody != ibody:
      r = rvars[jbody] - rvars[ibody]
      rcb = (np.sqrt(np.sum(r**2)))**3.0
      accel += k**2*mass[jbody]*r/rcb
  return accel
      
def vderiv2(rvars,vvars,mass,ibody,nbodies):
  #harrington's form in jacobi coords
  rab = rvars[0]
  RA = mass[0]/(mass[0]+mass[1])
  RB = mass[1]/(mass[0]+mass[1])
  rAC = rvars[1] + RB*rvars[0]
  rBC = rvars[1] - RA*rvars[0]
  rACcb = (np.sqrt(np.sum(rAC**2)))**3.0
  rBCcb = (np.sqrt(np.sum(rBC**2)))**3.0
  rABcb = (np.sqrt(np.sum(rvars[0]**2)))**3.0

  if ibody == 0:
    #inner body
    muAB = k**2*(mass[0]+mass[1])
    RC = mass[2]/(mass[0]+mass[1])
    accel = -muAB*(rvars[0]/rABcb+RC*(rAC/rACcb-rBC/rBCcb))
  else:
    #outer body
    muABC = k**2*(mass[0]+mass[1]+mass[2])
    accel = -muABC*(RA*rAC/rACcb + RB*rBC/rBCcb)    
  
  return accel

def rk4step(rvars, vvars, mass, nbodies, dt):
#rvars and vvars are the positions and speeds of bodies, with shape (nbodies,3)
#mass has shape (nbodies)
  kx1 = np.zeros((nbodies,3))
  kv1 = np.zeros_like(kx1)
  for ibody in range(nbodies):
    kx1[ibody] = rderiv(rvars,vvars,mass,ibody,nbodies)
    kv1[ibody] = vderiv(rvars,vvars,mass,ibody,nbodies)
  
  kx2 = np.zeros((nbodies,3))
  kv2 = np.zeros_like(kx2)
  for ibody in range(nbodies):
    kx2[ibody] = rderiv(rvars+0.5*kx1*dt,vvars+0.5*kv1*dt,mass,ibody,nbodies)
    kv2[ibody] = vderiv(rvars+0.5*kx1*dt,vvars+0.5*kv1*dt,mass,ibody,nbodies)
  
  kx3 = np.zeros((nbodies,3))
  kv3 = np.zeros_like(kx3)
  for ibody in range(nbodies):
    kx3[ibody] = rderiv(rvars+0.5*kx2*dt,vvars+0.5*kv2*dt,mass,ibody,nbodies)
    kv3[ibody] = vderiv(rvars+0.5*kx2*dt,vvars+0.5*kv2*dt,mass,ibody,nbodies)
    
  kx4 = np.zeros((nbodies,3))
  kv4 = np.zeros_like(kx4)
  for ibody in range(nbodies):
    kx4[ibody] = rderiv(rvars+kx3*dt,vvars+kv3*dt,mass,ibody,nbodies)
    kv4[ibody] = vderiv(rvars+kx3*dt,vvars+kv3*dt,mass,ibody,nbodies)
    
  for i in range(nbodies):
    dr = 1./6*(kx1[i] + 2*kx2[i] + 2*kx3[i] + kx4[i])*dt
    dv = 1./6*(kv1[i] + 2*kv2[i] + 2*kv3[i] + kv4[i])*dt
    
    rvars[i] += dr
    vvars[i] += dv
    
def output(rvars, vvars, nbodies, t):
  for ibody in range(nbodies):
    if os.path.exists('b%02d.out'%ibody): 
      f = open('b%02d.out'%ibody,'a')
    else:
      f = open('b%02d.out'%ibody,'w')
    f.write('%#.8f %#.8f %#.8f %#.8f %#.8f %#.8f %#.8f\n'%(t,rvars[ibody,0],rvars[ibody,1],rvars[ibody,2],vvars[ibody,0],vvars[ibody,1],vvars[ibody,2]))
    f.close()

def integrate(rvars,vvars,mass,dt,outdt,tfin):
  time = 0
  nbodies = len(mass)
  
  output(rvars,vvars,nbodies,time)
  nextout = time+outdt
  
  while time < tfin:
    rk4step(rvars,vvars,mass,nbodies,dt*365.25)
    time += dt
    if time >= nextout:
      output(rvars,vvars,nbodies,time)
      nextout = time+outdt

os.system('rm *.out')
# mass = np.array([1.0,3.7,2.07])
# # rvars = np.array([[0.32946639,-0.03194335,-0.02485509],[0.2966885,0.06010892,0.00443412], [-0.68947529,-0.0920094,0.00408158]])
# # vvars = np.array([[0.08906436,0.02777444,0.00653185],[-0.02295717,0.00057109,-0.0021237], [-0.00199171,-0.01443839,0.00064049]])
# # 
# # integrate(rvars,vvars,mass,0.0001,1,10)
# # 
# # test out harrington form in jacobi coords
# # yr = 365.25
# # P = 79.91*yr
# # 
# # ms = 1.
# # m = np.array([3.7,2.07])
# # # aAB = cm.per2semi(ms,m[0],P)
# # 
# a = np.array([0.100027283,1.04164747])
# e = np.array([0.08,0.269])
# inc = np.array([17.65,2.54])
# argp = np.array([0.0,270])
# longa = np.array([0,180])
# meana = np.array([99.9,66.6])
# 
# xBa, vBa = cm.osc2x(mass[0],mass[1],a[0],e[0],inc[0],argp[0],longa[0],meana[0])
# xB,vB,xA,vA = cm.astro2bary(mass[0],mass[1],xBa,vBa)
# xCoM = xBa*mass[1]/(mass[0]+mass[1])
# vCoM = vBa*mass[1]/(mass[0]+mass[1])
# 
# m2 = mass[2]*mass[0]/(mass[0]+mass[1])
# print(mass[0]+m2)
# xC, vC = cm.osc2x(mass[0],m2,a[1],e[1],inc[1],argp[1],longa[1],meana[1])
# 
# xa = np.hstack([xBa,xC])
# va = np.hstack([vBa,vC])
# 
# xb, vb, xs, vs = cm.astro2bary(mass[0],mass[1:],xa,va)
# 
# rvars = np.hstack([xs,xb]).T
# vvars = np.hstack([vs,vb]).T

mass = np.array([1.0,3.7,2.07])
a = np.array([0.100027283,1.04164747])
e = np.array([0.08,0.269])
inc = np.array([17.65,2.54])
argp = np.array([0.0,270])
longa = np.array([0,180])
meana = np.array([99.9,66.6])

xBa, vBa = cm.osc2x(mass[0],mass[1],a[0],e[0],inc[0],argp[0],longa[0],meana[0])
xB,vB,xA,vA = cm.astro2bary(mass[0],mass[1],xBa,vBa)
xCoM = xBa*mass[1]/(mass[0]+mass[1])
vCoM = vBa*mass[1]/(mass[0]+mass[1])

# m2 = mass[2]*mass[0]/(mass[0]+mass[1])
# print(mass[0]+m2)
xC, vC = cm.osc2x(mass[0]+mass[1],mass[2],a[1],e[1],inc[1],argp[1],longa[1],meana[1])

rCOM = (mass[0]*xA+mass[1]*xB+mass[2]*xC)/np.sum(mass)
vCOM = (mass[0]*vA+mass[1]*vB+mass[2]*vC)/np.sum(mass)

xb = np.hstack([xA,xB,xC])-rCOM
vb = np.hstack([vA,vB,vC])-vCOM

rvars = xb.T
vvars = vb.T
# import pdb; pdb.set_trace()
# plt.subplot(2,2,1)
# plt.plot(xb[0][1],xb[1][1],'b.')
# plt.plot(xb[0][0],xb[1][0],'k.')
# plt.plot(xb[0][2],xb[1][2],'m.')
# plt.plot([0],[0],'k+')
# 
# plt.subplot(2,2,2)
# plt.plot(xb[0][1],xb[2][1],'b.')
# plt.plot(xb[0][0],xb[2][0],'k.')
# plt.plot(xb[0][2],xb[2][2],'m.')
# plt.plot([0],[0],'k+')
# 
# plt.subplot(2,2,3)
# plt.plot(xb[1][1],xb[2][1],'b.')
# plt.plot(xb[1][0],xb[2][0],'k.')
# plt.plot(xb[1][2],xb[2][2],'m.')
# plt.plot([0],[0],'k+')
# plt.show()
integrate(rvars,vvars,mass,0.0001,1,500)

# m = np.array([9.54502e-4,2.857878e-4])
# a = np.array([5.20336301,9.53707032])
# e = np.array([0.04839266,0.05415060])
# inc = np.array([1.30530,2.48446])
# longa = np.array([100.556,113.715])
# argp = longa - np.array([14.75385,92.43194])
# meana = np.array([34.40438,49.94432])
# 
# xa, va = cm.osc2x(1,m,a,e,inc,argp,longa,meana)
# xb, vb, xs, vs = cm.astro2bary(1,m,xa,va)
# 
# mass = np.array([1,9.54502e-4,2.857878e-4])
# rvars = np.hstack([xs,xb]).T
# vvars = np.hstack([vs,vb]).T
# integrate(rvars,vvars,mass,0.1,1000,100000)


#fucking 2 body problem 
# mass = np.array([1.0,3.7])
# rvars = np.array([[0.32946639,-0.03194335,-0.02485509],[0.2966885,0.06010892,0.00443412]])
# vvars = np.array([[0.08906436,0.02777444,0.00653185],[-0.02295717,0.00057109,-0.0021237]])
# 
# integrate(rvars,vvars,mass,0.0001,1,10)