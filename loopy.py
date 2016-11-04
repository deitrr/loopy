import numpy as np
import os

k = 0.01720209895
def rderiv(rvars,vvars,mass,ibody,nbodies):
  return vvars[ibody,0],vvars[ibody,1],vvars[ibody,2]
  
def vderiv(rvars,vvars,mass,ibody,nbodies):
  accel = np.zeros(3)
  for jbody in range(nbodies):
    if jbody != ibody:
      r = rvars[ibody] - rvars[jbody]
      rcb = np.sum(r**3)
      accel += k**2*mass[jbody]*r/rcb
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
  output(rvars,vvars,len(mass),time)
  nextout = time+outdt
  
  while time < tfin:
    rk4step(rvars,vvars,mass,len(mass),dt*365.25)
    time += dt
    if time >= nextout:
      output(rvars,vvars,len(mass),time)
      nextout = time+outdt

mass = np.array([1.0,3.7,2.07])
rvars = np.array([[0.32946639,-0.03194335,-0.02485509],[0.2966885,0.06010892,0.00443412], [-0.68947529,-0.0920094,0.00408158]])
vvars = np.array([[0.08906436,0.02777444,0.00653185],[-0.02295717,0.00057109,-0.0021237], [-0.00199171,-0.01443839,0.00064049]])

integrate(rvars,vvars,mass,0.0001,1,10)
