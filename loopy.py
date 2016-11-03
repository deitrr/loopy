import numpy as np

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
    kx1[i] = rderiv(rvars,vvars,mass,ibody,nbodies)
    kv1[i] = vderiv(rvars,vvars,mass,ibody,nbodies)
  
  kx2 = np.zeros((nbodies,3))
  kv2 = np.zeros_like(kx2)
  for i in range(nbodies):
    kx2[i] = rderiv(rvars+0.5*kx1*dt,vvars+0.5*kv1*dt,mass,ibody,nbodies)
    kv2[i] = vderiv(rvars+0.5*kx1*dt,vvars+0.5*kv1*dt,mass,ibody,nbodies)
  
  kx3 = np.zeros((nbodies,3))
  kv3 = np.zeros_like(kx3)
  for i in range(nbodies):
    kx3[i] = rderiv(rvars+0.5*kx2*dt,vvars+0.5*kv2*dt,mass,ibody,nbodies)
    kv3[i] = vderiv(rvars+0.5*kx2*dt,vvars+0.5*kv2*dt,mass,ibody,nbodies)
    
  kx4 = np.zeros((nbodies,3))
  kv4 = np.zeros_like(kx4)
  for i in range(nbodies):
    kx4[i] = rderiv(rvars+kx3*dt,vvars+kv3*dt,mass,ibody,nbodies)
    kv4[i] = vderiv(rvars+kx3*dt,vvars+kv3*dt,mass,ibody,nbodies)
    
  for i in range(nbodies):
    dr = 1./6*(kx1[i] + 2*kx2[i] + 2*kx3[i] + kx4[i])*dt
    dv = 1./6*(kv1[i] + 2*kv2[i] + 2*kv3[i] + kv4[i])*dt
    
    rvars[i] += dr
    vvars[i] += dv