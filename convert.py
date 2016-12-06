import numpy as np
import celestmech as cm

ms = 1
m = np.array([3.7,2.07])
#read in bodycentric coords
tb, xb, yb, zb, ub, vb, wb = np.loadtxt('1.dat',unpack=True)
tc, xc, yc, zc, uc, vc, wc = np.loadtxt('2.dat',unpack=True)

f = open('1aei.dat','w')
g = open('2aei.dat','w')

#m2 = m[1]*ms/(ms+m[0])
for i in range(len(tb)):
  xB = np.array([xb[i],yb[i],zb[i]])
  vB = np.array([ub[i],vb[i],wb[i]])
  
  xCoM = xB*m[0]/(ms+m[0])
  vCoM = vB*m[0]/(ms+m[0])

  xC = np.array([xc[i],yc[i],zc[i]]) - xCoM
  vC = np.array([uc[i],vc[i],wc[i]]) - vCoM
  
  ab, eb, ib, apb, lab, mab = cm.x2osc(ms,m[0],xB,vB/365.25)
  ac, ec, ic, apc, lac, mac = cm.x2osc(ms+m[0],m[1],xC,vC/365.25)
  
  f.write('%#.1f %#.8f %#.8f %#.8f %#.8f %#.8f %#.8f\n'%(tb[i],ab,eb,ib,apb,lab,mab))
  g.write('%#.1f %#.8f %#.8f %#.8f %#.8f %#.8f %#.8f\n'%(tc[i],ac,ec,ic,apc,lac,mac))
    
f.close()
g.close()