import sys
import numpy as np
from scipy import optimize
import chorus_lib as cl
#most unstable
#linear gamma rate
exec(open('/Users/shawn/py_chorus/loadall.py').read())

compchor= chor[...,0] + chor[...,1]*1j
phase_chor = np.angle(compchor)
wall,klin,knl,kall = cl.calcWK(wl,gyro,np.sqrt(wp2),kmode,phase_chor,dT,zpos)

fsp='./data/'
np.save(fsp+"out_freq.npy",wall)
np.save(fsp+"out_klin.npy",klin)
np.save(fsp+"out_knl.npy",knl)
np.save(fsp+"out_ksim.npy",kall)
