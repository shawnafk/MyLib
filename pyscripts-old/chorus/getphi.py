import sys
import numpy as np
from scipy import optimize
import chorus_lib as cl
exec(open('/Users/shawn/py_chorus/loadall.py').read())
fsp='./data/'                                          
compchor= chor[...,0] + chor[...,1]*1j
compsour= sour[...,0] + sour[...,1]*1j
phase_chor = np.angle(compchor)
phase_sour = np.angle(compsour)

np.save(fsp+"out_phi_w0.npy",phase_chor)
np.save(fsp+"out_phi_h0.npy",phase_sour)

#data B phi of wave and hole and phi difference
#time phase should be continueous
phi_w =(np.unwrap(phase_chor[:,::-1],axis=1))[:,::-1]
phi_h =(np.unwrap(phase_sour[:,::-1],axis=1))[:,::-1]

np.save(fsp+"out_phi_w.npy",phi_w)
np.save(fsp+"out_phi_h.npy",phi_h)

#space start from right
phi_w =np.unwrap(phi_w[:,::-1],axis=1)
phi_h =np.unwrap(phi_h[:,::-1],axis=1)
#start from left at a given point
#lstart=350
#phi_w =np.unwrap(phi_w[:,lstart:],axis=1)
#phi_h =np.unwrap(phi_h[:,lstart:],axis=1)
#dphi: only unwrap the space direction
#phi_w =(np.unwrap(phase_chor[:,250:],axis=1))
#phi_h =(np.unwrap(phase_sour[:,250:],axis=1))
dphi = phi_w - phi_h

np.save(fsp+"out_dphi.npy",dphi)
