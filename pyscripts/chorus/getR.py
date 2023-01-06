import sys
import numpy as np
from scipy import optimize
import chorus_lib as cl
exec(open('/Users/shawn/py_chorus/loadall.py').read())
compchor= chor[...,0] + chor[...,1]*1j
compsour= sour[...,0] + sour[...,1]*1j
phase_chor = np.angle(compchor)
phase_sour = np.angle(compsour)
ampchor= np.linalg.norm(chor,axis=-1)
ampsour= np.linalg.norm(sour,axis=-1)

compensated_sour = compsour*np.exp(-1j*phase_chor)
phase_cpst = -np.angle(compensated_sour)
#JB/JE
#real is J_B while imagnary is J_E
#ratio=np.imag(compensated_sour)/np.real(compensated_sour)
ratio=np.real(compensated_sour)/np.imag(compensated_sour)
ratio[abs(ratio)>1000]=1000
ratio[abs(ratio)<0.01]=0.01
interp_R = cl.calcR_interp(abs(ratio))
fsp='./data/'                                          
np.save(fsp+"out_R.npy",interp_R)



