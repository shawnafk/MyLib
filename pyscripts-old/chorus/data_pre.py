import sys
import numpy as np
from scipy import optimize
import chorus_lib as cl

with open('loadall.py','r') as f:
    exec(f.read())

#np.save("out_dphi.npy",dphi)

#calc ratio
compensated_sour = compsour*np.exp(-1j*phase_chor)
phase_cpst = -np.angle(compensated_sour)
#JB/JE
#real is J_B while imagnary is J_E
#ratio=np.imag(compensated_sour)/np.real(compensated_sour)
ratio=np.real(compensated_sour)/np.imag(compensated_sour)
#np.save("out_J_ratio.npy",ratio)
#
##more complicited data
##R400=cl.calcR(ratio[400,:])
##R500=cl.calcR(ratio[500,:])
##R600=cl.calcR(ratio[600,:])
#ratio[abs(ratio)>1000]=1000
#ratio[abs(ratio)<0.01]=0.01
#interp_R = cl.calcR_interp(abs(ratio))
#np.save("out_R.npy",interp_R)
##np.save("../data/C_para_R500.npy",R500)
##np.save("../data/C_para_R600.npy",R600)
#
#
##plot non linear dp
#dphidt = (phi_w[1:,:] -  phi_w[:-1,:])/dT
#wall=dphidt+wl
#k = \partial phi/ \partial s




#calculate wave energy
#ampchor= np.linalg.norm(chor,axis=-1)

