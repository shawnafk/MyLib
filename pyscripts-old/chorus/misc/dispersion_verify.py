import numpy as np
wc=np.loadtxt("../data/ori/gyro.out")[1:]
wp2=np.loadtxt("../data/ori/wp2.out")[1:]
wall = np.load("../data/C_freq.npy")[:,1:-1]
k0 = np.load("../data/C_klin.npy")
knl = np.load("../data/C_knl.npy")
kall = np.load("../data/C_ksim.npy")
kall[kall<0.1]=np.nan
wl=0.061457
from matplotlib import pyplot as plt
#Omura's Dispersion relation
#from 29 and J_E B_w find omega_i
#set J_B = J_E; known Vg (must be cold, not cold, is from chired w, k in Vg is elimited) calculate, so right term only need omega_i and omega_w


chor= np.load("../data/T_chor.npy")
sour= np.load("../data/T_sour.npy")
compchor= chor[...,0] + chor[...,1]*1j
compsour= sour[...,0] + sour[...,1]*1j

ampchor= np.linalg.norm(chor,axis=-1)
#data B phi of wave and hole and phi difference
phi_w =(np.unwrap(np.angle(compchor),axis=0))
phi_h =(np.unwrap(np.angle(compsour),axis=0))
dphi = phi_w - phi_h
compensated_sour = compsour*np.exp(-1j*phi_w)
#JB/JE
#real is J_B while imagnary is J_E
#ratio=np.imag(compensated_sour)/np.real(compensated_sour)



JB=-np.real(compensated_sour) 
JE=-np.imag(compensated_sour) 
Vg0 = np.loadtxt("../data/ori/vg.out")

#first one k from cold plasma dispersion relation
k_cold = np.sqrt((wall**2 + wall * wp2/(wc - wall)))

#using JB to get k: 
k_JB = np.sqrt((wall**2 + wall * wp2/(wc - wall)) + (JB/ampchor)[:-1,1:-1])
#using JE to get k
k_JE = np.sqrt((wall**2 + wall * wp2/(wc - wall)) + (JE/ampchor)[:-1,1:-1])

#using omega_i to get k
#omega_i gamma linear
gm0 = np.loadtxt("../data/ori/gamma.out")[:-1]
gm=gm0
k_gm0p = np.sqrt(wall**2 * (1 + 2*gm/wall) + wall*wp2/(wc-wall) *(1+ wc/(wc-wall)*gm/wall))
k_gm0m = np.sqrt(wall**2 * (1 - 2*gm/wall) + wall*wp2/(wc-wall) *(1- wc/(wc-wall)*gm/wall))


#omega_i from JE and kall.Vg
omega_i = (JE/ampchor)[:-1,1:-1] /(2* wall + wc*wp2/(wc-wall)**2 ) 
gm=omega_i
k_Wip = np.sqrt(wall**2 * (1 + 2*gm/wall) + wall*wp2/(wc-wall) *(1+ wc/(wc-wall)*gm/wall))
k_Wim = np.sqrt(wall**2 * (1 - 2*gm/wall) + wall*wp2/(wc-wall) *(1- wc/(wc-wall)*gm/wall))


#verify K is consistent with this
right = kall**2 - wall**2 - wall*wp2/(wc-wall) 
left =  (JB/ampchor)[:-1,1:-1]
residual_sim = right + left

zpos = np.load("../data/T_zpos.npy")
zrange=np.arange(400,900)
ts=600
f,(ax1,ax2)=plt.subplots(2,1)
ax1.plot(zpos[zrange],k_cold[ts,zrange],label='cold plasma')
ax1.plot(zpos[zrange],knl[ts,zrange],label='nonlinear')
ax1.plot(zpos[zrange],kall[ts,zrange],label='simulation')
ax1.plot(zpos[zrange],k_JB[ts,zrange],label=r'$k^+~from~J_B$')

ax1.plot(zpos[zrange],k_gm0p[ts,zrange],'m',label=r'$k^+~from~\gamma_L$')
ax1.plot(zpos[zrange],k_gm0m[ts,zrange],'m--',label=r'$k^-~from~\gamma_L$')

ax1.plot(zpos[zrange],k_Wip[ts,zrange],'k',label=r'$k^+~from~J_E~$')
ax1.plot(zpos[zrange],k_Wim[ts,zrange],'k--',label=r'$k^-~from~J_E~$')
ax1.legend(fontsize=18)
#相对simulation
R1 = abs(k_cold[ts,zrange] - kall[ts,zrange])/kall[ts,zrange]
R2 = abs(knl[ts,zrange] - kall[ts,zrange])/kall[ts,zrange]
R3 = abs(k_JB[ts,zrange] - kall[ts,zrange])/kall[ts,zrange]

R4 = abs(k_gm0p[ts,zrange] - kall[ts,zrange])/kall[ts,zrange]
R5 = abs(k_gm0m[ts,zrange] - kall[ts,zrange])/kall[ts,zrange]

R6 = abs(k_Wip[ts,zrange] - kall[ts,zrange])/kall[ts,zrange]
R7 = abs(k_Wim[ts,zrange] - kall[ts,zrange])/kall[ts,zrange]

ax2.plot(zpos[zrange],R1,label='cold plasma')
ax2.plot(zpos[zrange],R2,label='nonlinear')
ax2.plot(zpos[zrange],R3,label=r'$k^+~from~J_B$')

ax2.plot(zpos[zrange],R4,'m',label=r'$k^+~from~\gamma_L$')
ax2.plot(zpos[zrange],R5,'m--',label=r'$k^-~from~\gamma_L$')

ax2.plot(zpos[zrange],R6,'k',label=r'$k^+~from~J_E~$')
ax2.plot(zpos[zrange],R7,'k--',label=r'$k^-~from~J_E~$')
ax2.legend(fontsize=20)
plt.show()

#plt.figure()
#ratio = abs(JE)/abs(JB)
#imshow(ratio[100:600,500:-10])
f,(ax1,ax2,ax3)=plt.subplots(3,1)
ax1.plot(abs(JB[200,:]),label=r'$J_B$')
ax1.plot(abs(JE[200,:]),label=r'$J_E$')
ax1.legend(fontsize='large')

ax2.plot(abs(JB[400,:]),label=r'$J_B$')
ax2.plot(abs(JE[400,:]),label=r'$J_E$')
ax2.legend(fontsize='large')

ax3.plot(abs(JB[600,:]),label=r'$J_B$')
ax3.plot(abs(JE[600,:]),label=r'$J_E$')
ax3.legend(fontsize='large')
