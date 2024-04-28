from chorus_lib import *
##visualization 
#wtr = 0.5
#S =  - 0.4
##figure1: phase structure
#fig,(ax1,ax2)=plt.subplots(2,1)
#
##plot potential cos(\xi) -S \xi
#xi=np.linspace(0,2*np.pi,100)
#ax1.plot(xi,np.cos(xi) - S *xi)
#ax1.set_xlim(0,2*np.pi)
#
##plot phase contour
#zeta = np.linspace(0,2*np.pi,100)
#theta = np.linspace(-1,1,50) * wtr * 2 
#Z,T = np.meshgrid(zeta,theta)
#phase = T**2 + 2* wtr**2 * (np.cos(Z) -  S *  Z)
#ax2.contour(Z,T,phase,50)
#
##plot X and O point
#zeta_o,zeta_x = find_zeta(wtr,S)
#distance = np.abs(zeta_o - zeta_x)
#ax2.scatter(zeta_x,0,s=100,marker='x',lw=1.5)
#ax2.scatter(zeta_o,0,s=100,marker='o',lw=1.5)
#
##plot upper and lower theta boundary of the trapping region
#zeta=np.sort([zeta_o,zeta_x])
#zetab = np.linspace(zeta[0],zeta[1],100)
#theta_b1 = wtr*np.sqrt(2*(np.cos(zeta_x) - np.cos(zetab) + S*(zetab-zeta_x)) )
#theta_b2 = -theta_b1
#ax2.plot(zetab,theta_b1,'r',lw=3)
#ax2.plot(zetab,theta_b2,'k',lw=3)
#
#
#plt.subplots_adjust(hspace=0)
#ax1.axvline(zeta_x,color='k')
#ax1.axvline(zeta_o,color='k')
#ax2.axvline(zeta_x,color='k')
#ax2.axvline(zeta_o,color='k')
#
##current E as function of S
#S= np.linspace(-1,0,100)
#J=np.array([resonant_J(_S,wtr) for _S in S])
#JE=J[:,0]
#JB=J[:,1]
#
#fig,ax=plt.subplots()
#ax.plot(S,-JE)
#ax.plot(S,-JB)
#plt.show()
S= np.linspace(-1,0,100)
J=np.array([resonant_J(_S) for _S in S])
JE=J[:,0]
JB=J[:,1]
Ratio1=(JE/JB)

mj=np.array([MJ_spx(_S) for _S in S])
jspx=mj[:,0]
mspx=mj[:,1]
Ratio2=(jspx/mspx*S)
plt.plot(S[10:],Ratio1[10:])
plt.plot(S[10:],Ratio2[10:])
plt.twinx()
plt.plot(S[10:],(Ratio2-Ratio1)[10:])
import numpy as np
from numpy.linalg import norm

chor=np.load("../data/T_chor.npy")
sour=np.load("../data/T_sour.npy")
#E =-( chor[1:,...] - chor[:-1,...])
#J= (sour[1:,...]+sour[:-1,...])/2
E=chor
J=sour
#Je = J \dot E / abs(E)
Je = norm(J*E,axis=-1)/norm(E,axis=-1)
#Jb = sqrt( J^2 - Je^2 )
Jb=np.sqrt(norm(J,axis=-1)**2 - Je**2)
ratio=Je/Jb
