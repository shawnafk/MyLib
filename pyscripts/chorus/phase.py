import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

wtr = 0.5;S =  -0.4
#(zeta0,zeta1,zeta2) : O point X point and the phase of the separatrix crossing theta = 0
zeta0 = np.pi+np.arcsin(S)
zeta1 = -np.arcsin(S) 
X = 2* wtr**2 * (np.cos(zeta1) - S *  zeta1)
#solve zeta2, initial value 20 is choosen far away from zeta1, should be different for different S and wtr
sol=optimize.root(lambda zeta:2*wtr**2*(np.cos(zeta)-S*zeta)-X,20)
zeta2 = sol.x[-1]


fig,(ax1,ax2)=plt.subplots(2,1)
plt.subplots_adjust(hspace=0)
zeta = np.linspace(0,2*np.pi,100)
ax1.plot(zeta,2*wtr**2*(np.cos(zeta)-S*zeta))
ax1.set_xticks([])
#phase contour
theta = np.linspace(-1,1,50) * wtr * 2 
Z,T = np.meshgrid(zeta,theta)
phase = T**2 + 2* wtr**2 * (np.cos(Z) - S *  Z)
ax2.contour(Z,T,phase,50)

#characteristic points
ax2.scatter(zeta0,0,s=100,marker='o',color='w',edgecolor='g')
ax2.scatter(zeta1,0,s=100,marker='x',lw=1.5)
ax2.scatter(zeta2,0,s=100,marker='o',lw=1.5)

#upper and lower theta boundary of the trapping region
zetab = np.linspace(zeta1,zeta2,100)
theta_b1 = wtr*np.sqrt(2*(np.cos(zeta1) - np.cos(zetab) + S*(zetab-zeta1)) )
theta_b2 = -theta_b1
ax2.plot(zetab,theta_b1,'r')
ax2.plot(zetab,theta_b2,'k')
plt.show()
