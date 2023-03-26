from sympy import symbols,diff,pi,lambdify,sqrt,exp
import numpy as np
from scipy.optimize import fsolve

def funcG(a,b,c,d):
    wl,wce0,vpar,vper,nh,beta = symbols('wl,wce0,vpar,vper,nh,beta')
    kl = sqrt(wl**2 + wl/(wce0 - wl))
    vg = 2*kl/(2*wl + wce0/(wce0-wl)**2)
    gamma = sqrt(2*pi)*wce0*vg*nh/(4*kl**2*vpar)*exp(-0.5*(wl-wce0)**2/kl**2/vpar**2)*((1+beta)*(vper/vpar)**2*(wce0-wl)/wce0 -1)
    dgamma = diff(gamma,wl)
    ddgamma = diff(dgamma,wl)
    fgamma = lambdify([wl,wce0],gamma.subs({vpar:a,vper:b,nh:c,beta:d}))
    fdgamma = lambdify([wl,wce0],dgamma.subs({vpar:a,vper:b,nh:c,beta:d}))
    fddgamma = lambdify([wl,wce0],ddgamma.subs({vpar:a,vper:b,nh:c,beta:d}))
    return fgamma,fdgamma,fddgamma

fgamma,fdgamma,fddgamma = funcG(0.15,0.3,0.002,0.3)
from scipy.optimize import root
N=100
fwce0=np.linspace(0.2,0.6,N)    
fomega=np.zeros(N)
def solv_w(f,jac,oce,ig):
    g = lambda x: f(x, oce)
    dg = lambda x: jac(x, oce)
    sol = root(g, ig, jac=dg, method= 'hybr')
    return sol.x
 
i=0
ig=0.05
for gyro0 in fwce0:
    fomega[i] = solv_w(fdgamma,fddgamma,gyro0,ig)
    ig=[fomega[i]]
    i=i+1
print(fomega)

wce = np.linspace(0.1,0.6,1000)
wl = np.linspace(0.01,0.4,800)
wce_r,wl_c = np.meshgrid(wce,wl)
G = fgamma(wl_c,wce_r)
