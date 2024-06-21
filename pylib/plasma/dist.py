import numpy as np
#std norm dist density
def std_norm_den(sig,vmax=1,vmid=0,N=1024):
    x= np.linspace(-vmax,vmax,N)
    f = np.exp(-(x-vmid)**2/2/sig**2)/sig/np.sqrt(2*np.pi)
    return x,f


def perp(x,y):
    return np.sqrt(x**2+y**2)
