import numpy as np
from fig_head import *

def density(R,Z,n0,R0):
    r=np.sqrt(R**2+Z**2)
    coslmd = (R/r)
    r0=r/coslmd**2
    dst_right = n0*(r0/R0)**-4
    sig=R0/8
    dst_left = n0*np.exp(-(r0-R0)**2/2/sig**2)
    mask_right=np.zeros(R.shape)
    mask_right[r0>R0] = 1
    mask_left=np.zeros(R.shape)
    mask_left[r0<R0] = 1
    dst = (dst_left * mask_left + dst_right * mask_right)*coslmd
    return dst

_r=np.linspace(0.05,1,200)
_z=np.linspace(-0.3,0.3,100)

R,Z=np.meshgrid(_r,_z)

dst = density(R,Z,1,0.55)
fig,ax=plt.subplots()
ax.pcolormesh(R,Z,dst,shading='gouraud')
ax1=ax.twinx()
ax1.plot(_r,dst[50,:],'r')
plt.show()
