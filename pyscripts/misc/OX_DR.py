import numpy as np
import PlasmaDR as pdr
from sympic_io.tools import unit
import matplotlib.pyplot as plt

L=1e-3
U=unit.gen_unit(L)
case=pdr.Plasma(1,1,U)



Ow,Xw=case.oxdr(1e19,2,np.logspace(-2,1,1000))
Lw,Rw=case.lrdr(1e19,2,np.logspace(-2,1,1000))
plt.scatter(Ow[0],Ow[1],s=0.2,label='O')
plt.scatter(Xw[0],Xw[1],s=0.2,label='X')
plt.scatter(Lw[0],Lw[1],s=0.2,label='L')
plt.scatter(Rw[0],Rw[1],s=0.2,label='R')
plt.legend()
plt.show()
