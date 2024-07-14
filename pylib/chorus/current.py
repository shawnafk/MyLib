import numpy as np
#djc/dt = wpe^2/4/pi * E + Jc cross wce
#divide by 4pi due to the difference of wpe definition bwtween cgs and SI unit
def ColdCurrentFromFluid(Ex,Ey,dT,wpe,wce):
    tsteps,sgrids=Ex.shape
    Jx = np.zeros((tsteps,sgrids))
    Jy = np.zeros((tsteps,sgrids))
    for t in range(tsteps-1):
        Ax = -Jx[t,:]*(wce*dT/2 - 2/wce/dT) + 2 * Jy[t,:]+ wpe**2/4/np.pi * (2/wce*Ex[t,:] + dT*Ey[t,:])
        Ay = -Jy[t,:]*(wce*dT/2 - 2/wce/dT) - 2 * Jx[t,:]+ wpe**2/4/np.pi * (2/wce*Ey[t,:] - dT*Ex[t,:])
        #Ax = -Jx[t,:]*(wce*dT/2 - 2/wce/dT) + 2 * Jy[t,:]+ wpe**2 * (2/wce*Ex[t,:] + dT*Ey[t,:])
        #Ay = -Jy[t,:]*(wce*dT/2 - 2/wce/dT) - 2 * Jx[t,:]+ wpe**2 * (2/wce*Ey[t,:] - dT*Ex[t,:])
        a = wce*dT/2 + 2/wce/dT
        Jx[t+1,:] = Ax/a
        Jy[t+1,:] = Ay/a
    return Jx,Jy

#cgs unit, j = -i w / 4/pi * chi *E
#chi = -wpe^2/w/(w-wce) : for electron
def ColdCurrentFromDP(Ex,Ey,w,wpe,wce):
    chi = wpe**2/w/(wce-w)
    J = -1j*w/4/np.pi*chi*(Ex+1j*Ey)
    return np.real(J),np.imag(J)

def fluid_cgs(Ex,Ey,dT,wpe,wce):
    tsteps,sgrids=Ex.shape
    Jx = np.zeros((tsteps,sgrids))
    Jy = np.zeros((tsteps,sgrids))
    for t in range(tsteps-1):
        Ax = -Jx[t,:]*(wce*dT/2 - 2/wce/dT) + 2 * Jy[t,:]+ wpe**2/4/np.pi * (2/wce*Ex[t,:] + dT*Ey[t,:])
        Ay = -Jy[t,:]*(wce*dT/2 - 2/wce/dT) - 2 * Jx[t,:]+ wpe**2/4/np.pi * (2/wce*Ey[t,:] - dT*Ex[t,:])
        #Ax = -Jx[t,:]*(wce*dT/2 - 2/wce/dT) + 2 * Jy[t,:]+ wpe**2 * (2/wce*Ex[t,:] + dT*Ey[t,:])
        #Ay = -Jy[t,:]*(wce*dT/2 - 2/wce/dT) - 2 * Jx[t,:]+ wpe**2 * (2/wce*Ey[t,:] - dT*Ex[t,:])
        a = wce*dT/2 + 2/wce/dT
        Jx[t+1,:] = Ax/a
        Jy[t+1,:] = Ay/a
    return Jx,Jy

def fluid_si(Ex,Ey,dT,wpe,wce):
    tsteps,sgrids=Ex.shape
    Jx = np.zeros((tsteps,sgrids))
    Jy = np.zeros((tsteps,sgrids))
    for t in range(tsteps-1):
        Ax = -Jx[t,:]*(wce*dT/2 - 2/wce/dT) + 2 * Jy[t,:]+ wpe**2 * (2/wce*Ex[t,:] + dT*Ey[t,:])
        Ay = -Jy[t,:]*(wce*dT/2 - 2/wce/dT) - 2 * Jx[t,:]+ wpe**2 * (2/wce*Ey[t,:] - dT*Ex[t,:])
        #Ax = -Jx[t,:]*(wce*dT/2 - 2/wce/dT) + 2 * Jy[t,:]+ wpe**2 * (2/wce*Ex[t,:] + dT*Ey[t,:])
        #Ay = -Jy[t,:]*(wce*dT/2 - 2/wce/dT) - 2 * Jx[t,:]+ wpe**2 * (2/wce*Ey[t,:] - dT*Ex[t,:])
        a = wce*dT/2 + 2/wce/dT
        Jx[t+1,:] = Ax/a
        Jy[t+1,:] = Ay/a
    return Jx,Jy
