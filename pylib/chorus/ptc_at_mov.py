import numpy as np
from scipy.integrate import solve_ivp
from chebpy import chebfun

#functions related to the Hamiltonian
#H
def Kamitonian(xi,Omg,args):
    J,Pi,gyro,dgyro,kl,a = args
    mu = J + Omg + Pi
    drg = 1/kl*(J-Pi/2)*dgyro
    T = kl**2/2*Omg**2 
    V = np.sqrt(2*gyro*mu) * np.real(a * np.exp(-1j*xi)) + drg * xi
    return T+V 

#phase space flows
def phaseflow(t,cor,*args):
    xi=cor[0]
    Omg = cor[1]
    J,Pi,gyro,dgyro,kl,a = args
    mu = J + Omg + Pi
    dxi = kl**2 * Omg + gyro/np.sqrt(2*gyro*mu) * np.real(a * np.exp(-1j*xi))
    drg = 1/kl*(J-Pi/2)*dgyro
    ref= np.sqrt(2*gyro*mu) * np.imag(a*np.exp(-1j*xi))
    dO = - ref - drg
    return dxi,dO

#solves dynamics *arg take all below as a list, ** args take all as dict?
def interpY(znew,idxL,zpos,Y):
    return np.interp(znew, (zpos[idxL],zpos[idxL+1]), (Y[idxL],Y[idxL+1]))

def ivp_ptc(iv,s0,t0,iters,para_list,subT=48.9):
    z,J,Pi,b,db,k,a,vr = para_list
    ivs=np.zeros((iters+1,4))
    ivs[0,:2] =  iv
    ivs[0,2:] =  [z[s0],vr[s0]]
    _J,_Pi,_b,_db,_k,_a,_vr = (J[s0],Pi[s0],b[s0],db[s0],k[s0],a[t0,s0],vr[s0])
    for i in range(iters):
        t = t0 + i
        #sol = solve_ivp(phaseflow,[0,subT],ivs[i,:2],args=[J[s],Pi[s],b[s],db[s],k[s],a[t,s]])
        sol = solve_ivp(phaseflow,[0,subT],ivs[i,:2],args=[_J,_Pi,_b,_db,_k,_a])
        #calculate vnew, use last Omega
        vnew = sol.y[1][-1]*_k + _vr
        #calculate snew = s0 + vnew * dT
        snew = ivs[i,2] + vnew * subT
        #interp new parameter
        idxL = np.searchsorted(z,snew)
        _J=interpY(snew,idxL,z,J)
        _Pi=interpY(snew,idxL,z,Pi)
        _b=interpY(snew,idxL,z,b)
        _db=interpY(snew,idxL,z,db)
        _k=interpY(snew,idxL,z,k)
        _vr=interpY(snew,idxL,z,vr)
        _a=interpY(snew,idxL,z,a[t,:])
        #ivs[i+1,:] = [np.mod(sol.y[0][-1],2*np.pi),sol.y[1][-1]]
        ivs[i+1,:2] = [sol.y[0][-1],sol.y[1][-1]]
        ivs[i+1,2:] = [snew,vnew]
    return ivs

#----------- iter with in one cell
def ivp_single(ivs,s0,t0,para_list,subT=48.9,N=100):
    z,J,Pi,b,db,k,a = para_list
    singleT=subT*N
    dT=5
    NT = int(singleT/dT)
    tsave=np.linspace(0,singleT,NT)
    Ncell = ivs.shape[0]-1
    ptc=np.zeros((Ncell,NT,2))
    _J,_Pi,_b,_db,_k,_a = (J[s0],Pi[s0],b[s0],db[s0],k[s0],a[t0,s0])
    for i in np.arange(Ncell):
        t=t0 + i
        sol = solve_ivp(phaseflow,[0,singleT],ivs[i,:2],args=[_J,_Pi,_b,_db,_k,_a],t_eval=tsave)
        snew = ivs[i,2]
        idxL = np.searchsorted(z,snew)
        _J=interpY(snew,idxL,z,J)
        _Pi=interpY(snew,idxL,z,Pi)
        _b=interpY(snew,idxL,z,b)
        _db=interpY(snew,idxL,z,db)
        _k=interpY(snew,idxL,z,k)
        _a=interpY(snew,idxL,z,a[t,:])
        ptc[i,:,0]=sol.y[0]
        ptc[i,:,1]=sol.y[1]
    return ptc

#phase space derivatives dOmega/dxi
def dO_dxi(xi,Omg,*args):
    J,Pi,gyro,dgyro,kl,a = args
    mu = J + Omg + Pi
    dxi = kl**2 * Omg + gyro/np.sqrt(2*gyro*mu) * np.real(a * np.exp(-1j*xi))
    drg = 1/kl*(J-Pi/2)*dgyro
    ref= np.sqrt(2*gyro*mu) * np.imag(a*np.exp(-1j*xi))
    dO = - ref - drg
    return dO/dxi


# --- ----------  integrate procedure
#H as function of xi, with fixed parameter and Omega = 0 (approximation) 
def Kam_xi(K0,*args):
    J,Pi,gyro,dgyro,kl,a = args
    mu = J +  Pi
    drg = 1/kl*(J-Pi/2)*dgyro
    return lambda xi: np.sqrt(2*gyro*mu) * np.real(a * np.exp(-1j*xi)) + drg * xi - K0



#given initial condtion and paramters, calculate the area of trajecory
def get_area(iv,ptc_xi,_J,_Pi,_b,_db,_k,_a):
    #get H0
    K0 = Kamitonian(iv[0],iv[1],[_J,_Pi,_b,_db,_k,_a])
    #solve the boundary, two end points of xi, Omega approximately 0
    f_k_xi =  Kam_xi(K0,_J,_Pi,_b,_db,_k,_a)
    #solroot = chebfun(f_k_xi,[0,2*np.pi])
    solroot = chebfun(f_k_xi,[np.min(ptc_xi)-0.1,np.max(ptc_xi)+0.1])
    if len(solroot.roots()) > 1:
        xil=solroot.roots()[0]
        xim=solroot.roots()[1]
        ts=np.linspace(xil,xim,1000)
        boundary=solve_ivp(dO_dxi,[xil,xim],[0],args=[_J,_Pi,_b,_db,_k,_a],t_eval=ts)
        I = np.abs(np.trapz(boundary.y[0],ts)*2)
    else:
         I = 0
    return I

def get_I(ivs,ptc_xi,s0,t0,para_list):
    z,J,Pi,b,db,k,a = para_list
    Nc=ivs.shape[0]-1
    I=np.zeros(Nc)
    _J,_Pi,_b,_db,_k,_a = (J[s0],Pi[s0],b[s0],db[s0],k[s0],a[t0,s0])
    for i in range(Nc):
        #print(i)
        t=t0 + i
        I[i] = get_area(ivs[i,:],ptc_xi[i,:],_J,_Pi,_b,_db,_k,_a)
        snew = ivs[i,2]
        idxL = np.searchsorted(z,snew)
        _J=interpY(snew,idxL,z,J)
        _Pi=interpY(snew,idxL,z,Pi)
        _b=interpY(snew,idxL,z,b)
        _db=interpY(snew,idxL,z,db)
        _k=interpY(snew,idxL,z,k)
        _a=interpY(snew,idxL,z,a[t,:])
    return I

def get_E(ivs,s0,t0,para_list):
    z,J,Pi,b,db,k,a = para_list
    Nc = ivs.shape[0]-1
    E0 = np.zeros(Nc)
    _J,_Pi,_b,_db,_k,_a = (J[s0],Pi[s0],b[s0],db[s0],k[s0],a[t0,s0])
    for i in range(Nc):
        t = t0 + i
        E0[i] = Kamitonian(ivs[i,0],ivs[i,1],[_J,_Pi,_b,_db,_k,_a])
        snew = ivs[i,2]
        idxL = np.searchsorted(z,snew)
        _J=interpY(snew,idxL,z,J)
        _Pi=interpY(snew,idxL,z,Pi)
        _b=interpY(snew,idxL,z,b)
        _db=interpY(snew,idxL,z,db)
        _k=interpY(snew,idxL,z,k)
        _a=interpY(snew,idxL,z,a[t,:])
    return E0

def local_iter(iv,s0,t0,iters,para_list,subT=48.9):
    J,Pi,b,db,k = para_list
    ivs=np.zeros((iters+1,2))
    ivs[0,:2] =  iv
    _J,_Pi,_b,_db,_k = (J[s0],Pi[s0],b[s0],db[s0],k[s0])
    for i in range(iters):
        a = get_a(i*subT)
        sol = solve_ivp(phaseflow,[0,subT],ivs[i,:2],args=[_J,_Pi,_b,_db,_k,a])
        ivs[i+1,:] = [sol.y[0][-1],sol.y[1][-1]]
    return ivs


