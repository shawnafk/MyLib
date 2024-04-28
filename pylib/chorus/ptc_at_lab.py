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

#in one DT (48.9) which is much smaller than the period, a is static
#a is changing while moving along field line
#for static location, 
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
def ivp_ptc(iv,s0,t0,iters,para_list,subT=48.9):
    J,Pi,b,db,k,a = para_list
    ivs=np.zeros((iters+1,2))
    ivs[0,:] =  iv
    for i in range(iters):
        s = s0 - i
        t = t0 + i
        sol = solve_ivp(phaseflow,[0,subT],ivs[i,:],args=[J[s],Pi[s],b[s],db[s],k[s],a[t,s]])
        #ivs[i+1,:] = [np.mod(sol.y[0][-1],2*np.pi),sol.y[1][-1]]
        ivs[i+1,:] = [sol.y[0][-1],sol.y[1][-1]]
    return ivs

#----------- iter with in one cell
def ivp_single(ivs,s0,t0,para_list,subT=48.9,N=100):
    J,Pi,b,db,k,a = para_list
    singleT=subT*N
    dT=5
    NT = int(singleT/dT)
    tsave=np.linspace(0,singleT,NT)
    Ncell = ivs.shape[0]-1
    ptc=np.zeros((Ncell,NT,2))
    for i in np.arange(Ncell):
        s=s0 - i
        t=t0 + i
        sol = solve_ivp(phaseflow,[0,singleT],ivs[i,:],args=[J[s],Pi[s],b[s],db[s],k[s],a[t,s]],t_eval=tsave)
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
    J,Pi,b,db,k,a = para_list
    Nc=ivs.shape[0]-1
    I=np.zeros(Nc)
    for i in range(Nc):
        #print(i)
        s=s0 - i
        t=t0 + i
        I[i] = get_area(ivs[i,:],ptc_xi[i,:],J[s],Pi[s],b[s],db[s],k[s],a[t,s])
    return I

