from scipy import integrate
from scipy import optimize
from chebpy import chebfun
import numpy as np
#define potential
def potential(S,zeta,wtr=1,k=1):
    return wtr**2/k**2 * (np.cos(zeta) + S *  zeta)

#find zeta is O point or X point
#energy O is greater then X
def min_potential_zeta(S,zeta1,zeta2):
    if(potential(S,zeta1) > potential(S,zeta2)):
        return zeta2, zeta1
    else:
        return zeta1, zeta2

def find_zeta_ox(S):
    #wtr can be any value
	#X point, energy first order derivate = 0; - sin z + S  = 0 
    zeta_1 = np.mod(np.arcsin(S),2*np.pi)
    zeta_2 = np.pi-np.arcsin(S)
    zeta_o,zeta_x = min_potential_zeta(S,zeta_1,zeta_2)
    return zeta_o,zeta_x

def find_zeta_ox2(S):
    #wtr can be any value
	#X point, energy first order derivate = 0; - sin z + S  = 0 
    zeta_1 = np.mod(np.arcsin(S),2*np.pi)
    zeta_2 = np.pi-np.arcsin(S)
    if (-np.cos(zeta_1)<0):
        zeta_x=zeta_1
        zeta_o=zeta_2
    else:
        zeta_x=zeta_2
        zeta_o=zeta_1
    return zeta_o,zeta_x


#take wtr S and return (zeta0,zeta1,zeta2) : lowest potential point, X point, and O point
def find_zeta(S):
    zeta_o,zeta_x = find_zeta_ox(S)
    const_E = potential(S,zeta_x)
    #use chebfun function to find the other intersection (zeta_o) of "E=const_E" line and "E = potential" curve
    sol=chebfun(lambda zeta: (np.cos(zeta) + S*zeta)-(const_E),[0,2*np.pi])
    #allzeta result may contains both zeta_o and zeta_x, choose the one with largest distance to our known zeta_x as zeta_o
    allzeta = []
    #filter nan
    for r in sol.roots():
        if np.isnan(r):
            continue
        else:
            allzeta.append(r)
    allzeta = np.array(allzeta)
    if (len(allzeta)!=1):
        d2zeta_x = np.abs(allzeta - zeta_x)
        zeta_c = allzeta[np.where(d2zeta_x==np.max(d2zeta_x))][-1]
    else:
        zeta_c=allzeta[0]
        zeta_c = np.array(zeta_c).flatten()[-1]
        zeta_x = np.array(zeta_x).flatten()[-1]
    return zeta_c,zeta_x,zeta_o

def hole_width(S,wtr,k):
    zeta_o,zeta_x = find_zeta_ox(S)
    const_E = potential(S,zeta_x,wtr,k)
    dO = np.sqrt(2*(const_E - potential(S,zeta_o,wtr,k)))/k
    return dO

##omura
def resonant_J(S):
    zeta_c,zeta_x,_ = find_zeta(S)
    #zeta_1 is the X point
    def currentB(zeta,zeta1,S):
        return (np.cos(zeta1) - np.cos(zeta) - S*(zeta - zeta1))**0.5 * np.cos(zeta)
    def currentE(zeta,zeta1,S):
        return (np.cos(zeta1) - np.cos(zeta) - S*(zeta - zeta1))**0.5 * np.sin(zeta)
    JE=-integrate.quad(currentE,zeta_x,zeta_c,args=(zeta_x,S))[0]
    JB=integrate.quad(currentB,zeta_x,zeta_c,args=(zeta_x,S))[0]
    return JE,JB

def MJ_spx(S):
    zeta_c,zeta_x,_ = find_zeta(S)
    #zeta_1 is the X point
    def Mspx(zeta,zeta1,S):
        return (np.cos(zeta1) + S*zeta1 - (np.cos(zeta) + S*zeta) )**0.5 * np.cos(zeta)
    def Jspx(zeta,zeta1,S):
        return (np.cos(zeta1) + S*zeta1 - (np.cos(zeta) + S*zeta) )**0.5
    mspx=integrate.quad(Mspx,zeta_x,zeta_c,args=(zeta_x,S))[0]
    jspx=integrate.quad(Jspx,zeta_x,zeta_c,args=(zeta_x,S))[0]
    return jspx,mspx 

def test_MJ_spx(S):
    zeta_c,zeta_x = find_zeta(S)
    #zeta_1 is the X point
    def Mspx(zeta,zeta1,S):
        return (np.cos(zeta1) + S*zeta1 - (np.cos(zeta) + S*zeta) )**0.5 * np.cos(zeta)
    def Jspx(zeta,zeta1,S):
        return (np.cos(zeta1) + S*zeta1 - (np.cos(zeta) + S*zeta) )**0.5
    mspx=integrate.quad(Mspx,zeta_c,zeta_x,args=(zeta_x,S))[0]
    jspx=integrate.quad(Jspx,zeta_x,zeta_c,args=(zeta_x,S))[0]
    return S*jspx,mspx 


def ratioJ_w(S):
    #JB/JE
    jspx,mspx = MJ_spx(S)
    return mspx/jspx/S 

def ratioJ_o(S):
    Je,Jb = resonant_J(S) 
    return Jb/Je


#the nonlinear equation is
#jspx/mspx(S) * S - const = 0
#this constant is current ratio

#return function: take ratio give R
def interpR():
	ANS_R = np.logspace(-4,-0.0000001,1000)
	ANS_RATIO_J = np.array([ratioJ_w(r) for r in ANS_R])
	from scipy import interpolate
	calcR_interp = interpolate.interp1d(abs(ANS_RATIO_J), ANS_R)
	return calcR_interp

def calcR_w(ratio):
    sol = []
    for r in ratio:
        if(abs(r) > 0.5):
            ig = 0.2
        else:
            ig = 0.99
        print(r,ig)
        _sol = optimize.root(lambda x:ratioJ_w(x) - r, ig)
        sol.append(_sol.x[0])
    return np.array(sol)

def calcR_o(ratio):
    sol = []
    for r in ratio:
        if(abs(r) > 0.5):
            ig = 0.2
        else:
            ig = 0.99
        print(r,ig)
        _sol = optimize.root(lambda x:ratioJ_o(x) - r, ig)
        sol.append(_sol.x[0])
    return np.array(sol)

#20230501
#calculate by definition
def f_vperp(gyro,J,pcr,Omega=0):
    return np.sqrt(2*gyro*(J+pcr+Omega))


def f_Jact(omega,gyro,kmode,gyro0,vperp,kmode0):
    pcr0 = (omega - gyro0)/kmode0**2
    Jact0 = vperp*vperp/gyro0 - pcr0
    const=(omega-gyro0)*Jact0+0.5*(omega-gyro0)**2/kmode0**2
    return (const-0.5*(omega-gyro)**2/kmode**2)/(omega-gyro)

def f_mu(Jact,gyro,freq,k):
    mu = Jact + (freq - gyro)/k/k
    return mu

def f_wb2(gyro,J,a,freq,k):
    #mu=J+Omega+PI = J+ p_para/k
    mu = J + (freq - gyro)/k/k
    vperp = np.sqrt(2*gyro*mu)
    wtr2= vperp * k**2 * a
    return wtr2
def f_do(dw,dk,vr,k0):
    dO = (dw - vr *dk)/k0**2
    return dO

#J is a time constant

#linear dispersion + k by dispersion
def f_alpha_t(freq,freq_t,gyro,gyro_s,wb2,k,vr,vg,J):
    mu = J + (freq - gyro)/k/k
	#linear vg and vr
    alp = 1/wb2*(1-vr/vg)**2 * freq_t + (k*mu - 3*vr/2)/wb2 * gyro_s
    return alp
#nonlinear dispersion + k by dis
    #nonlinear
def f_alpha_w(freq_t,gyro_s,wb2,k,kl,J,PI):
    #nonlinear dp: vg = - vr
    alp = 4/wb2 * freq_t + (J + (kl/2/k - 1)*PI)*k/wb2 * gyro_s
    return alp
#all by defination
def f_alpha(freq_t,gyro_s,wb2,k,k_s,vr,vg,J):
	# vg and vr imsimulation
    alp = 1/wb2 *(1-2*vr/vg) * freq_t + vr**2*k_s *(J)*k/wb2 * gyro_s
    return alp

#chirping rate at source
def c_rate_s(freq,gyro,gyro_s,k,vr,vg,J):
    mu = J + (freq - gyro)/k/k
    c_rate = (1-vr/vg)**-2 * (k*mu - 3*vr/2) * gyro_s
    return c_rate

def qk_wb(B,mu,w,vp):
    return np.sqrt(mu*w*vp*B)

def f_feq(vperp,vpara,gyro0,gyro,omega,kl,Jact,Omega,loss_cone=1e-16):
    small=1e-10
    PI=(omega-gyro)/kl/kl
    mu = Jact+PI+Omega
    mu[mu<0] = small
    feq = gyro0/((2.0*np.pi)**1.5*vperp*vperp*vpara) *np.exp(-mu*gyro0/vperp**2)*np.exp(-mu*(gyro-gyro0)/vpara**2)*np.exp(-0.5*(kl*(PI+Omega)/vpara)**2)
    dfeq =feq * (-gyro0/vperp**2 - (gyro-gyro0)/vpara**2 - kl**2*(PI+Omega)/vpara**2)
    return (feq,dfeq)
