import os
import numpy as np
from scipy import optimize
from scipy import integrate
from chebpy import chebfun
import matplotlib.pyplot as plt

#define potential
def potential(wtr,S,zeta):
    return 2* wtr**2 * (np.cos(zeta) - S *  zeta)

#find zeta is O point or X point
def min_potential_zeta(wtr,S,zeta1,zeta2):
    if(potential(wtr,S,zeta1) > potential(wtr,S,zeta2)):
        return zeta2
    else:
        return zeta1

#take wtr S and return (zeta0,zeta1,zeta2) : lowest potential point, X point, and O point
def find_zeta(S):
    #wtr can be any value
    wtr=1
    zeta_x=np.mod(-np.arcsin(S),2*np.pi)
    const_E = potential(wtr,S,zeta_x) 
    #use chebfun function to find the other intersection (zeta_o) of "E=const_E" line and "E = potential" curve
    sol=chebfun(lambda zeta:2*wtr**2*(np.cos(zeta) - S*zeta)-(const_E),[0,2*np.pi])
    allzeta = sol.roots()
    #allzeta result may contains both zeta_o and zeta_x, choose the one with largest distance to our known zeta_x as zeta_o 
    if (len(allzeta)!=1):
        d2zeta_x = np.abs(allzeta - zeta_x)
        zeta_o = allzeta[np.where(d2zeta_x==np.max(d2zeta_x))][-1]
    else:
        zeta_o=allzeta[0]
        zeta_o = np.array(zeta_o).flatten()[-1]
        zeta_x = np.array(zeta_x).flatten()[-1]
    return np.array((zeta_o,zeta_x))

##omura
def resonant_J(S):
    zeta_o,zeta_x = find_zeta(S)
    #zeta_1 is the X point
    def currentB(zeta,zeta1,S):
        return (np.cos(zeta1) - np.cos(zeta) + S*(zeta - zeta1))**0.5 * np.cos(zeta)
    def currentE(zeta,zeta1,S):
        return (np.cos(zeta1) - np.cos(zeta) + S*(zeta - zeta1))**0.5 * np.sin(zeta)
    JE=-integrate.quad(currentE,zeta_x,zeta_o,args=(zeta_x,S))[0]
    JB=integrate.quad(currentB,zeta_x,zeta_o,args=(zeta_x,S))[0]
    return JE,JB

def MJ_spx(S):
    zeta_o,zeta_x = find_zeta(S)
    #zeta_1 is the X point
    def Mspx(zeta,zeta1,S):
        return (np.cos(zeta1) - S*zeta1 - (np.cos(zeta) - S*zeta) )**0.5 * np.cos(zeta)
    def Jspx(zeta,zeta1,S):
        return (np.cos(zeta1) - S*zeta1 - (np.cos(zeta) - S*zeta) )**0.5
    mspx=integrate.quad(Mspx,zeta_x,zeta_o,args=(zeta_x,S))[0]
    jspx=integrate.quad(Jspx,zeta_x,zeta_o,args=(zeta_x,S))[0]
    return jspx,mspx 

def ratioJ_w(S):
    #JB/JE
    jspx,mspx = MJ_spx(S)
    return mspx/jspx/S 

def ratioJ_o(S):
    Je,Jb = resonant_J(S) 
    return Jb/Je

def pull_phase(compA,dim=0):
    return np.unwrap(np.angle(compA),axis=dim)

#the nonlinear equation is
#jspx/mspx(S) * S - const = 0
#this constant is current ratio

#saved_f = "./tmpRJ.npy"
#ANS_R = np.logspace(-4,-0.0000001,1000)
#if(os.path.exists(saved_f)):
#    ANS_RATIO_J = np.load(saved_f)
#else:
#    ANS_RATIO_J = np.array([ratioJ_w(r) for r in ANS_R])
#    np.save(saved_f,ANS_RATIO_J)
#from scipy import interpolate
#calcR_interp = interpolate.interp1d(abs(ANS_RATIO_J), ANS_R)

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


def calcWK(wl,wc,wp,kl,phi_w,dT,zpos):
    #phase unwrap along dim 0; make continueous w.r.t time
    #omega = \partial phi/ \partial t
    phi_w_t = np.unwrap(phi_w[:,:],axis=0)

    #phi_w_t = phi_w[:,:]
    dphi_t = (phi_w_t[1:,:] -  phi_w_t[:-1,:])
    wsim = dphi_t/dT

    #phase unwrap along dim 1; make continueous w.r.t space
    #k = -\partial phi/ \partial s
    #phi_w_s = phi_w[:,:]
    phi_w_s = np.unwrap(phi_w[:,:],axis=1)
    dphi_z = (phi_w_s[:,1:] -  phi_w_s[:,:-1])/1
    dz = (zpos[np.newaxis,1:] - zpos[np.newaxis,:-1])
    ksim=-(dphi_z[:,:-1]/dz)

    wall=wsim+wl
    klin = np.sqrt(wall[:,:-1]**2+wp**2*wall[:,:-1]/(wc-wall[:,:-1]))

    #kl2 = kmode, gap comes from freq
    #left aligned
    knl = (wc-wl)/(wc - wall[:,:-1])*kl[np.newaxis,:]
    kall = kl[:-1]+ ksim[:-1,:]
    return wall,klin,knl,kall

def calc_ko(wall,wc,wp,gm):
    R = gm/wall
    ko_p = np.sqrt(wall**2*(1 + 2*R) + wp**2*wall/(wc-wall)*(1 + wc/(wc-wall) * R) )
    ko_m = np.sqrt(wall**2*(1 - 2*R) + wp**2*wall/(wc-wall)*(1 - wc/(wc-wall) * R) )

    return ko_p,ko_m


def calcVG(wce,wp2,wall,ksim,zpos,dT):
    dz=zpos[1:] - zpos[:-1]
    dfdz=(wall[:,1:-1] - wall[:,:-2])/dz[np.newaxis,:]
    dfdt=(wall[1:,:] - wall[:-1,:])/dT
    vgsim =-dfdt[:,:-2]/dfdz[1:,:]
    vg0sim = 2*ksim/(2*wall[:,:-2] + wce[:-1]*wp2[:-1]/(wce[:-1]-wall[:,:-2])**2)
    vrsim = (wall[:,1:-1] - wce[:-1]) /ksim
    return vg0sim,vrsim,vgsim

def calcAdv(amp,zpos,vg,dT):
    da_t = amp[1:,:] -  amp[:-1,:]
    dadt = da_t / dT
    da_s = amp[:,1:] -  amp[:,:-1]
    dz = zpos[np.newaxis,1:] - zpos[np.newaxis,:-1]
    dads=da_s[:,1:]/dz
    adv = dadt[:,1:-1] + dads[:-1,:]*vg[np.newaxis,:-1]
    return adv
def Lshell(xi,a0=4.5,B0=0.312,RE=6378e3):
    c=2.9997e8
    omega_ce_surf = 1.76e7*B0
    return 1/(a0/(xi*c**-2*omega_ce_surf**2*RE**2))**(1/4)

def cvg(w,k,wce,wpe):
    c=1
    vg = 2*c**2*k/(2*w + wce*wpe**2/(wce-w)**2)
    return vg

def cvr(w,k,wce,gamma):
    return (w-wce/gamma)/k

def gamma(v):
    return 1/np.sqrt(1-v**2)

def D(w,k,wpe,wce):
    return k**2 - w**2 -w*wpe**2/(wce-w)

def time_ratio(R,w,k,wce,wpe,gamma):
    c=1
    vg = cvg(w,k,wce,wpe)
    vr = cvr(w,k,wce,gamma)
    return 2*np.pi*R/gamma/(1-vr/vg)**2

#checked by our results
def GM_L(wl,nh,v_para,v_perp,wpe,wce,wpe0,wce0,beta):
    aniso = (v_perp/v_para)**2
    kl=linear_k(wl,wpe,wce)
    vg=cvg(wl,kl,wce,wpe)
    gml = np.sqrt(2*np.pi) * wce*vg*nh*wpe0**2/(4*kl**2*v_para)*np.exp(-(wl-wce)**2/(2*kl**2*v_para**2)) * (1+aniso*(wce-wce0)/wce0)**-2 *  (1+beta*aniso*(wce-wce0)/wce0)**-2 * (beta*aniso**2*((wce+wce0-2*wl)*(wce-wce0))/wce0**2 + (1+beta)*aniso*(wce0-wl)/wce0 -1)
    return gml

def linear_k(w,wpe,wce):
    return  np.sqrt(w**2 + w*wpe**2/(wce-w))

#checked by Tao 2017
def OMG_Tr(w,wpe,wce,v_perp,B1,delta,lorentz_f):
    #calculate k v_1 in unit of gyro frequency
    #k from linear dispersion relation, given omega in Omega_0 return k in Omega_0/c unit
    #v_1 in c
    kl=linear_k(w,wpe,wce)
    return np.sqrt(kl*v_perp*B1*wce)*delta*lorentz_f**-0.5

#our way to calculate trapping frequency
#J is a constant load from file, O is approximate 0, Pcr is constant
#in unit wpe only
def OMG_Tr2(wl,wpe,wce,Jact,a):
    kl=linear_k(wl,wpe,wce)
    pcr= (wl - wce)/kl/kl
    return np.sqrt(np.sqrt(2*wce*(Jact + pcr)) * kl**2 * a)

def fit_gm(B):
    a1=0.0371
    a2=0.76
    return (B/a1)**(1/a2)

def readchor(file,length,dim):
    chor = np.loadtxt(file)
    if dim>1:
        newc = chor.reshape([int(chor.shape[0]/length),length,dim])
    else:
        newc = chor.reshape([int(chor.shape[0]/length),length])
    return newc

## calculate basic data on a magnetic field line of Earth's dipole field.
def zpos():
    return 0
def wce(zpos):
    return 0
def wpe(zpos):
    return 0

#def ord2(zpos,z0,dkdt,vr,vg):
#    return 0.5*(-1/vr-1/vg)*dkdt*(zpos-z0)**2
#exp2nd = ord2(zpos,zpos[512],dkdt[:,:],vr,vg)
