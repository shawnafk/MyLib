import numpy as np
#f_ definition

#background profile
#lbd=(z/(Lshell*RE))
def f_gyro(lbd,R_a):
    fout = (1.0 + R_a*lbd**2)
    return fout

#density aka wpe**2
def f_wp2(lbd,R_b):
    fout = 1.0 + R_b*lbd**2
    return fout

def f_k(w,wce,wpe):
    c=1
    k=np.sqrt(w**2+w*wpe**2/(wce-w))/c
    return k

def f_w(k,wce,wpe):
    c=1
    w=wce/(wpe**2/c**2/k**2+1)
    return w


#knl = (wc-wl)/(wc - wall)*kl

def f_vg(w,wce,wpe):
    c=1
    k=f_k(w,wce,wpe)
    vg = 2*c**2*k/(2*w + wce*wpe**2/(wce-w)**2)
    return vg

def f_vr(w,wce,wpe):
    k=f_k(w,wce,wpe)
    vr = (w - wce) /k
    return vr

def f_Omega(w,wl,wce,wpe):
    vr_0 = f_vr(wl,wce,wpe)
    k_0 = f_k(wl,wce,wpe)
    PI = vr_0/k_0 
    vr = f_vr(w,wce,wpe)
    k = f_k(w,wce,wpe)
    return vr/k - PI

def calcWK(phi_w,dT,dz):
    #axis 0: t, 1: z
    #phase unwrap along dim 0; make continueous w.r.t time
    #omega = \partial phi/ \partial t
    phi_w_t = np.unwrap(phi_w[:,:],axis=0)
    #phi_w_t = phi_w[:,:]
    dphi_t = np.gradient(phi_w_t,axis=0)
    deltaw = dphi_t/dT
    #phase unwrap along dim 1; make continueous w.r.t space
    #k = -\partial phi/ \partial s
    #phi_w_s = phi_w[:,:]
    phi_w_s = np.unwrap(phi_w[:,:],axis=1)
    #dphi_z = (phi_w_s[:,1:] -  phi_w_s[:,:-1])/1
    dphi_z = np.gradient(phi_w_s,axis=1)
    deltak=-(dphi_z/dz)
    return deltaw,deltak

def calcVg(wall,kall,dT,dz):
    dfdt = np.gradient(wall,axis=0)/dT
    dfdz= np.gradient(wall,axis=1)/dz
    vgconsw =-dfdt/dfdz
    return vgconsw

#----------------------------------------- above are revised --------
def calc_ko(wall,wc,wp,gm):
    R = gm/wall
    ko_p = np.sqrt(wall**2*(1 + 2*R) + wp**2*wall/(wc-wall)*(1 + wc/(wc-wall) * R) )
    ko_m = np.sqrt(wall**2*(1 - 2*R) + wp**2*wall/(wc-wall)*(1 - wc/(wc-wall) * R) )

    return ko_p,ko_m


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
    vg=cvg(wl,kl,wce,wpe) # type: ignore
    gml = np.sqrt(2*np.pi) * wce*vg*nh*wpe0**2/(4*kl**2*v_para)*np.exp(-(wl-wce)**2/(2*kl**2*v_para**2)) * (1+aniso*(wce-wce0)/wce0)**-2 *  (1+beta*aniso*(wce-wce0)/wce0)**-2 * (beta*aniso**2*((wce+wce0-2*wl)*(wce-wce0))/wce0**2 + (1+beta)*aniso*(wce0-wl)/wce0 -1)
    return gml

#checked by Tao 2017
def OMG_Tr(w,wpe,wce,v_perp,B1,delta,lorentz_f):
    #calculate k v_1 in unit of gyro frequency
    #k from linear dispersion relation, given omega in Omega_0 return k in Omega_0/c unit
    #v_1 in c
    kl=f_k(w,wpe,wce)
    return np.sqrt(kl*v_perp*B1*wce)*delta*lorentz_f**-0.5

#our way to calculate trapping frequency
#J is a constant load from file, O is approximate 0, Pcr is constant
#in unit wpe only
#wrong not from linear
def OMG_Tr2(wl,wpe,wce,Jact,a):
    kl=linear_k(wl,wpe,wce)
    pcr= (wl - wce)/kl/kl
    return np.sqrt(np.sqrt(2*wce*(Jact + pcr)) * kl**2 * a)

def fit_gm(B):
    a1=0.0371
    a2=0.76
    return (B/a1)**(1/a2)

## calculate basic data on a magnetic field line of Earth's dipole field.
'''
def zpos():

    return 0
def wce(zpos):
    return 0
def wpe(zpos):
    return 0

#def ord2(zpos,z0,dkdt,vr,vg):
#    return 0.5*(-1/vr-1/vg)*dkdt*(zpos-z0)**2
#exp2nd = ord2(zpos,zpos[512],dkdt[:,:],vr,vg)


#equilibrium function
#  subroutine set_feq(eon, chorus)
#    ! energetic electron equilibrium distribution for one macro-particle
#    implicit none
#    type(energetic_eon), intent(inout) :: eon
#    type(chorus_mode), intent(in) :: chorus
#    real(fp) :: mu, Jdist1, Jdist2
#    integer :: j
#    do j = 0, Np
#       mu = max(eon%Jpos+chorus%pcr+eon%pcor(j), small)
#       Jdist1 = exp(-mu*gyro0/(vperp*vperp))  
#       Jdist2 = exp(-mu*gyro0/(loss_cone*vperp*vperp))
#       eon%feq(j) = gyro0/((2.0_fp*pi)**1.5_fp*vperp*vperp*vll)/(1.0_fp-loss_cone)        &
#                  * exp(-0.5_fp*(chorus%kmode*(eon%pcor(j)+chorus%pcr)/vll)**2)           &
#                  * exp(-mu*(chorus%gyro-gyro0)/(vll*vll))*(Jdist1-Jdist2)
#       eon%dfeq(j) = -eon%feq(j)*((chorus%kmode**2*eon%pcor(j)+chorus%omega-gyro0)        &
#                   / (vll*vll)+gyro0/(loss_cone*vperp*vperp)                              &
#                   * (loss_cone*Jdist1-Jdist2)/(Jdist1-Jdist2)) 
#    end do
#    ! delta f 
#    eon%feon = 0.0_fp
#    return
#  end subroutine set_feq



#verify consistency
#def get_feq(vperp,vpara,gyro0,gyro,omega,kl,Jact,Omega,loss_cone=1e-16):
#    small=1e-10
#    PI=(omega-gyro)/kl/kl
#    mu = Jact+PI+Omega
#    mu[mu<0] = small
#    feq = gyro0/((2.0*np.pi)**1.5*vperp*vperp*vpara) *np.exp(-mu*gyro0/vperp**2)*np.exp(-mu*(gyro-gyro0)/vpara**2)*np.exp(-0.5*(kl*(PI+Omega)/vpara)**2)
#    dfeq =feq * (-gyro0/vperp**2 - (gyro-gyro0)/vpara**2 - kl**2*(PI+Omega)/vpara**2)
#    Jdist1 = np.exp(-mu*gyro0/(vperp*vperp)) 
#    Jdist2 = np.exp(-mu*gyro0/(loss_cone*vperp*vperp))
#    feq1= gyro0/((2.0*np.pi)**1.5*vperp*vperp*vpara)/(1.0-loss_cone) * np.exp(-0.5*(kl*(Omega+PI)/vpara)**2) * np.exp(-mu*(gyro-gyro0)/(vpara*vpara))*(Jdist1-Jdist2) 
#    dfeq1 = - feq1*((kl**2*Omega+omega-gyro0) / (vpara*vpara)+gyro0/(loss_cone*vperp*vperp)*(loss_cone*Jdist1-Jdist2)/(Jdist1-Jdist2)) 
#    Jdist2 = 0.0
#    feq2 = gyro0/((2*np.pi)**1.5*vperp*vperp*vpara)* np.exp(-0.5*(kl*(Omega+PI)/vpara)**2)*np.exp(-mu*(gyro-gyro0)/(vpara*vpara))*(Jdist1-Jdist2)
#    dfeq2 = -feq2 *((kl**2*Omega+omega-gyro0)/(vpara*vpara)+gyro0/(vperp*vperp))
#    return (feq,dfeq),(feq1,dfeq1),(feq2,dfeq2)



def get_feq(vperp,vpara,gyro0,gyro,omega,kl,Jact,Omega,loss_cone=1e-16):
    small=1e-10
    PI=(omega-gyro)/kl/kl
    mu = Jact+PI+Omega
    mu[mu<0] = small
    feq = gyro0/((2.0*np.pi)**1.5*vperp*vperp*vpara) *np.exp(-mu*gyro0/vperp**2)*np.exp(-mu*(gyro-gyro0)/vpara**2)*np.exp(-0.5*(kl*(PI+Omega)/vpara)**2)
    dfeq =feq * (-gyro0/vperp**2 - (gyro-gyro0)/vpara**2 - kl**2*(PI+Omega)/vpara**2)
    return (feq,dfeq)
'''



# ---------------------------------------------------------------------------------------- wave growth rate analysis
import sympy as sp
from scipy.optimize import fsolve
#all frequency normalized to wpe

#construct gamma and its derivative function
#in code
#a parallel velocity
#a perpendicular velocity
#c density at equator
#d loss cone
def funcG(a,b,c,d):
    wl,wpe,wce,wce0,vpar,vper,nh,beta = sp.symbols('wl,wpe,wce,wce0,vpar,vper,nh,beta')
    kl = sp.sqrt(wl**2 + wpe**2*wl/(wce - wl))
    vg = 2*kl/(2*wl + wpe**2*wce/(wce-wl)**2)
    temp =  (vper/vpar)**2
    T1 = sp.sqrt(2*np.pi)*wce*vg*nh/(4*kl**2*vpar)
    T2 = sp.exp(-0.5*(wl-wce)**2/kl**2/vpar**2)
    T3 = (1+temp*(wce-wce0)/wce0)**-2
    T4 = (1+beta*temp*(wce-wce0)/wce0)**-2
    T5 = (beta*temp**2*(wce+wce0-2*wl)*(wce-wce0)/wce0**2)
    T6 = (1+beta)*temp*(wce0-wl)/wce0 
    T7 = -1
    gamma = T1*T2*T3*T4*(T5+T6+T7)
    dgamma = sp.diff(gamma,wl)
    ddgamma = sp.diff(dgamma,wl)
    fgamma = sp.lambdify([wl,wce,wce0,wpe],gamma.subs({vpar:a,vper:b,nh:c,beta:d}))
    fdgamma = sp.lambdify([wl,wce,wce0,wpe],dgamma.subs({vpar:a,vper:b,nh:c,beta:d}))
    fddgamma = sp.lambdify([wl,wce,wce0,wpe],ddgamma.subs({vpar:a,vper:b,nh:c,beta:d}))
    return fgamma,fdgamma,fddgamma    

#find the most unstable frequency
from scipy import optimize 
def solv_w(f,jac,wce,wce0,wpe,ig):
    g = lambda x: f(x, wce, wce0, wpe)
    dg = lambda x: jac(x, wce, wce0, wpe)
    sol = optimize.root(g, ig, jac=dg, method= 'hybr')
    return sol.x

def growth(fce,fce0,fpe,vpa,vpe,nh,beta=0,ig=0.06):
    fgamma,fdgamma,fddgamma = funcG(vpa,vpe,nh,beta)
    f_m = solv_w(fdgamma,fddgamma,fce,fce0,fpe,ig)
    k_m = f_k(f_m,fce,fpe)
    gm = fgamma(f_m,fce,fce0,fpe)
    return f_m,k_m,gm

#sololy gamma
def f_gamma(wl,wce0,wce,wpe,vpar,vper,nh,beta=0):
    kl = np.sqrt(wl**2 + wpe**2*wl/(wce - wl))
    vg = 2*kl/(2*wl + wpe**2*wce/(wce-wl)**2)
    temp =  (vper/vpar)**2
    T1 = np.sqrt(2*np.pi)*wce*vg*nh/(4*kl**2*vpar)
    T2 = np.exp(-0.5*(wl-wce)**2/kl**2/vpar**2)
    T3 = (1+temp*(wce-wce0)/wce0)**-2
    T4 = (1+beta*temp*(wce-wce0)/wce0)**-2
    T5 = (beta*temp**2*(wce+wce0-2*wl)*(wce-wce0)/wce0**2)
    T6 = (1+beta)*temp*(wce0-wl)/wce0 
    T7 = -1
    gamma = T1*T2*T3*T4*(T5+T6+T7)
    return gamma


def B_evlp(compa,dz,kmode):
    paps=np.gradient(compa,axis=1)/dz
    ika=1j*kmode*compa
    B = 1j*(paps - ika)
    return B

def E_evlp(compa,dT,w):
    papt=np.gradient(compa,axis=0)/dT
    ioa=1j*w*compa
    E = -(papt + ioa)
    return E

def Jh_evlp(comps):
    return -comps/4/np.pi

def Jc_evlp(E,w,wpe,wce):
    chi = wpe**2/w/(wce-w)
    J = -1j*w*chi*E/4/np.pi
    return J
from chorus import convert


def cpst_j(j,c):
    phase_w = convert.g_phase(c)
    cpst_j = j * np.exp(-phase_w*1j)
    #real is J_B while imagnary is J_E
    return cpst_j

#da^2/dt
def power_trans(j,a,vg,kl):
    return 4 * np.pi * vg / kl * np.real(-1j*j*np.conj(a))

def delta_phi(j,a):
    phi_j = np.angle(j)
    phi_a = np.angle(a)
    return phi_j - phi_a
