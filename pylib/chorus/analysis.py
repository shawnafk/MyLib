from chorus import wave
from chorus import io
from chorus import ptc
import numpy as np
#wfile chor..
#ifile zpos, kl, gyro, wp2
def _mid(arr):
    return (arr[1:] + arr[:-1])/2
def mid(s,ax):
    return np.apply_along_axis(_mid,ax,s)
def comp(s):
    return s[...,0]+s[...,1]*1j

def g_phase(s):
    return np.angle(s)

#return omega,wavenumber,phasew,phases,vg,wb2,wb2_cplx,alpha,drg
def g_all(res):
    #complex
    compchor= mid(comp(res['c']),1)
    compsour= mid(comp(res['s']),1)
    #phase
    phasew = g_phase(compchor)
    phases = g_phase(compsour)
    #w and k
    dz = np.gradient(res['z'])
    dw,dk = wave.calcWK(phasew,res['dt'],dz)
    omega=res['w']+dw
    wavenumber=res['k']+dk
    # vg
    vg = wave.calcVg(omega,wavenumber,res['dt'],dz)
    #particle
    dgyro = np.gradient(res['g'])/np.gradient(res['z'])
    dkmode = np.gradient(res['k'])/np.gradient(res['z'])
    dwdt =np.gradient(omega,axis=0)/res['dt']
    ampchor=np.abs(compchor)
    PI = (res['w'] - res['g'])/res['k']**2
    #wb2 amp and vect
    wb2 = ptc.f_wb2(res['g'],res['j'],ampchor,omega,wavenumber)
    wb2_cplx = ptc.f_wb2(res['g'],res['j'],ampchor,omega,wavenumber)
    #alpha and drg
    alpha = ptc.f_alpha_w(dwdt,dgyro,wb2,wavenumber,res['k'],res['j'],PI)
    #alpha1 = ptc.f_alpha_w(0,dgyro,wb2,wavenumber,res['j'],res['j'],PI)
    drg=res['j']/res['k'] * dgyro - (res['w']-res['g'])**2/res['k']**4 * dkmode
    #return omega,wavenumber,phasew,phases,vg,wb2,wb2_cplx,alpha,alpha1,drg
    return omega,wavenumber,phasew,phases,vg,wb2,wb2_cplx,alpha,drg

#interp_R=ptc.init_R_ratio()
def project(chor,sour):
    phase_chor = np.angle(compchor)
    #phase_sour = np.angle(compsour)
    #ampchor= np.abs(compchor)
    #ampsour= np.abs(compsour)
    compensated_sour = compsour*np.exp(-1j*phase_chor)
    #phase_cpst = -np.angle(compensated_sour)
    #JB/JE
    #real is J_B while imagnary is J_E
    #ratio=np.imag(compensated_sour)/np.real(compensated_sour)
    return compensated_sour
def adiabatic_R(compensated_sour,interp_R):
    ratio=np.real(compensated_sour)/np.imag(compensated_sour)
    ratio[abs(ratio)>1000]=1000
    ratio[abs(ratio)<0.01]=0.01
    R=interp_R(abs(ratio))
    return R

#frequency in unit of B
def calc_wb(w,deltaB,vp=0.6,wce=1,wpe=5):
    k = wave.f_k(w,wce,wpe)
    c=1
    mu = c*k/w
    return np.sqrt(mu* w/wce * vp * deltaB)

def cold_current(compa,dT,w,wpe,wce):
    E=-np.gradient(compa,axis=0)/dT
    chi = wpe**2/w/(wce-w)
    J = -1j*w*chi*E/4/np.pi
    return J

def B_evlp(compa,zpos,kmode):
    paps=np.gradient(compa,axis=1)/np.gradient(zpos)
    ika=1j*kmode*compa
    B = 1j*(paps - ika)
    return B

def E_evlp(compa,dT,w):
    papt=np.gradient(compa,axis=0)/dT
    ioa=1j*w*compa
    E = -(papt + ioa)
    return E

def hot_current_evlp(comps):
    return -comps/4/np.pi

def cold_current_evlp(E,w,wpe,wce):
    chi = wpe**2/w/(wce-w)
    J = -1j*w*chi*E/4/np.pi
    return J

def source_term(gm,compa,vg,kl):
    return -1j * gamma * compa / 2/pi/vg*kl
