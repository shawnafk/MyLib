from chorus import wave
from chorus import ptc
import numpy as np

from scipy.signal import spectrogram
def draw_spectrum(ax,seq,fs,ws,novlp):
    f,t,S=spectrogram(seq,fs,nperseg=ws,noverlap=novlp,scaling='spectrum')
    T,F=np.meshgrid(t,f)
    im = ax.contourf(T,wl+F*2*np.pi,np.sqrt(S[:,:]),levels=150,cmap='jet')
    return im


def draw_boundary(ax,S,wtr,kapa,shift,offset=0):
    zeta_c,zeta_x,_ = ptc.find_zeta(S) 
    zeta=np.sort([zeta_c,zeta_x])
    zetab = np.linspace(zeta[0],zeta[1],1000)
    const_E = np.cos(zeta_x) + S * zeta_x
    Omega_b1 = wtr*np.sqrt(2 * (const_E - (np.cos(zetab) + S*zetab)) )/kapa**2
    Omega_b2 = - Omega_b1
    ax.plot(zetab+offset,Omega_b1+shift,'k--',lw=1.5)
    ax.plot(zetab+offset,Omega_b2+shift,'k--',lw=1.5)
    zeta_x1 = zeta_x + offset
    #zeta_o1 = zeta_o + offset
    zeta_c1 = zeta_c + offset
    ax.scatter(zeta_x1,shift,s=150,marker='x',color='k',lw=1.5)
    ax.scatter(zeta_c1,shift,s=150,marker='^',facecolors='none', edgecolors='r',color='k',lw=1.5)
    return 0

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



def source_term(gm,compa,vg,kl):
    return -1j * gm * compa / 2/np.pi/vg*k
