from matplotlib import pyplot as plt
import numpy as np
def wk(a,dt,ds):
    f,ax=plt.subplots()
    n,s = a.shape
    f_u = np.pi/dt
    f_d = -f_u
    k_u = np.pi/ds
    k_d = -k_u
    f_l = int(n/4)
    f_r = int(3*n/4)
    k_l = int(s/4)
    k_r = int(3*s/4)
    im = ax.imshow(np.log(np.abs(np.fft.fftshift(np.fft.fftn(a[f_l:f_r,k_l:k_r])))),extent= [k_d,k_u,f_d,f_u],aspect='auto')
    plt.colorbar(im)

from scipy.signal import spectrogram
def draw_spectrum(seq,fs,ws,novlp,low=-20,high=-10,wl=0.0606):
    f,ax=plt.subplots()
    f,t,S=spectrogram(seq,fs,nperseg=ws,noverlap=novlp,scaling='density')
    T,F=np.meshgrid(t,f)
    im = ax.contour(T,wl+F*2*np.pi,np.log(np.abs(S[:,:])),np.linspace(low,high,200),cmap='jet')
    plt.colorbar(im)
    return im
