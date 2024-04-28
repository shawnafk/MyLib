import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import spectrogram
waveamp=np.load("chor.npy")[:,:,0]
f,t,S=spectrogram(waveamp[:,531],nperseg=256,noverlap=255)
dT = 49.6
time,freq = np.meshgrid(t*dT,f/dT*2*np.pi)
plt.pcolormesh(time,freq,(S),shading='gouraud')
