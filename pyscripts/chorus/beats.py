import numpy as np
from scipy.optimize import root 
from matplotlib import pyplot as plt
#assume chorus and whistler have same k
#given whistler frequency w0 and a frequency dw
#we can elimnate k in the dispersion relations get beat frequency as an implcit function of frequency w0 
def codp(dw,w0,wm,wce,wpe):
    km = np.sqrt(wm**2 + wpe**2*wm/(wce-wm))
    return ((wce-wm)/(wce-(w0+dw))*km)**2 - w0**2 - wpe**2*w0/(wce-w0)

#assume chorus frequency is w1, whistler is w1  - dw
def codp2(dw,w1,wm,wce,wpe):
    km = np.sqrt(wm**2 + wpe**2*wm/(wce-wm))
    return ((wce-wm)/(wce-w1)*km)**2 - (w1-dw)**2 - wpe**2*(w1-dw)/(wce-(w1-dw))


def codp3(dw,w0,wm,wce,wpe):
    km = np.sqrt(wm**2 + wpe**2*wm/(wce-wm))
    kc = km*(wce-wm)/(wce-w0)
    vr = (w0-wce)/kc
    kw = np.sqrt((w0-dw)**2+wpe**2*(w0-dw)/(wce-w0+dw))
    return vr-(w0-dw-wce)/kw

def beat_w_c(freq,wm,wce,wpe):
    root_dw = root(codp3,0,(freq,wm,wce,wpe)).x[0]
    return freq/root_dw
 
wce0=2900
#wce = fce/fce0
wce = 3310/wce0
#wpe = ratio * wce
wpe = 100 * wce
wm = .246

beats = []
freqs= np.linspace(0.247,0.5,1000)
for freq in freqs:
    beats.append(beat_w_c(freq,wm,wce,wpe))


w_c = 1/(73e-3/55)/2900
w_d = 1/(21e-3/19)/2900
w_e = 1/(70*0.02/93/15)/2900
w_g = 1/(35*0.02/115/8)/2900


#w_c = 0.25
#w_d = 0.25
#w_e = 0.36
#w_g = 0.425
plt.semilogy(freqs,beats)
plt.scatter(w_c,55,s=50,marker='+')
plt.scatter(w_d,19,s=50,marker='+')
plt.scatter(w_e,15,s=50,marker='+')
plt.scatter(w_g,8,s=50,marker='+')
plt.ylabel('BEATS')
plt.xlabel('Frequency')
plt.show()
