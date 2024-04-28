import numpy as np
from scipy.optimize import root

#given frequency f at one preiod, calculate k from linear whistler dispersion relation
#then regard that k as the same k for the nonlinear wave, since the source is the same
#calculated nonlinear frequency with the most unstablity omega_l which is the start frequency and the corresponding k_l
#and get the frequency differency

# w: wave frequency at the center of that subpacket, or using averaged frequency by counting the wave periods
# wm: the start frequency, which is the most unstable frequency
# wce and wpe: using local one

def delta_f(wl,wm,wce,wpe):
    kl = np.sqrt(wl**2 + wpe**2*wl/(wce-wl))
    knl = kl
    km= np.sqrt(wm**2 + wpe**2*wm/(wce-wm))
    wnl = wce - (wce-wm)*km/knl
    df = wnl-wl
    num  = wnl/df
    return df, num

def D(w,k,wce,wpe):
    return k**2 - w**2 - wpe**2 * w/(wce -w )
def delta_f2(wnl,wm,wce,wpe):
    km= np.sqrt(wm**2 + wpe**2*wm/(wce-wm))
    knl = km * (wce - wm)/(wce - wnl)
    kl=knl
    root_w = root(D,wm,(kl,wce,wpe)).x
    wl = root_w[0]
    print(wl)
    #if(len(root_w)==1):
    #    wl = root_w[0]
    #    print(wl)
    #else:
    #    print(root_w)
    #    exit(1)
    df = wnl-wl
    num  = wnl/df
    return df, num





#for 22.48.20

wce0=2900
#wce = fce/fce0
wce = 3310/wce0
#wpe = ratio * wce
wpe = 9.51 * wce
wm = .246

def check_dw(duration,num,wce0=2900):
    #center frequency
    w_ave = 1/(duration/num)/wce0
    #delta freq from numbers in modulation
    dw = w_ave/num
    #df_c,nc = delta_f2 (w_ave,wm,wce,wpe)
    dw_t = delta_f(w_ave-2*dw,wm,wce,wpe)
    return dw/dw_t
#C
check_dw(73e-3,55)
#D
check_dw(21e-3,19)
#w_c = (0.246 + 0.277)/2
#w_c = 1/(73e-3/55)/2900
#df_c,nc = delta_f(w_c,wm,wce,wpe)
##df2_c,nc2 = delta_f2(w_c,wm,wce,wpe)
#
#w_d = (0.258 + 0.352)/2
#w_d = 1/(21e-3/19)/2900
#df_d,nd = delta_f(w_d,wm,wce,wpe)
##df2_d,nd2 = delta_f2(w_d,wm,wce,wpe)
#
#
#w_g = 1/(35*0.02/115/8)/2900
#df_g,ng = delta_f(w_g,wm,wce,wpe)
#
#
#
#w_e = 1/(70*0.02/93/15)/2900
#df_e,ne = delta_f(w_e,wm,wce,wpe)
#

#assume chorus and whistler have same k
#given whistler frequency w0 and a frequency dw
#we can elimnate k in the dispersion relations get beat frequency as an implcit function of frequency w0 
def codp(dw,w0,wm,wce,wpe):
    km = np.sqrt(wm**2 + wpe**2*wm/(wce-wm))
    return ((wce-wm)/(wce-(w0+dw))*km)**2 - w0**2 - wpe**2*w0/(wce-w0)

def beat_w_c(freq,wm,wce,wpe):
    root_dw = root(codp,0,(freq,wm,wce,wpe)).x[0]
    return freq/root_dw
   
beats = []
freqs= np.linspace(0.258,0.5,1000)
for freq in freqs:
    beats.append(beat_w_c(freq,wm+0.01,wce,wpe))


w_c = 1/(73e-3/55)/2900
w_d = 1/(21e-3/19)/2900
w_e = 1/(70*0.02/93/15)/2900
w_g = 1/(35*0.02/115/8)/2900

plot(freqs,beats)
scatter(w_c,55,s=30,marker='+')
scatter(w_d,19,s=30,marker='+')
scatter(w_e,15,s=30,marker='+')
scatter(w_g,8,s=30,marker='+')

