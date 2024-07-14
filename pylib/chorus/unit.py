import numpy as np
#unit used in Vlasov code
from plasma import nrl
#cgs constant, contained in nrl later
c_cgs = 2.9979e10
qe_cgs = 4.8032e-10
me_cgs = 9.1084e-28
#r ratio wpe/wce
#given B0 determine wce and wpe by given ratio, then determine ne to calculatet nqv
def gen_cgs(wce0,wpe0):
    units={}
    #B00 = 0.312
    #B0 = B00/L**3
    #wce0 = 1.76e7*B0
    
    #wpe = (4*np.pi * ne * qe**2 / me) **0.5 
    #ne = (wpe0/5.64e4)**2
    ne = wpe0**2*me_cgs / qe_cgs**2 /4/np.pi
    #wpe0=r*wce0
    units['w'] = wpe0
    units['t'] = 1/units['w']
    units['k'] = wpe0/c_cgs
    units['l'] = 1/units['k']
    units['ne'] = ne
    units['b'] = wpe0/qe_cgs*me_cgs*c_cgs
    #e and b are same unit
    units['e'] = units['b']
    units['a'] = me_cgs*c_cgs**2/qe_cgs
    units['j'] = wpe0 * units['e']
    consts={}
    #j/j0 = 4 pi
    consts['j0'] = ne * qe_cgs * c_cgs
    consts['b0'] = wce0/qe_cgs*me_cgs*c_cgs
    return units,consts

def gen_si(wce0,wpe0):
    units={}
    #B00 = 0.312
    #B0 = B00/L**3
    #wce0 = 1.76e7*B0
    #wpe = (4*np.pi * ne * qe**2 / me) **0.5 
    #ne = (wpe0/5.64e4)**2
    ne = wpe0**2*me_cgs / qe_cgs**2 /4/np.pi
    #wpe0=r*wce0
    units['w'] = wpe0
    units['t'] = 1/units['w']
    units['k'] = wpe0/c_cgs*(1e2)
    units['l'] = 1/units['k']
    units['ne'] = ne*1e6
    units['b'] = wpe0/qe_cgs*me_cgs*c_cgs/1e4
    #e and b are same unit
    units['e'] = units['b']/(1/3*1e-4)
    units['a'] = me_cgs*c_cgs**2/qe_cgs/1e6
    units['j'] = wpe0 * units['e']/(3e5)
    consts={}
    #j/j0 = 4 pi
    consts['j0'] = ne * qe_cgs * c_cgs/(3e5)
    consts['b0'] = wce0/qe_cgs*me_cgs*c_cgs/1e4
    return units,consts


Uchar={
't':"$~/~\omega_{\mathrm{pe0}}^{-1}$",
'f':"$~/~\omega_{\mathrm{pe0}}$",
's':"$~\cdot~\omega_{\mathrm{pe0}}/\mathrm{c}$",
'k':"$~\cdot~\mathrm{c}/\omega_{\mathrm{pe0}}$",
'v':"$~/~\mathrm{c}$"
}
