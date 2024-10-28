import numpy as np
#physical constant in SI unit
mu0=4*np.pi*1e-7
epsi0=8.8542e-12
c=299792458.
me=9.10000000000000006e-31
#Elementary charge in C
qe = 1.66022e-19

#n in cm**-3
#omega in Rad/s
def n_wpe(n):
    return 5.64e4*np.sqrt(n)

def wpe_n(wpe):
    return (wpe/5.64e4)**2

#B in gauss
def B_wce(b):
    return 1.76e7*b

def wce_B(wce):
    return wce/1.76e7

#m = mi/mp
#q charge state
def n_wpi(n,q=1,m=1):
    return 1.32e3*q*(m**-0.5)*np.sqrt(n)

def wpi_n(wpi,q=1,m=1):
    return (wpi/(1.32e3*q*(m**-0.5)))**2

#B in gauss
def B_wci(b,q=1,m=1):
    return 9.58e3*q*(m**-1)*b

def wci_B(wce,q=1,m=1):
    return wce/(9.58e3*q*(m**-1))



#n in cm**-3, T in eV, B in gauss
def nTB_beta(n,T,B):
    #T * qe in Joules
    #pressure in Pa
    pP = (n * 1e6) * T * qe 
    pB = (B/1e4)**2/2/mu0
    return pP/pB

def nBb_T(n,B,beta):
    #T * qe in Joules
    #pressure in Pa
    pB = (B/1e4)**2/2/mu0
    pP = beta * pB
    T = pP / (n*1e6) /qe
    return T
