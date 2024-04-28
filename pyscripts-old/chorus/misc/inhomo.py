import numpy as np
def earth(L):
    #gauss
    B_g = 3.12e-1
    B_0 = B_g/L**3
    wce = (1.76e7*B_0)
    #m/s
    c = 299792458
    R_e = 6.375e6
    #unit_L
    u_L = c/wce
    R_e_dl = R_e/u_L
    return 4.5/(L*R_e_dl)
    
    
