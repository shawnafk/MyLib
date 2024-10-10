import numpy as np
import matplotlib.pyplot  as plt
def plot_earth(ax):
    t = np.linspace(0,2*np.pi,128)
    x = np.cos(t)
    y = np.sin(t)
    ax.plot(x,y,'k')
    return 0

#need bow shock, magnetospause boundary
def BS_args(B_z,beta,M_ms,D_p):
    a_1=11.1266
    a_2=0.0010
    a_3=-0.0005
    a_4=2.5966
    a_5=0.8182
    a_6=-0.0170
    a_7=-0.0122
    a_8=1.3007
    a_9=-0.0049
    a_10=-0.0328
    a_11=6.047
    a_12=1.029
    a_13=0.0231
    a_14=-0.002
    if B_z >0:
        r_0=a_1 * (1+a_2*B_z)*(1+a_9 * beta ) * (1+a_4 * ((a_8-1)*M_ms**2+2)/ (a_8+1)/M_ms**2) * D_p**(-1 / a_11)
        alpha = a_5 * (1+a_13 * B_z)*(1+ a_7 * D_p)*(1+a_10*np.log(1+beta)) * (1+a_14 * M_ms)
    else:
        r_0=a_1 * (1+a_3*B_z)*(1+a_9 * beta ) * (1+a_4 * ((a_8-1)*M_ms**2+2)/ (a_8+1)/M_ms**2) * D_p**(-1/a_11)
        alpha=a_5*(1+a_6 * B_z)*(1+ a_7 * D_p)*(1+a_10*np.log(1+beta)) * (1+a_14 * M_ms)
    epsi = a_12
    return r_0, alpha,epsi

def BS(B_z,beta,M_ms,D_p):
    r0,alp,epsi =  BS_args(B_z,beta,M_ms,D_p)
    theta = np.linspace(0,2*np.pi,128)
    r = r0 * ( (1+epsi)/(1+epsi * np.cos(theta) ) ) ** alp
    x = r * np.cos(theta)
    R = r * np.sin(theta)
    return x,R

def MP_args(B_z,D_p):
    a_1 = 11.646
    a_2 = 0.216
    a_3 = 0.122
    a_4 = 6.215
    a_5 = 0.578
    a_6 = -0.009
    a_7 = 0.012
    if B_z > 0:
        r0  = a_1 * D_p**(-1/a_4)
    elif B_z < 0 and B_z > -8:
        r0 = (a_1 + a_2* B_z) * D_p**(-1/a_4)
    else:
        r0 =( a_1 + 8*a_3 - 8 * a_2 + a_3 * B_z ) * D_p**(-1/a_4)
    alpha = (a_5 + a_6 * B_z) * (1+ a_7 * D_p)
    return r0, alpha

def MP(B_z,D_p):
    r0,alp = MP_args(B_z,D_p)
    theta = np.linspace(0,2*np.pi,128)
    r = r0 * ( 2/(1+np.cos(theta) ) ) ** alp
    x = r * np.cos(theta)
    R = r * np.sin(theta)
    return x,R
