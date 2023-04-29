import numpy as np
from scipy import optimize
from scipy import integrate
from chebpy import chebfun

# define wave potential
def potential(wtr, S, zeta):
    return 2 * wtr**2 * (np.cos(zeta) - S * zeta)

# take wtr S and return (zeta0,zeta1,zeta2) which are the lowest potential point, X point, and O point
def find_zeta(S):
    # wtr can be any value
    wtr = 1
    zeta_x = np.mod(-np.arcsin(S), 2*np.pi)
    const_E = potential(wtr, S, zeta_x)
    # use chebfun function to find the other intersection (zeta_o) of "E=const_E" line and "E = potential" curve
    sol = chebfun(lambda zeta: 2*wtr**2 *
                  (np.cos(zeta) - S*zeta)-(const_E), [0, 2*np.pi])
    allzeta = sol.roots()
    # allzeta result may contain both zeta_o and zeta_x, choose the one with largest distance to our known zeta_x as zeta_o
    if (len(allzeta) != 1):
        d2zeta_x = np.abs(allzeta - zeta_x)
        zeta_o = allzeta[np.where(d2zeta_x == np.max(d2zeta_x))][-1]
    else:
        zeta_o = allzeta[0]
        zeta_o = np.array(zeta_o).flatten()[-1]
        zeta_x = np.array(zeta_x).flatten()[-1]
    return zeta_o, zeta_x

# integrate between O and X points
# mspx = int d\xi cos \xi \sqrt(E_spx - cos \xi - \alpha \xi)
# and
# jspx  = int d\xi \sqrt(E_spx - cos \xi - \alpha \xi)
# E_spx = cos(zeta_1) \pm S*zeta1
# zeta_1 is the X point, pm depends on the sign of S
def MJ_spx(S):
    zeta_o, zeta_x = find_zeta(S)

    def Mspx(zeta, zeta1, S):
        return (np.cos(zeta1) - S*zeta1 - (np.cos(zeta) - S*zeta))**0.5 * np.cos(zeta)

    def Jspx(zeta, zeta1, S):
        return (np.cos(zeta1) - S*zeta1 - (np.cos(zeta) - S*zeta))**0.5
    mspx = integrate.quad(Mspx, zeta_x, zeta_o, args=(zeta_x, S))[0]
    jspx = integrate.quad(Jspx, zeta_x, zeta_o, args=(zeta_x, S))[0]
    return jspx, mspx

# The current ratio JB/JE, JB is imag(Source) and JE is real(Source), the phase of the source has been substracted by phase of the wave.
def ratioJ_w(S):
    jspx, mspx = MJ_spx(S)
    return mspx/jspx/S

#solve alpha directly
def calcR_w(ratio):
    sol = []
    for r in ratio:
        if(abs(r) > 0.5):
            ig = -0.2
        else:
            ig = -0.99
        print(r,ig)
        _sol = optimize.root(lambda x:ratioJ_w(x) - r, ig)
        sol.append(_sol.x[0])
    return np.array(sol)

def calcR(J):
    sol = []
    def equ(x):
        j,m=MJ_spx(x)
        return x*j*np.imag(J) - m*np.real(J)
    _sol = optimize.root(lambda x:equ(x), 0.2)
    sol.append(_sol.x[0])
    return np.array(sol)

#or alternatively we can define an interpolation function between alpha and current ratio
def init_R_ratio():
    # alpha (ANS_R) is chosen to be from 1e-4 to 9.9999...
    ANS_R = np.logspace(-4, -0.0000001, 1000)
    # calculate the current ratio (ANS_RATIO_J) using ANS_R
    ANS_RATIO_J = np.array([ratioJ_w(r) for r in ANS_R])
    from scipy import interpolate
    #construct a interpolation function based on ANS_R = f(ANS_RATIO_J)
    calcR_interp = interpolate.interp1d(abs(ANS_RATIO_J), ANS_R)
    return calcR_interp

#initialize the interpolation function
clcR = init_R_ratio()

chor=np.load("chor.npy")
sour=np.load("sour.npy")
compchor= chor[...,0] + chor[...,1]*1j
compsour= sour[...,0] + sour[...,1]*1j
phase_chor = np.angle(compchor)
phase_sour = np.angle(compsour)
ampchor= np.linalg.norm(chor,axis=-1)
ampsour= np.linalg.norm(sour,axis=-1)

compensated_sour = compsour*np.exp(-1j*phase_chor)
phase_cpst = -np.angle(compensated_sour)
#JB/JE
#real is J_B while imagnary is J_E
#ratio=np.imag(compensated_sour)/np.real(compensated_sour)
ratio=np.real(compensated_sour)/np.imag(compensated_sour)
#eleminate too large/small values
ratio[abs(ratio)>1000]=1000
ratio[abs(ratio)<0.01]=0.01
interp_R = clcR(abs(ratio))


#calculate alpha by definition
def getJact(omega,gyro,kmode,gyro0,vperp,kmode0):
    pcr0 = (wl - gyro0)/kmode0**2
    Jact0 = vperp*vperp/gyro0 - pcr0
    const=(omega-gyro0)*Jact0+0.5*(omega-gyro0)**2/kmode0**2
    return (const-0.5*(omega-gyro)**2/kmode**2)/(omega-gyro)

def alpha(gyro,gyro_s,kl,J,PI,a,k,freq,freq_t):
    #mu=J+Omega+PI = J+ p_para/k
    mu = J + (freq - gyro)/k/k
    vperp = np.sqrt(2*gyro*mu)
    wb2= vperp * k**2 * a
    alp = 4/wb2 * freq_t + (J + (kl/2/k - 1)*PI)*k/wb2 * gyro_s
    return alp

wl = 0.06065765033627572
kl=np.load('kmode.npy')
zpos=np.load('zpos.npy')
gyro=np.load('gyro.npy')
dgyro = (gyro[1:]-gyro[:-1])/(zpos[1:]-zpos[:-1])

Jact = getJact(wl,gyro,kl,0.2,0.3,kl[int(kl.size/2)])
PI = (wl - gyro)/kl**2

t=800
wall=np.load("wall.npy")
kall=np.load('kall.npy')
ALP=[]
dT=48.9
dwall=(wall[1:,:]-wall[:-1,:])/dT
for z in np.arange(5500,9000,100):
    ALP.append(alpha(gyro[z],dgyro[z],kl[z],Jact[z],PI[z],ampchor[t,z],kall[t,z],wall[t,z],dwall[t,z]))