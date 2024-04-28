zpos=np.load("zpos.out")
kl=np.load("kmode.out")
gyro=np.load("gyro.out")
gm=np.load("gamma.out")
wp2=np.load("wp2.out")
vr = np.load("vr.out")
vg = np.load("vg.out")
chor= np.load("chor.npy")
sour= np.load("sour.npy")
#fixed data
wl = 0.30328825168137874/5
gamma_l = 5.3357799456513026E-003
dT0=48.898864632961590
chor=np.load("chor.npy")
compchor= chor[...,0] + chor[...,1]*1j
phase_chor = np.angle(compchor)
wall,klin,knl,kall,vg1t,vg0t,vrt = cl.calcWK(wl,gyro,wp2,kl,phase_chor[:,1:],dT0,zpos)
dgyro = np.gradient(gyro)/np.gradient(zpos)
dk = np.gradient(kall,axis=1)/np.gradient(zpos)
Jact = cl.f_Jact(wl,gyro,kl,0.2,0.3,kl[int(kl.size/2)])
PI = (wl - gyro)/kl**2
dT=48.9
dwall=np.gradient(wall,axis=0)/dT0
ampchor=np.linalg.norm(chor,axis=-1)
#construct 998 1000
_a = ampchor[:,1:]
_wb2 = cl.f_wb2(gyro,Jact,_a,wall,kall)

holepos=(wl-gyro)/kall/kall - (wl - gyro)/kl/kl
ALP_w = cl.f_alpha_w(dwall,dgyro,_wb2,kall,kl,Jact,PI)

def f_alpha_w(freq_t,gyro_s,wb2,k,kl,J,PI):
    #nonlinear dp: vg = - vr
    alp = 1/wb2 *( 4 * freq_t + (J + (kl/2/k - 1)*PI)*k * gyro_s)
    return alp
def f_alpha_t(freq,freq_t,gyro,gyro_s,wb2,k,vr,vg,J):
    mu = J + (freq - gyro)/k/k
	#linear vg and vr
    alp = 1/wb2*( (1-vr/vg)**2 * freq_t + (k*mu - 3*vr/2) * gyro_s)
    return alp
def f_alpha(freq_t,gyro_s,wb2,k,k_s,vr,vg,J):
	# vg and vr imsimulation
    alp = 1/wb2 * ( (1-2*vr/vg) * freq_t - vr**2*k_s + (J)*k * gyro_s )
    return alp

ALP_t= f_alpha_t(wall,dwall,gyro,dgyro,_wb2,kall,vrt,vg0t,Jact)
ALP_d= f_alpha(dwall,dgyro,_wb2,kall,dk,vrt,-vrt,Jact)
fs = signal.savgol_filter(wall[:,:],169,1,axis=0)
dfs=np.gradient(fs,axis=0)/dT0
ALP_t1= f_alpha_t(fs,dfs,gyro,dgyro,_wb2,kall,vrt,vg0t,Jact)
ALP_d1= f_alpha(dfs,dgyro,_wb2,kall,dk,vrt,-vrt,Jact)
def c_rate_t(freq,gyro,gyro_s,k,vr,vg,J):
    mu = J + (freq - gyro)/k/k
	c_rate = (1-vr/vg)**-2 * (k*mu - 3*vr/2) * gyro_s
    return c_rate


def c_rate_full(alp,wb2,k,k_s,gyro_s,vr,vg,J):
	c_rate = (alp*wb2 + vr**2*k_s - J*k* gyro_s)/(1-2*vr/vg)
    return c_rate
c_rate_t(wall,gyro,dgyro,kall,vrt,vg0t,Jact)

c_rate_full(0.1,_wb2,kall,dk,dgyro,vrt,-vrt,Jact)
#wb2
#upshift omega
#upshift=(_f-_g)/_k - (wl - _g)/_kl
#holepos=(wl-_g)/_k - (wl - _g)/_kl
upshift=(wall-gyro)/kall - (wl - gyro)/kl
holepos=(wl-gyro)/kall - (wl - gyro)/kl

sq2 = np.sqrt(2)
_fp = 1/sq2*_f + (1-1/sq2)*_g + (1-1/sq2)*_k*_k * _PI

