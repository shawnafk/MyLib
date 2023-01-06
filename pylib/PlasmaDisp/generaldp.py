from sympic_io.tools import unit
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
mpl.rcParams['legend.edgecolor']= '1'
mpl.rcParams['legend.framealpha']= '0'
mpl.rcParams['xtick.direction']= 'in'
mpl.rcParams['ytick.direction']= 'in'
mpl.rcParams['xtick.labelsize']=15
mpl.rcParams['ytick.labelsize']=15
#mpl.rc('text', usetex=True)


def nkv(case,n_range,B,kp,w):
    n=[]
    kv=[]
    for __n in n_range:
        for __kv in case.omegakv(__n,B,kp,w):
            n.append(__n)
            kv.append(__kv)
    DP=zip(kv,n)
    res=sorted(DP,key= lambda v:v[1])
    kv,n=zip(*res)
    return (kv,n)
def nkv_hot(case,n_range,B,kp,w,vth):
    n=[]
    kv=[]
    for __n in n_range:
        for __kv in case.omegakv_hot(__n,B,kp,w,vth):
            n.append(__n)
            kv.append(__kv)
    DP=zip(kv,n)
    res=sorted(DP,key= lambda v:v[1])
    kv,n=zip(*res)
    return (kv,n)

def bkv(case,n,B_range,kp,w):
    b=[]
    kv=[]
    for __b in B_range:
        for __kv in case.omegakv(n,__b,kp,w):
            b.append(__b)
            kv.append(__kv)
    DP=zip(kv,b)
    res=sorted(DP,key= lambda v:v[1])
    kv,n=zip(*res)
    return (kv,b)

def bkv_hot(case,n,B_range,kp,w,vt):
    b=[]
    kv=[]
    for __b in B_range:
        for __kv in case.omegakv_hot(n,__b,kp,w,vt):
            b.append(__b)
            kv.append(__kv)
    DP=zip(kv,b)
    res=sorted(DP,key= lambda v:v[1])
    kv,n=zip(*res)
    return (kv,b)


def nbkv_slow(case,n_range,B_range,kp,w):
    ln=len(n_range)
    lb=len(B_range)
    sk=np.zeros((ln,lb))
    for idn in range(ln):
        for idb in range(lb):
            res=case.omegakv(n_range[idn],B_range[idb],kp,w)
            if(res):
                sk[idn,idb]=max(res)
    return sk

def wkv_hot(case,n,B,kp,w_range,vt):
    w=[]
    kv=[]
    for __w in w_range:
        for __kv in case.omegakv_hot(n,B,kp,__w,vt):
            w.append(__w)
            kv.append(__kv)
    DP=zip(kv,w)
    res=sorted(DP,key= lambda v:v[0])
    kv,w=zip(*res)
    return (kv,w)

L=1e-3
U=unit.gen_unit(L)
case=lht.LH_Dp(1,1,U)
'''
kp=2*np.pi/384
#kv1,w1=wkv(case,1e17,0.4,kp,np.linspace(0.00001,0.012,10000))
kv1,w1=wkv(case,1e18,0.2,kp,np.linspace(0.00001,0.012,10000))
kv2,w2=wkv(case,1e18,0.4,kp,np.linspace(0.00001,0.012,10000))
kv3,w3=wkv(case,1e18,3,kp,np.linspace(0.00001,0.02,10000))
f1,ax1=plt.subplots()
ax1.semilogx(kv1,w1,label='0.2,1e18')
ax1.semilogx(kv2,w2,label='0.4,1e18')
ax1.semilogx(kv3,w3,label='3,1e18')
ax1.legend()

kvn1,n1=nkv(case,np.linspace(0.001,1.5,1000)*1e18,0.3,kp,0.0037)
kvn2,n2=nkv(case,np.linspace(0.001,1.5,1000)*1e18,0.4,kp,0.0037)
kvn3,n3=nkv(case,np.linspace(0.001,1.5,1000)*1e18,1.6,kp,0.0037)
#kvb,b=bkv(case,1e18,np.linspace(0.2,2,100),kp,0.0037)

f2,ax2=plt.subplots()
ax2.scatter(n1,np.log(kvn1),1.5,label='0.3')
ax2.scatter(n2,np.log(kvn2),1.5,label='0.4')
ax2.scatter(n3,np.log(kvn3),1.5,label='1.6')

ax2.legend()

plt.show()
'''

'''
kw=np.load('kv_w_pulse.npy')
x,y=kw.shape
x=int(x/2)
y=int(y/2)
X=np.linspace(0,np.pi,x)
Y=np.linspace(0,np.pi/(0.5*200),y)

K,W=np.meshgrid(X,Y)
plt.contour(K,W,abs(np.fft.fftshift(kw)[x:2*x,y:2*y]).transpose(),levels=3)

w=0.0043
dx=1e-3;U=gen_unit(dx)
mu=1;Z=1
#ne_pork=1e18*np.linspace(1e-2,1,1000)
B0=0.23;
B=np.linspace(B0,2.5,100)
kp=2*pi/(384)
wl=wz(w,kp,1e18,B,mu,Z,U)
#wl2=wz(w,kp,2*ne_pork,B,mu,Z,U)

plot(B,wl)
show()
'''

def k2w(n,B,k,ct,m,q,U):
    wce=unit.Oce(B,U)
    wci=unit.Oci(m,q,B,U)
    wpe=unit.Ope(n,U)
    wpi=unit.Opi(m,q,n,U)
    wp=(wpe**2 +wpi**2)**0.5
    c=1
    c10 = 1
    c8  = 2*c**2 * k**2 + wce**2 + wci**2 + 3*wp**2 
    
    c6  =   c**4 * k**4 + (2* c**2 * k**2 + wp**2 )*(wce**2 + wci**2 + 2*wp**2) + (wp**2 + wce*wci)**2
    
    c4  =   c**4 * k**4 * (wce**2 + wci**2 + wp**2) + 2*c**2 * k**2 * (wp**2 + wci*wce)**2 + \
            c**2 * k**2 * wp**2 * (wce**2 + wci**2 -wci*wce)*(1+ct**2) +  wp**2 * (wp**2+wci*wce)**2
    
    c2  =   c**4 * k**4 * (wp**2* (wce**2 + wci**2 -wci*wce)*ct**2 + wce*wci*(wp**2 + wci*wce)) +\
            c**2 * k**2 * wp**2*wci*wce*(wp**2 + wci*wce)*(1+ct**2)
    
    c0  =   c**4 * k**4 * wce**2 *wci**2 * wp**2 *ct**2
    w= np.sqrt(np.roots([c10,-c8,c6,-c4,c2,-c0]))
    return w
#fix k
#find theta by match omega
#get k//
kx=8*np.pi/160
kz=5/4*np.pi/512
k=(kx**2+kz**2)**0.5
N=1000
ret_w=np.zeros((5,N))
i=0


theta= np.linspace(0,np.pi/2,N) 
for __theta in theta: 
    ret=k2w(1e18,1,k,np.cos(__theta),1,1,U)
    ret_w[:,i]=np.sort(ret)
    i+=1

#fix k// find kvert and verify
Nv=1000
kp=2*np.pi/1024
kv=np.linspace(0.001,np.pi/2,Nv)
ret_w=np.zeros((5,Nv))
i=0
for __kv in kv: 
    k=(kp**2+__kv**2)**0.5
    ret=k2w(1e18,1.5,k,kp/k,1,1,U)
    ret_w[:,i]=np.sort(ret)
    i=i+1


'''
Np=100
Nv=100
kp=np.linspace(0.001,np.pi/2,Np)
kv=np.linspace(0.001,np.pi/2,Nv)
ret_w=np.zeros((5,Np,Nv))
i=0
for __kv in kv: 
    j=0
    for __kp in kp: 
        k=(__kp**2+__kv**2)**0.5
        ret=k2w(1e18,1.5,k,__kp/k,1,1,U)
        ret_w[:,i,j]=np.sort(ret)
        j+=1
    i=i+1
'''
#kv1,w1=wkv(case,1e18,1.6,kp,np.linspace(0.00001,0.05,1000))
