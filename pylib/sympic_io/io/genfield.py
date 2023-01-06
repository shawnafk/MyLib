import scipy.io as sio
from sympic_io.tools import unit 
import scipy.special as sp
from numpy import *
import numpy as np
def is_inside(z,r):
    return (z>r[0] and z<r[1]) or z==r[0] or z==r[1] 
def preload(fname,info,DimArray):
    fid = open(fname,'wb+');
    info.tofile(fid)
    DimArray.tofile(fid)
    return fid
def siglecoil_field1(I0,N_coil,r0,x,y,z):
    c0=4e-7; # mu0/pi=4e-7
    I0=I0*N_coil; 
    # current, unit: A
    C=c0*I0;
    # Set field for single coil
    rho=(x**2+y**2)**0.5
    r=(x**2+y**2+z**2)**0.5
    alpha=(r0**2 + r**2 - 2*r0*rho)**0.5
    beta=(r0**2 + r**2 + 2*r0*rho)**0.5
    k2 = 1 - alpha**2/beta**2
    (K,E)=(sp.ellipk(k2),sp.ellipe(k2))
    Bx=C*x*z/(2*alpha**2*beta*rho**2) * ( (r0**2+r**2)*E - alpha**2 * K)
    By=C*y*z/(2*alpha**2*beta*rho**2) * ( (r0**2+r**2)*E - alpha**2 * K)
    Bz=C/(2*alpha**2*beta) * ((r0**2 - r**2)*E + alpha**2 * K)
    return (Bx,By,Bz)
def siglecoil_fB(I0,N_coil,r0,x,y,z):
    c0=4e-7; # mu0/pi=4e-7
    I0=I0*N_coil; 
    # current, unit: A
    C=c0*I0
    # Set field for single coil
    rho=(x**2+y**2)**0.5
    r=(x**2+y**2+z**2)**0.5
    alpha=(r0**2 + r**2 - 2*r0*rho)**0.5
    beta=(r0**2 + r**2 + 2*r0*rho)**0.5
    gamma=x**2-y**2
    k2 = 1 - alpha**2/beta**2
    (K,E)=(sp.ellipk(k2),sp.ellipe(k2))
    t1=r0**4*(-gamma*(3*z**2+r0**2) + rho**2*(8*x**2-y**2))
    t2=r0**2*(rho**4*(5*x**2+y**2)  - 2*rho**2*z**2*(2*x**2+y**2)
            +3*z**4*gamma) 
    t3=r**4*(2*x**4+gamma*(y**2+z**2))
    t4=r0**2*( gamma*(r0**2+2*z**2) - rho**2*(3*x**2-2*y**2) )
    t5=r**2*( 2*x**4 + gamma*(y**2 + z**2) )
    Bxx=C*z/(2*alpha**4*beta**3*rho**4)* ( (t1-t2-t3)* E + (t4+t5)*alpha**2*K)
    t1=r0**4*(gamma*(3*z**2+r0**2) + rho**2*(8*y**2-x**2))
    t2=r0**2*(rho**4*(5*y**2+x**2)  - 2*rho**2*z**2*(2*y**2+x**2)
            -3*z**4*gamma) 
    t3=r**4*(2*y**4-gamma*(x**2+z**2))
    t4=r0**2*( -gamma*(r0**2+2*z**2) - rho**2*(3*y**2-2*x**2) )
    t5=r**2*( 2*y**4 - gamma*(x**2 + z**2) )
    Byy=C*z/(2*alpha**4*beta**3*rho**4)* ( (t1-t2-t3)* E + (t4+t5)*alpha**2*K )
    
    Bzz=C*z/(2*alpha**4*beta**3)* ( (6*r0**2*(rho**2-z**2) -7*r0**4 + r**4 )* E + (r0**2 - r**2)*alpha**2*K )
    

    return (Bxx,Byy,Bzz)
#
#num I N R Z
def gen_coils(num_coils_in_group,current,coils,radi,zloc):
    I=ones(sum(num_coils_in_group))
    N=ones(sum(num_coils_in_group))
    R=ones(sum(num_coils_in_group))
    Z=ones(sum(num_coils_in_group))
    offset=0
    patch=0
    for n in num_coils_in_group:
        #i
        if len(current[patch]) == 1 or len(current[patch]) ==n:
            I[offset:offset+n]=current[patch]
        else:
            I[offset:offset+n]=np.linspace(current[patch][0],current[patch][1],n)
       #coils 
        if len(coils[patch]) == 1 or len(coils[patch]) ==n:
            N[offset:offset+n]=coils[patch]
        else:
            N[offset:offset+n]=np.linspace(coils[patch][0],coils[patch][1],n)
        #r
        if len(radi[patch]) == 1 or len(radi[patch]) ==n:
            R[offset:offset+n]=radi[patch]
        else:
            R[offset:offset+n]=np.linspace(radi[patch][0],radi[patch][1],n)
        #z
        if (len(zloc[patch]) == 1 or len(zloc[patch]) == n):
            Z[offset:offset+n]=zloc[patch]
        else:
            Z[offset:offset+n]=np.linspace(zloc[patch][0],zloc[patch][1],n)
        patch+=1
        offset+=n
    return (I,N,R,Z)
def saveBfield(Bx,By,Bz,dim,dl,U):
    #single coil test
    #Bxp,Byp,Bzp=siglecoil_field1(100,10,0.2,X,Y,Z-2);
    (Nx,Ny,Nz)=dim
    B=np.zeros((3,Nx,Ny,Nz))
    B[0,...]=Bx*(dl[0]*dl[2])
    B[1,...]=By*(dl[2]*dl[1])
    B[2,...]=Bz*(dl[0]*dl[1])
    Version = 0;
    Type = 7;
    Dim= 4;
    info=np.array([Version,Type,Dim])
    DimArray =np.array([3,Nx,Ny,Nz])
    B=B/U[0]['B'];
    Data=B.flatten(order='F');
    savefid = open("B_file",'wb');
    info.tofile(savefid)
    DimArray.tofile(savefid)
    Data.tofile(savefid)
    savefid.close()

def savefilter(name,Bx,By,Bz,dim):
    #single coil test
    #Bxp,Byp,Bzp=siglecoil_field1(100,10,0.2,X,Y,Z-2);
    (Nx,Ny,Nz)=dim
    B=np.zeros((3,Nx,Ny,Nz))
    B[0,...]=Bx
    B[1,...]=By
    B[2,...]=Bz
    Version = 0;
    Type = 7;
    Dim= 4;
    info=np.array([Version,Type,Dim])
    DimArray =np.array([3,Nx,Ny,Nz])
    Data=B.flatten(order='F');
    savefid = open(name,'wb');
    info.tofile(savefid)
    DimArray.tofile(savefid)
    Data.tofile(savefid)
    savefid.close()
def saveNei(name,Bx,By,dim):
    (Nx,Ny,Nz)=dim
    B=np.zeros((Nx,Ny,Nz,2))
    B[...,0]=Bx
    B[...,1]=By
    Version = 0;
    Type = 7;
    Dim= 4;
    info=np.array([Version,Type,Dim])
    DimArray =np.array([Nx,Ny,Nz,2])
    Data=B.flatten(order='F');
    savefid = open(name,'wb');
    info.tofile(savefid)
    DimArray.tofile(savefid)
    Data.tofile(savefid)
    savefid.close()

'''
#3
B[0,Nx:,:,:Nz]=np.flip(Bxp,axis=2)
B[2,Nx:,:,:Nz]=np.flip(Bzp,axis=2)
#2
B[0,:Nx,:,:Nz]=np.flip(Bxp,axis=(0,2))
B[2,:Nx,:,:Nz]=np.flip(Bzp,axis=(0,2))

#4
B[0,Nx:,:,Nz:]=Bxp
B[2,Nx:,:,Nz:]=Bzp

#1
B[0,:Nx,:,Nz:]=np.flip(Bxp,axis=0) 
B[2,:Nx,:,Nz:]=np.flip(Bzp,axis=0) 
'''

def tanhdist(a,b,c,d,x):
    return a*np.tanh(b*(x-c))+d
