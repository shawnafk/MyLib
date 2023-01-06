#!/usr/bin/env python3
# coding=utf-8
from GAPS_IO import GAPS_IO as gio 
import numpy as np
def powerspectrum(ftEk,omega,dt,steps):
	time=np.arange(0,dt*steps,dt)
	P=abs(np.average(exp(-1j*omega*time)*ftEk[...,:steps],axis=-1))**2
	return P

#Sy=-Ex*Bz
#return S
#time average is equivlant to X average.
def ref(E,B,Rz,X,N):
	Sr=np.cross(E[:,X,0,Rz[0]:Rz[1],:],B[:,X,0,Rz[0]:Rz[1],:],axis=0)
	Sa=(np.average(Sr[0,:,:],axis=(0)))
	Weight=np.ones(N)/N
	swr=np.convolve(Weight,Sa)
	return swr

def ref2(E,B,X):
        Sr=np.cross(E[:,X,0,:,:],B[:,X,0,:,:],axis=0)
	Sa=(np.average(Sr[0,:,:],axis=(0)))
	return Sa
path=['1/','2/','3/','4/','5/']
num=5;
tmpE=[]
tmpB=[]
for _ in path:
	EB=gio.gload(_+'EB');
	tmpE.append(EB['Data'][...,0::2])
	tmpB.append(EB['Data'][...,1::2])

N=20
Loc=81
p=[]
for i in range(num):
	p.append(ref(tmpE[i],tmpB[i],(0,512),Loc,N))
np.save('no_gap.npy',p)


