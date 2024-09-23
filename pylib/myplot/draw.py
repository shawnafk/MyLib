#!/usr/bin/env python3
# coding=utf-8
#regradless our normalzation unit, this DT below is normalized by 1/wpe, Dx is normalized with L inertial
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
def gen_w(DT,N):
	w_max = 2*np.pi/2/DT
	w = np.linspace(0,w_max,int(N/2))
	return w

def gen_k(DX,N):
	k_max = 2*np.pi/2/DX
	k = np.linspace(0,k_max,int(N/2))
	return k
#maybe high dimensional

def cartoon1d_low(x,yt_s,Vmin=0,Vmax=0,dt=1,units='s',*args,**kwargs):
	fig=plt.figure()
	ax=plt.axes()
	time = ax.annotate(0,xy=(0, 0),xycoords='figure fraction')
	enter=0
	#draw first
	def animate(i):
		#if set ylim	
		if (Vmin!=Vmax and enter==0):
			enter=1
			ax.set_ylim(Vmin,Vmax)
		ax.clear()
		text='Frames = '+str(i) + 'Time =' +str( round(i*dt,2)) + units
		for yt in yt_s:
			ax.plot(x,yt[:,i])
		ax.annotate(text, xy=(0.85, 0.95), xycoords='axes fraction')
	anim = animation.FuncAnimation(fig, animate, *args,**kwargs)
	return anim

def cartoon1d_fixed(x,yt_s,Vmin,Vmax,*args,**kwargs):
    fig=plt.figure()
    ax=plt.axes()
    time = ax.annotate(0,xy=(0.2, 0.9),xycoords='figure fraction')
    #draw first
    lines = []
    lt=0
    for yt in yt_s:
        lines.append(ax.plot(x, yt[:, 0],lw=2)[0])
        lt=lt+1
    def animate(i):
        text='Frames = '+str(i)
        time.set_text(text)
        l_th = 0
        for yt in yt_s:
            lines[l_th].set_ydata(yt[:,i])
            l_th=l_th+1
        ax.set_ylim(Vmin,Vmax)
    anim = animation.FuncAnimation(fig, animate, *args,**kwargs)
    return anim

def cartoon1d_dyn(x,yt_s,*args,**kwargs):
    fig=plt.figure()
    ax=plt.axes()
    time = ax.annotate(0,xy=(0.2, 0.9),xycoords='figure fraction')
    #draw first
    lines = []
    lt=0
    for yt in yt_s:
        lines.append(ax.plot(x, yt[:, 0],lw=2)[0])
        lt=lt+1
    def animate(i):
        text='Frames = '+str(i)
        time.set_text(text)
        l_th = 0
        y_min=[]
        y_max=[]
        for yt in yt_s:
            y_min.append(np.min(yt[:,i]))
            y_max.append(np.max(yt[:,i]))
        ymin,ymax=np.min(y_min),np.max(y_max)
        yrange = ymax - ymin
        ypadding = yrange * 0.1 # add 10% padding to the range
        for yt in yt_s:
            lines[l_th].set_ydata(yt[:,i])
            l_th=l_th+1
        ax.set_ylim(ymin - ypadding, ymax + ypadding)
    anim = animation.FuncAnimation(fig, animate, *args,**kwargs)
    return anim

def cartoon2d_dyn(xyextent,array,*args,**kwargs):
	fig=plt.figure()
	ax=plt.axes()
	time = ax.annotate(0,xy=(0.2, 0.9),xycoords='figure fraction')
	im = ax.imshow(array[...,0],aspect='auto',extent=xyextent,cmap='jet',interpolation='lanczos')
	cbar = fig.colorbar(im)
	def animate(i):
		Z=array[:,:,i]
		im.set_array(Z)
		im_min, im_max = np.min(Z), np.max(Z)
		im.set_clim(im_min, im_max)
		time.set_text('Frames = %d' %i)
	anim = animation.FuncAnimation(fig, animate, *args,**kwargs)
	return anim

def dispersion(ax,data,k,w,level=10):
	#def dispersion(ax,data,k,w):
	grids=data.shape[0]
	steps=data.shape[1]
	X,Y = np.meshgrid(k,w)
	fdata=np.fft.fftshift(np.fft.fft2(data))
	#even
	if grids % 2 == 0:
		if steps % 2 == 0:
			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2)-1:grids-1,int(steps/2)-1:steps-1])),level)
#			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2)-1:grids-1,int(steps/2)-1:steps-1])))
		else:
			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2)-1:grids-1,int(steps/2):steps-1])),level)
#			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2)-1:grids-1,int(steps/2):steps-1])))
	else:	
		if steps % 2 == 0:
			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2):grids-1,int(steps/2)-1:steps-1])),level)
#			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2):grids-1,int(steps/2)-1:steps-1])))
		else:
			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2):grids-1,int(steps/2):steps-1])),level)
#			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2):grids-1,int(steps/2):steps-1])))

	return 0

def dispersionf(ax,data,k,w,level=10):
	#def dispersion(ax,data,k,w):
	grids=data.shape[0]
	steps=data.shape[1]
	X,Y = np.meshgrid(k,w)
	fdata=np.fft.fftshift(np.fft.fft2(data))
	#even
	if grids % 2 == 0:
		if steps % 2 == 0:
			ax.contourf(X,Y,np.transpose(abs(fdata[int(grids/2)-1:grids-1,int(steps/2)-1:steps-1])),level)
#			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2)-1:grids-1,int(steps/2)-1:steps-1])))
		else:
			ax.contourf(X,Y,np.transpose(abs(fdata[int(grids/2)-1:grids-1,int(steps/2):steps-1])),level)
#			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2)-1:grids-1,int(steps/2):steps-1])))
	else:	
		if steps % 2 == 0:
			ax.contourf(X,Y,np.transpose(abs(fdata[int(grids/2):grids-1,int(steps/2)-1:steps-1])),level)
#			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2):grids-1,int(steps/2)-1:steps-1])))
		else:
			ax.contourf(X,Y,np.transpose(abs(fdata[int(grids/2):grids-1,int(steps/2):steps-1])),level)
#			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2):grids-1,int(steps/2):steps-1])))

	return 0



#c = 1.0;
#vth = 0.00442;
#x1 = np.linspace(0,k_max/4.,20)
#x2 = np.linspace(0,k_max,80)
#w_em = np.sqrt(1+1.5*x1**2*c**2)
#w_l = np.sqrt(1+1.5*x2**2*vth**2)
#ax.plot(x1,w_em,'r.',ms=1.2)
#ax.plot(x2,w_l,'k+',ms=1.2)
#w_ce = 1./(3.21*10**-3*10**4*10**-1)
#ax.plot((0,k_max),(w_ce,w_ce),'k-')
#ax.plot((0,w_max),(0,w_max),'r-')
#ax.plot(x2,w_l,'g.',ms=1.2)
#w right cyclotron
#w_Rcutoff=w_ce/2 + np.sqrt(1+w_ce**2/4)
#w_Lcutoff=-w_ce/2 + np.sqrt(1+w_ce**2/4)

#w_rc1 = np.linspace(0.001,w_ce-0.001,200)
#w_rc2 = np.linspace(w_Rcutoff,w_max,200)
#k_rc1 = np.sqrt(w_rc1**2 - 1/(1-w_ce/w_rc1));
#k_rc2 = np.sqrt(w_rc2**2 - 1/(1-w_ce/w_rc2));

#ax.plot(k_rc1,w_rc1,'g-',alpha=0.7)
#ax.plot(k_rc2,w_rc2,'g-')

#w_lc = np.linspace(w_Lcutoff,w_max,200)
#k_lc = np.sqrt(w_lc**2 - 1/(1+w_ce/w_lc));

#ax.plot(k_lc,w_lc,'k-',alpha=0.7)

def powerspectrum(ftEk,omega,dt,steps):
	time=np.arange(0,dt*steps,dt)
	P=abs(average(exp(-1j*omega*time)*ftEk[...,:steps],axis=-1))**2
	return P

#Sy=-Ex*Bz
def ref(E,B,Rz,X,N):
	Sr=np.cross(E[:,X,0,Rz[0]:Rz[1],:],B[:,X,0,Rz[0]:Rz[1],:],axis=0)
	Sa=(average(Sr[0,:,:],axis=(0)))
	Weight=np.ones(N)/N
	swr=np.convolve(Weight,S_a)
	return swr
