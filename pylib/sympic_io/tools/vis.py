#!/usr/bin/env python3
# coding=utf-8
#regradless our normalzation unit, this DT below is normalized by 1/wpe, Dx is normalized with L inertial
from numpy.fft import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation


def gen_w(DT, N):
	w_max = 2*np.pi/2/DT
	w = np.linspace(0, w_max, int(N/2))
	return w


def gen_k(DX, N):
	k_max = 2*np.pi/2/DX
	k = np.linspace(0, k_max, int(N/2))
	return k
#maybe high dimensional


def cartoon_2l(array, num, Vmin, Vmax, dt=0.5, *args, **kwargs):
	fig = plt.figure()
	ax = plt.axes()
	time = ax.annotate(0, xy=(0, 0), xycoords='figure fraction')
	ax.set_ylim(Vmin, Vmax)
	#draw first
	lines = [ax.plot(array[_][:, 0])[0] for _ in range(num)]

	def animate(i):
		for _, line in enumerate(lines):
			line.set_ydata(array[_][:, i])
		time.set_text('Frames = %d,time = %0.2f ns' % (i, i*dt))
		#using function pointer
	anim = animation.FuncAnimation(fig, animate, *args, **kwargs)
	return anim

#fig and ax, data, dt,dim, and animation args


def cartoon(array, Vmin=0, Vmax=0, dt=0.5, flag='1d', *args, **kwargs):
	fig = plt.figure()
	ax = plt.axes()
	time = ax.annotate(0, xy=(0, 0), xycoords='figure fraction')
	if '1d' == flag:
		ax.set_ylim(Vmin, Vmax)
		#draw first
		line = []
		line = ax.plot(array[:, 0], color='k', lw=2)[0]

		def animate(i):
			line.set_ydata(array[:, i])
			time.set_text('Frames = %d,time = %0.2f ns' % (i, i*dt))
			#using function pointer
	elif '2d' == flag:
		#draw first
		#		cax = ax.pcolormesh(array[:-1, :-1, 0], vmin=-1, vmax=1, cmap='jet')
		#cax = ax.pcolormesh(array[:, :, 0], cmap='jet ',shading='gouraud')
		#cax = ax.pcolormesh(array[:, :, 0],vmin=-Vrange,vmax=Vrange, shading='gouraud')
		im = plt.imshow(array[:, :, 0], vmin=Vmin, vmax=Vmax,
							interpolation='lanczos', aspect='auto')
		fig.colorbar(im)
		def animate(i):
		#im.set_array(array[:,:,i].flatten())
			im.set_array(array[:, :, i])
			#plt.imshow(array[:, :, i],aspect='auto')
			time.set_text('Frames = %d,time = %0.2f ns' % (i, i*dt))
	anim = animation.FuncAnimation(fig, animate, *args, **kwargs)
	return anim

def cartoon1d(array,*args,**kwargs):
	fig = plt.figure()
	ax = plt.axes()
	time = ax.annotate(0, xy=(0, 0), xycoords='figure fraction')
	#draw first
	line = []
	line = ax.plot(array[:, 0], color='k', lw=2)[0]
	def animate(i):
		y=array[:, i]
		line.set_ydata(y)
		ax.set_ylim(min(y),max(y))
		time.set_text('Frames = %d'%i)
		#using function pointer
	anim = animation.FuncAnimation(fig, animate, *args, **kwargs)
	return anim

def dispersion(ax, data, k, w, level=10):
	#def dispersion(ax,data,k,w):
	grids = data.shape[0]
	steps = data.shape[1]
	X, Y = np.meshgrid(k, w)
	fdata = np.fft.fftshift(np.fft.fft2(data))
	#even
	if grids % 2 == 0:
		if steps % 2 == 0:
			ax.contour(X, Y, np.transpose(
				abs(fdata[int(grids/2)-1:grids-1, int(steps/2)-1:steps-1])), level)
#			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2)-1:grids-1,int(steps/2)-1:steps-1])))
		else:
			ax.contour(X, Y, np.transpose(
				abs(fdata[int(grids/2)-1:grids-1, int(steps/2):steps-1])), level)
#			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2)-1:grids-1,int(steps/2):steps-1])))
	else:
		if steps % 2 == 0:
			ax.contour(X, Y, np.transpose(
				abs(fdata[int(grids/2):grids-1, int(steps/2)-1:steps-1])), level)
#			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2):grids-1,int(steps/2)-1:steps-1])))
		else:
			ax.contour(X, Y, np.transpose(
				abs(fdata[int(grids/2):grids-1, int(steps/2):steps-1])), level)
#			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2):grids-1,int(steps/2):steps-1])))

	return 0


def dispersionf(ax, data, k, w, level=10):
	#def dispersion(ax,data,k,w):
	grids = data.shape[0]
	steps = data.shape[1]
	X, Y = np.meshgrid(k, w)
	fdata = np.fft.fftshift(np.fft.fft2(data))
	#even
	if grids % 2 == 0:
		if steps % 2 == 0:
			ax.contourf(X, Y, np.transpose(
				abs(fdata[int(grids/2)-1:grids-1, int(steps/2)-1:steps-1])), level)
#			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2)-1:grids-1,int(steps/2)-1:steps-1])))
		else:
			ax.contourf(X, Y, np.transpose(
				abs(fdata[int(grids/2)-1:grids-1, int(steps/2):steps-1])), level)
#			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2)-1:grids-1,int(steps/2):steps-1])))
	else:
		if steps % 2 == 0:
			ax.contourf(X, Y, np.transpose(
				abs(fdata[int(grids/2):grids-1, int(steps/2)-1:steps-1])), level)
#			ax.contour(X,Y,np.transpose(abs(fdata[int(grids/2):grids-1,int(steps/2)-1:steps-1])))
		else:
			ax.contourf(X, Y, np.transpose(
				abs(fdata[int(grids/2):grids-1, int(steps/2):steps-1])), level)
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

def powerspectrum(ftEk, omega, dt, steps):
	time = np.arange(0, dt*steps, dt)
	P = abs(average(exp(-1j*omega*time)*ftEk[..., :steps], axis=-1))**2
	return P

#Sy=-Ex*Bz


def ref(E, B, Rz, X, N):
	Sr = np.cross(E[:, X, 0, Rz[0]:Rz[1], :], B[:, X, 0, Rz[0]:Rz[1], :], axis=0)
	Sa = (average(Sr[0, :, :], axis=(0)))
	Weight = np.ones(N)/N
	swr = np.convolve(Weight, S_a)
	return swr


def intE(E, ax, start):
	subE = np.split(E, [start], axis=ax)
	pl = -np.cumsum(np.flip(subE[0], axis=ax), axis=ax)
	pl = -np.flip(pl, axis=ax)
	pr = -np.cumsum(subE[1], axis=ax)
	return np.concatenate((pl, pr), axis=ax)


def poisson2d(q, gx, gz):
	ax = [0, 1]
	x = np.linspace(-np.pi, np.pi, gx)
	z = np.linspace(-np.pi, np.pi, gz)
	kx, kz = np.meshgrid(z, x, sparse=True)
	ksq = (kx**2 + kz**2)
	ftq = fftn(q, axes=ax)
	phi = np.fft.fftn(
		fftshift((fftshift(ftq, axes=ax)/-ksq[..., np.newaxis]), axes=ax), axes=ax)
	return np.real(phi)
