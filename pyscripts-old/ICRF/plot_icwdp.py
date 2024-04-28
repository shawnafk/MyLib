from icw import *
#Examples
den=1e13
mag=1e4
ic=1
im=1
mu = 1
z = 1
wci = Oci(im,ic,mag)
Va = Alfven_v(den,mag,im,ic)/c

rf=1.1
angles=np.linspace(0,2*pi,1000)
n2 = np.array([general_n2t(den, mag, im,ic,t,rf*wci) for t in angles])

#saw_n = SAW(den,mag,ic,im,angles,rf*wci)
#mswf_n,msws_n= MSW(0.3,den,mag,ic,im,angles,rf*wci)

#only correct at wci
icwf_n = fast(den,mag,ic,im,angles,rf*wci)
icws_n = icw(den,mag,ic,im,angles,rf*wci)



fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(angles,1/np.sqrt(n2[:,0]),'k')
ax.plot(angles,1/np.sqrt(n2[:,1]),'k')
#ax.plot(angles,1/saw_n,'r')
#ax.plot(angles,1/msws_n,'g')
#ax.plot(angles,1/mswf_n,'b')

#ax.plot(angles,1/msws_n,'b')
ax.plot(angles,1/icws_n,'r--')
ax.plot(angles,1/icwf_n,'g--')


plt.show()
plt.savefig(str(rf)+'.pdf',bbox_inches='tight')
fig.set_size_inches(12,8)

#icw_vert
rf= np.logspace(-0.5,2,1000) 
freqs=rf* wci
npara=10
nv_stix=np.array([icwf_stix_nv(den, mag,im,ic,npara,f) for f in freqs])
nv_full=(cplx_filter(np.array([general_n2vp(den, mag, im,ic,npara,f) for f in freqs])))
 
fig, ax = plt.subplots()
#ax.plot(np.sqrt(nv_full[:,0])/c*freqs,rf,'k-.')
ax.plot(np.sqrt(nv_full[:,1])/c*freqs,rf,'g-')
ax.plot(np.sqrt(nv_stix)/c*freqs,rf,'r-.')
va=Alfven_v(den,mag,im,ic)
ax.plot(freqs/va,rf,'k')

rds=np.linspace(0.1,10,1000)
nv_full=(cplx_filter(np.array([general_n2vp(rd*den,mag,im,ic,npara,2*wci) for rd in rds])))