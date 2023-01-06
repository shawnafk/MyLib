'''
amp = np.linalg.norm(chor[:-1,1:-1,...],axis=-1)
wbsq = np.sqrt(np.sqrt(2*wc*np.average(J+pcr,axis=-1))[np.newaxis,:-1]*kall**2)*amp
np.save("../data/C_wbsp.npy",wbsq)
def getfreq(freq,dT,trange):
    from PyEMD import EMD
    emd = EMD()
    imfs=emd(freq)
    dfdt = (imfs[-1,1:] - imfs[-1,:-1])/dT
    ave=np.average(dfdt[trange])  
    return ave

coef = (1 - vr/vg)**-2
fcr_c=[]
zloc = np.arange(520,600,10)
for loc in zloc:
    fcr_c.append(getfreq(wall[np.arange(200,600),loc],dT,np.arange(250,350)))

R1 = np.array(fcr_c)/(coef[zloc]*wbsq[500,zloc])    
R0=[]
trange=np.arange(450,600)
for loc in zloc:
    R0.append(cl.calcR(ratio[trange,loc]))
fcr_t = []
fcr_ave_spt = []
fcr_ave = []
for _R,loc in zip(R0,zloc):
    fcr = _R * coef[loc]*wbsq[trange,loc]
    fcr_t.append(fcr)
    fcr_ave.append(np.average(fcr))
    fcr_ave_spt.append(np.average(_R) * coef[loc]*np.average(wbsq[trange,loc]))
#\frac{\partial \omega}{\partial t}=\frac{1}{2}\left(1-\frac{v_{r}}{v_{g}}\right)^{-2} \omega_{t r}^{2}

def D_alpha(wbsq,fcr,J,k,dgds):
    return 4/wbsq*fcr + J*k/wbsq*dgds
dgds=(gyro[1:]-gyro[:-1])/(zpos[1:]-zpos[:-1])
Alpha=[]
for i,loc in zip(range(8),zloc):
    wbsqave = np.average(wbsq[trange,loc])
    Jave=np.average(J[loc,:])
    kave=np.average(kall[trange,loc])
    Alpha.append(D_alpha(wbsqave,fcr_c[i],Jave,kave,dgds[loc]))
'''
