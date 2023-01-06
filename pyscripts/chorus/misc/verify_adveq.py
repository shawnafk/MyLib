import numpy as np
zpos=np.load("../data/T_zpos.npy")
vg = np.loadtxt("../data/ori/vg.out")

kl=np.loadtxt("../data/ori/kmode.out")
knl = np.load("../data/C_knl.npy")
#JB/JE
ratio=np.load("../data/C_ratio.npy")
chor= np.load("../data/T_chor.npy")
amp = np.linalg.norm(chor,axis=-1)
dt = 49.6
ds = zpos[np.newaxis,1:] - zpos[np.newaxis,:-1]
dampdt = (amp[1:,:] - amp[:-1,:])/dt
dampds = (amp[:,1:-1] - amp[:,:-2])/ds

#time cut 1, space cut 2
real_cmp = 1/amp[1:,1:-1] * (dampdt[:,1:-1] + vg[np.newaxis,1:]*dampds[1:,:])


#imag_cmp = vg[np.newaxis,:] * (knl**2 - kl[np.newaxis,:]**2)/kl[np.newaxis,:]/2
#left = real_cmp
#right = imag_cmp/ratio[:-1,:-1]


imag_cmp = vg[np.newaxis,:] * (knl - kl[np.newaxis,:])
left = real_cmp/imag_cmp[:,:-1]
right =-1/ratio[:-1,:-1]


from fig_head import *
fig,axes = plt.subplots(2,2)
t = [[200,300],[400,500]]
ord = [['(a)','(b)'],['(c)','(d)']]
lgd_loc = [['upper right','upper right'],['lower left','lower left']]
zrange=np.arange(500,900)
for i in range(2):
    for j in range(2):
        axes[i][j].plot(zpos[zrange],left[t[i][j],zrange],lw=1.5,color='r',label='left')
        axes[i][j].plot(zpos[zrange],right[t[i][j],zrange],lw=1.5,color='k',label='right')
        #axes[i][j].legend()
        #axes[i][j].legend(loc=lgd_loc[i][j])
        #axes[i][j].set_xlim(zrange[0],zrange[-1])
        #axes[i][j].ticklabel_format(axis='y', style='scientific', scilimits=(0,1))
        plt.text(0.1, 0.9, ord[i][j],fontsize=18,horizontalalignment='center',verticalalignment='center',transform=axes[i][j].transAxes)
plt.subplots_adjust(hspace=0.25)
fig.text(0.5, 0.04, r"$Z/ \frac{\rm{c}}{\omega_{\rm{p0}}}$", ha='center',font=labelfont)
fig.text(0.04, 0.5, 'Amplitude', va='center', rotation='vertical',font=labelfont)
fig.set_size_inches(12,8)
plt.show()
plt.savefig('adv_term.pdf',bbox_inches='tight')

