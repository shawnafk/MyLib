#base data
def readchor(file,length):
     chor = np.loadtxt(file)
     newc = chor.reshape([int(chor.shape[0]/length),length,2])
     return newc
#phase space
dat=np.loadtxt("../data/ori/dist.out")
phase_space=dat.reshape([7,1001,31,401])
np.save("../data/T_phase.npy",phase_space)
