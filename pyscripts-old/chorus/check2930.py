def ave_d1(array):
    return (array[1:,...] + array[:-1,...])/2

def ave_d2(array):
    return (array[:,1:,...] + array[:,:-1,...])/2
chor=np.load("chor.npy")

ampchor=np.linalg.norm(chor,axis=-1)
dT = 48

papt = (ampchor[1:,:]-ampchor[:-1,:])/dT
zpos=np.load("zpos.npy")
dS = zpos[1:] - zpos[:-1]
da_s = (ampchor[:,1:]-ampchor[:,:-1])
paps = ave_d2(da_s)/dS

#simulation apply the cold group velocity, so am I
vg=np.load("vg.npy")

knl = np.load("./data/out_ksim.npy")
kl = np.load("./data/out_klin.npy")

#match dimensions
right1 = (papt[:,1:-1] + paps[1:,:] * vg[1:])
right2 = (knl**2 - kl[:,1:]**2) * ampchor[1:,1:-1]
source = np.load("sour.npy")

left1 = source[1:,1:,0]*vg/2/kl[:,:]
left2 = source[:,:,1]


L = (papt[:,1:-1] + paps[1:,:] * vg[1:])/ampchor[1:,1:-1]
R = (knl**2 - kl[:,1:]**2) *vg[1:]/2/kl[:,1:]  * source[1:,1:-1,0] / source[1:,1:-1,1] 