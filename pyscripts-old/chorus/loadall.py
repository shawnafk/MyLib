import numpy as np
#input parameter
case=2
if case == 1: 
    #half
    dT = 89.714891963013827
    wl= 0.30328825168137874
    gamma_l = 6.6697249320641283E-003
elif case == 2:
    dT = 48.898864632961590
    wl = 0.30328825168137874/5
    gamma_l = 5.3357799456513026E-003
elif case == 3:
    #6
    dT = 34.447285418928061
    wl = 0.30728461554653536/5
    gamma_l = 7.7344273377194586E-003
elif case == 4:
    #7
    dT = 25.308209695538984
    wl = 0.30728461554653536/5
    gamma_l = 7.7344273377194586E-003
#fixed data
zpos=np.load("zpos.npy")
kmode=np.load("kmode.npy")
gyro=np.load("gyro.npy")
gm=np.load("gamma.npy")
wp2=np.load("wp2.npy")
vr = np.load("vr.npy")
vg = np.load("vg.npy")

#data A complex
chor= np.load("chor.npy")
sour= np.load("sour.npy")

