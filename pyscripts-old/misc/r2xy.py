import numpy as np
r0=256
Nr=256
Nt=512
Nz=640
Tr=np.zeros(Nr,Nz)
idr=0
Tiv=np.load("Tiv.npy")
for r in range(1,Nr):
    x=[r+1]
    y=[0]
    for theta in np.linspace(0,2*np.pi,Nt):
        rx=np.floor(r*np.cos(theta))
        ry=np.floor(r*np.sin(theta))
        count=0
        if(rx==x[-1] and ry==y[-1]):
            print("same grid, skip")
        else:
            sumT += Tiv[rx,ry,:]
            count+=1
    Tr[r,:] = sumT/count
