import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import scipy.special as sp

Nx=5000
Ny=1;
Nz=5000
#此处确定图像显示范围的参数，各坐标轴单位均为1m
def siglecoil_field1(I0,N_coil,r0,x,y,z):
    c0=4e-7; # mu0/pi=4e-7
    I0=I0*N_coil;
    # current, unit: A
    C=c0*I0;
    # Set field for single coil
    rho=(x**2+y**2)**0.5
    r=(x**2+y**2+z**2)**0.5
    alpha=(r0**2 + r**2 - 2*r0*rho)**0.5
    beta=(r0**2 + r**2 + 2*r0*rho)**0.5
    k2 = 1 - alpha**2/beta**2
    (K,E)=(sp.ellipk(k2),sp.ellipe(k2))
    Bx=C*x*z/(2*alpha**2*beta*rho**2) * ( (r0**2+r**2)*E - alpha**2 * K)
    By=C*y*z/(2*alpha**2*beta*rho**2) * ( (r0**2+r**2)*E - alpha**2 * K)
    Bz=C/(2*alpha**2*beta) * ((r0**2 - r**2)*E + alpha**2 * K)
    return (Bx,By,Bz)

dl=np.array((1,1,1))#对修改单位影响巨大，应该是步长一类的东西
#print(dl)
deltax=1e-3
dx,dy,dz=deltax*dl
x=np.linspace(-Nx/2,Nx/2,Nx)*dx
z=np.linspace(-Nz/2,Nz/2,Nz)*dz
#print(x.shape,z.shape)
#print(x,z)
Z,X=np.meshgrid(z,x)
#print(X.shape,Z.shape)
#print(X,Z)

I=1
B=np.zeros((3,Nx,1,Nz))
A=np.zeros((3,3,1,2))#生成三行三列的分块矩阵，矩阵中的每个元素是一行二列的小矩阵
N=1
R0=1
#print(B.shape)
#print(A)
Bxr,Byr,Bzr=siglecoil_field1(I,N,R0,X,0,Z)
#print(Byr.shape)
B[0,:,0,:]+=(Bxr)
B[1,:,0,:]+=(Byr)
B[2,:,0,:]+=(Bzr)
#print(__z,B.shape)
#print(B)
B[0,...]*=(dx*dz)/deltax**2
B[1,...]*=(dz*dy)/deltax**2
B[2,...]*=(dx*dy)/deltax**2
[sx,sz]=np.meshgrid(np.arange(0,Nz),np.arange(0,Nx));
plt.figure(1)
plt.title('Siglecoil_Field')
plt.xlabel('Z/m')
plt.ylabel('X/m')
xticks=[0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000]
xticklabels=['-2.5','-2','-1.5','-1','-0.5','0','0.5','1','1.5','2','2.5']
plt.xticks(xticks,xticklabels)

yticks=[0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000]
yticklabels=['-2.5','-2','-1.5','-1','-0.5','0','0.5','1','1.5','2','2.5']
plt.yticks(yticks,yticklabels)

plt.streamplot(sx, sz, B[2,:,0,:], B[0,:,0,:], density=[2,2])
#plt.savefig('Siglecoil_Field.png')
plt.show()
