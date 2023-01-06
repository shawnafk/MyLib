#!/opt/anaconda3/bin/python3
import scipy.io as sio
import numpy as np
def py2mat(dat,outfile):
    sio.savemat(outfile+'.mat',{'Dat':dat})

#PEVTK
'''
from pyevtk.hl import gridToVTK
# Dimensions
dx, dy, dz = 1e-3, 1e-3, 12e-3
Lx, Ly, Lz = nx*dx, ny*dy, nz*dz
ncells = nx * ny * nz
xs=-Lx/2
xe=Lx/2

ys=-Ly/2
ye=Ly/2

zs=-Lz/2
ze=Lz/2

npoints = (nx + 1) * (ny + 1) * (nz + 1)

X=np.linspace(xs,xe,nx);
Y=np.linspace(ys,ye,ny);
Z=np.linspace(zs,ze,nz);
y,x,z=meshgrid(X,Y,Z);

gridToVTK("./structured", x, y, z, cellData = {"pressure" : pressure}, pointData = {"temp" : temp})
'''
from tvtk.api import tvtk, write_data
# Unpack velocity information
def py2vtkv(dat,name,ext):
    nx,ny,nz=dat.shape
    # Generate the grid,points
    xx,yy,zz=np.mgrid[0:nx,0:ny,0:nz]
    pts = np.empty((nx,ny,nz,3), dtype=int)
    pts[..., 0] = xx
    pts[..., 1] = yy
    pts[..., 2] = zz
    # We reorder the points and vectors so this is as per VTK's
    # requirement of x first, y next and z last.
    pts = pts.transpose(2, 1, 0, 3).copy()
    pts.shape = pts.size // 3, 3
    vectors = dat.transpose(3, 2, 1, 0).copy()
    vectors.shape = vectors.size // 3, 3
    sg = tvtk.StructuredGrid(dimensions=xx.shape, points=pts)
    sg.point_data.vectors = vectors
    sg.point_data.vectors.name = name
    write_data(sg, name+ext)

def py2vtks(dat,name,ext):
     nx,ny,nz=dat.shape
     pts = np.empty((nx,ny,nz,3), dtype=int)
     xx,yy,zz=np.mgrid[0:nx,0:ny,0:nz]
     pts[..., 0] = xx
     pts[..., 1] = yy
     pts[..., 2] = zz
     pts = pts.transpose(2, 1, 0, 3).copy()
     pts.shape = pts.size // 3, 3
     sg = tvtk.StructuredGrid(dimensions=xx.shape, points=pts)
     scalars = dat.T.copy()
     sg.point_data.scalars = scalars.ravel()
     sg.point_data.scalars.name = name
     write_data(sg, name+ext)


