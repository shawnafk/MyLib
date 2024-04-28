#!/usr/bin/env python3
# coding=utf-8
import numpy as np
from numpy import *
nx=128
ny=128
theta=np.zeros((nx,ny))
#meshgrid which is which
[sy,sx]=meshgrid(arange(-int(nx/2),int(nx/2)),arange(-int(ny/2),int(ny/2)));
r=sqrt(sx**2+sy**2)
theta[:,int(ny/2):]=arccos(sx[:,int(ny/2):]/r[:,int(ny/2):])
theta[:,:int(ny/2)]=2*pi-arccos(sx[:,:int(ny/2)]/r[:,:int(ny/2)])
