#!/usr/bin/env python3
# coding=utf-8
# load using np.fromfile instance, can partially load some data
# gmap using np.memmap instance, more flexible, need more work to substract data we need.

import numpy as np


# maintain a dict
def GAPS_IO_GetType(inp):
    if inp == 0:
        y = 'int32'
        Bs = 4
    elif inp == 1:
        y = 'int64'
        Bs = 8
    elif inp == 2:
        y = 'int16'
        Bs = 2
    elif inp == 3:
        y = 'uint32'
        Bs = 4
    elif inp == 4:
        y = 'uint64'
        Bs = 8
    elif inp == 5:
        y = 'uint16'
        Bs = 2
    elif inp == 6:
        y = 'float32'
        Bs = 4
    elif inp == 7:
        y = 'float64'
        Bs = 8
    elif inp == 8:
        y = 'd'
        Bs = 8
    elif inp == 9:
        y = 'd'
        Bs = 8
    elif inp == 10:
        y = 'd'
        Bs = 8
    elif inp == 11:
        y = 'uint8'
        Bs = 1
    elif inp == 12:
        y = 'uint8'
        Bs = 1
    else:
        print('Input must be a integar from 0 to 12')
    return (y, Bs)


def load_info(filename):
    default_type = 'int64'
    fileID = open(filename, 'rb')
    # the first three int64 in GAPS_IO are version Type and DIM, Version can be 0 1 2, type means data type double, int, float etc, DIM is (data on grid) plus grid dimension.
    (Version, Type, Dim) = np.fromfile(fileID, dtype=default_type, count=3)
    DimArray = np.fromfile(fileID, dtype=default_type, count=Dim)
    numperstep = 1
    for i in DimArray:
        numperstep = numperstep * i
    position = [fileID.tell()]
    # 0 means abslutely from the begainning, 1 means from the first place, 2 means eof,move pointer to eof
    fileID.seek(0, 2)
    position.append(fileID.tell())
    # using position difference to calculate numsteps
    (precision, Bs) = GAPS_IO_GetType(Type)
    NumSteps = int((position[1] - position[0]) / (numperstep * Bs))
    # add time steps
    DimArray = np.append(DimArray, NumSteps)
    y = {};
    y['Version'] = Version
    y['Type'] = Type
    y['Dim'] = Dim
    y['DimArray'] = DimArray
    y['NumPerStep'] = numperstep
    y['Precision'] = precision
    y['Bs'] = Bs

    fileID.close()

    return y


# give filename and fileinfo
def seek_data(filename, y, ith):
    fileID = open(filename, 'rb')
    intp = 8;
    # data start
    fileID.seek(intp * (y['Dim'] + 3), 0)
    # offset
    numperstep = y['NumPerStep']
    Bs = y['Bs']
    fileID.seek(numperstep * Bs * ith, 1)
    return fileID


def gload(filename):
    # GAPS_IO_Load Load data from GAPS-IO standard file into workspace
    #   S = GAPS_IO_Load(FILENAME) loads the data from a GAPS-IO standard file
    #   named FILENAME into a struture S
    #   Data in S:
    #       S.Version: version information of the GAPS-IO file.
    #       S.Type: data type index of GAPS-IO file.
    #       S.Dim: dimentions of data.
    #       S.DimArray: a row array with S.Dim+1 elements, which gives
    #       numbers of data in each dimention. The last element of
    #       S.DimArray is number of time steps.
    #       S.NumPerStep: number of data per step.
    #       S.Data: a (S.Dim+1)-dimentional data

    # filename=varargin{1}

    # default_type='int64'
    y = load_info(filename)
    # move to data
    fileID = seek_data(filename, y, 0)
    NumSteps = y['DimArray'][-1];
    numperstep = y['NumPerStep']
    n = NumSteps * numperstep
    precision = y['Precision']
    DimArray = y['DimArray']
    #   #fit matlab fortran type
    Data = np.reshape(np.fromfile(fileID, dtype=precision, count=n), DimArray, order="F")
    y['Data'] = Data
    fileID.close()
    return y


def pload(filename, ith, NumSteps):
    # start and end are int type, i-th steps
    # default_type='int64'
    y = load_info(filename)
    # move to data
    fileID = seek_data(filename, y, ith)
    numperstep = y['NumPerStep']
    n = NumSteps * numperstep
    #   #fit matlab fortran type
    precision = y['Precision']
    DimArray = y['DimArray']
    DimArray[-1] = NumSteps
    Data = np.reshape(np.fromfile(fileID, dtype=precision, count=n), DimArray, order="F")
    y['Data'] = Data
    fileID.close()
    return y


# pload and gload canbe merged to one function
# careful that manipulate(cpu access) on memmaped file may cause page missing and thus very high latency. sys will consistently access disk, cause it use disk as cache.
# our solution is load all chunk into memory, so that when cpu load that many small part of data, although it is still insequential, this part of cache missing occurs at DRAM level which is much faster than disk io.
# conclusion, memmap should not access insequential data
def gmap(filename):
    y = load_info(filename)
    dtype = y['Precision']
    DimArray = y['DimArray']
    intp = 8;
	# offset skip date info segments
    offset = intp * (y['Dim'] + 3)
    Data = np.memmap(filename, dtype, 'r', offset, order='F', shape=tuple(DimArray))
    return Data


# For EM data
def map_E(filename):
    EB = gmap(filename)
    E = EB[..., 0::2]
    return E


def map_B(filename):
    EB = gmap(filename)
    B = EB[..., 1::2]
    return B


# for EN
def map_EK(filename):
    EN = gmap(filename)
    EK = EN[0:3, ...]
    return EK


# the 0th components storages all species particle infomation, the data are organized like
# 1 double: Ns (number of Sth species) + grid_cache_len * 6 double : sth species partilce infomation
# s species, int
# xb yb zb boundray on each directions,tuple
def load_particle(dat, grid_cache_len, s, xb, yb, zb):
    # ns storage particle num of a cell (ijk).In dat, each cell (:::) has (1+grid_cache_len)*s elements
    # within first ele denotes num particle of this species
    ns = dat[(grid_cache_len * 6 + 1) * s, :, :, :, 0]
    tmp = [];
    # ergodic all grid and extract thir
    for i in range(xb[0], xb[1]):
        for j in range(yb[0], yb[1]):
            for k in range(zb[0], zb[1]):
                if ns[i, j, k] == 0:
                    continue
                else:
                    particle_start = (grid_cache_len * 6 + 1) * s + 1
                    particle_end = particle_start + int(ns[i, j, k]) * 6
                    # python slice will choose particle_start to particle_end - 1 which are int_ns particles
                    # tmp.append( dat[1:int(ne[i,j,k])*6+1,i,j,k,0].reshape((int(ne[i,j,k]),6)))
                    tmp.append(dat[particle_start:particle_end, i, j, k, 0].reshape((int(ns[i, j, k]), 6)))

    Part = np.concatenate((tmp))
    return Part


def load_tail(dat, grid_cache_len, s,vcrt):
    # ns storage particle num of a cell (ijk).In dat, each cell (:::) has (1+grid_cache_len)*s elements
    # within first ele denotes num particle of this species
    # Jx_array=np.zeros(dat.shape[1:-1])
    # Jy_array=np.zeros(dat.shape[1:-1])
    Jz_array = np.zeros(dat.shape[1:-1])
    # ergodic all grid and extract thir
    for spec in range(s):
        print(spec)
        ns = dat[(grid_cache_len * 6 + 1) * spec, :, :, :, 0]
        print(dat.shape)
        I = dat.shape[1] - 1 if dat.shape[1] != 1 else 1
        J = dat.shape[2] - 1 if dat.shape[2] != 1 else 1
        K = dat.shape[3] - 1 if dat.shape[3] != 1 else 1
        for i in range(I):
            for j in range(J):
                for k in range(K):
                    if ns[i, j, k] == 0:
                        continue
            else:
                for p in range(int(ns[i, j, k])):
                    p_offset = (grid_cache_len * 6 + 1) * spec + 1
                    xx = dat[p_offset + 6 * p + 0, i, j, k, 0] % 1.
                    xy = dat[p_offset + 6 * p + 1, i, j, k, 0] % 1.
                    xz = dat[p_offset + 6 * p + 2, i, j, k, 0] % 1.
                    vx = dat[p_offset + 6 * p + 3, i, j, k, 0]
                    vy = dat[p_offset + 6 * p + 4, i, j, k, 0]
                    vz = dat[p_offset + 6 * p + 5, i, j, k, 0]
                    if (abs(vz) > vcrt):
                        s1 = (1 - xx) * (1 - xy)
                        s2 = xx * (1 - xy)
                        s3 = xx * xy
                        s4 = xy * (1 - xx)
                        Jz_array[i, j, k] += s1 * vz
                        Jz_array[i + 1, j, k] += s2 * vz
                        Jz_array[i, j, k + 1] += s4 * vz
                        Jz_array[i + 1, j, k + 1] += s3 * vz
    return Jz_array

