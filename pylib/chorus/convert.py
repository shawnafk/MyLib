import numpy as np
def _mid(arr):
    return (arr[1:] + arr[:-1])/2
def mid(s,ax):
    return np.apply_along_axis(_mid,ax,s)
def comp(s):
    return s[...,0]+s[...,1]*1j

def g_phase(s):
#should be - or +
    return np.angle(s)
