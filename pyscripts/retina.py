import numpy as np
#assume human eye angle limit is 2*pi/360 / 60 
def least_ppi(L):
    A_eye=2*np.pi/360 / 60
    return 1/(2*L*39.37*np.tan(A_eye/2))

def ppi(px,size,rx=16,ry=9):
    py=px/rx*ry
    pd=np.sqrt(px**2+py**2)
    return pd/size
