import numpy as np
c=1
ope=5
oce=1
omega = np.linspace(0.1,0.5,100)
k = np.sqrt(omega**2 + omega*ope**2/(oce-omega))
#assume vr_perp = 0.45
v=0.65
gm = 1/np.sqrt(1-(v/c)**2)
vr = (omega - oce/gm)/k

vg = 2*c**2*k/(2*omega + oce*ope**2/(oce-omega)**2)
delta2 = 1-(omega/k)**2
de = (1-vr/vg)**2*gm/delta2
R=0.25
coef = 2*np.pi*R/de