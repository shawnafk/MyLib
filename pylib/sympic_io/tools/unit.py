#we using me and 
import numpy as np
#physical constant in SI unit
mu0=4*np.pi*1e-7
epsi0=8.8542e-12
c=299792458.
me=9.10000000000000006e-31
def gen_unit(REAL_DX):
	#base units 4
	unit={};
	#m
	unit['LENGTH']=REAL_DX        
	#H/m
	unit['PERMEABILITY']=mu0
	#m/s
	unit['VELOCITY']=c
	#kg
	unit['MASS']=me
	#derived units
	#s
	unit['TIME']=unit['LENGTH']/unit['VELOCITY']
	#c
	unit['CHARGE']=np.sqrt(unit['MASS']*unit['LENGTH']/unit['PERMEABILITY'])
	#A/m**2
	unit['CURRENT_DENSITY']=unit['CHARGE']/unit['LENGTH']**2/unit['TIME']
	#T
	unit['B']=unit['MASS']/unit['TIME']/unit['CHARGE']
	#V/m
	unit['E']=unit['MASS']*unit['VELOCITY']/unit['TIME']/unit['CHARGE']
	#V
	unit['POTENTIAL']=unit['E']*unit['LENGTH']
	#W/m**3
	unit['POWER_DENSITY']=unit['MASS']/unit['LENGTH']/unit['TIME']**3
	#J
	unit['ENERGY']=unit['MASS']*unit['VELOCITY']**2
	#...
	unit['PERMITTIVITY']=(unit['TIME']**2*unit['CHARGE']**2)/(unit['MASS']*unit['LENGTH']**3)
	#constant
	const={}
	const['mu0']=1.
	const['epsi0']=1.
	const['c']=1.
	return unit,const
