#we using me and 
import numpy as np
#physical constant in SI unit
mu0=4*np.pi*1e-7
epsi0=8.8542e-12
c=299792458.
me=9.10000000000000006e-31
def gen_si(REAL_DX):
	#base units 4
	unit={};
	#m
	unit['l']=REAL_DX        
	#H/m
	unit['PERMEABILITY']=mu0
	#m/s
	unit['VELOCITY']=c
	#kg
	unit['MASS']=me
	#derived units
	#s
	unit['t']=unit['l']/unit['VELOCITY']
	#c
	unit['CHARGE']=np.sqrt(unit['MASS']*unit['l']/unit['PERMEABILITY'])
	#A/m**2
	unit['j']=unit['CHARGE']/unit['l']**2/unit['t']
	#T
	unit['b']=unit['MASS']/unit['t']/unit['CHARGE']
	#V/m
	unit['e']=unit['MASS']*unit['VELOCITY']/unit['t']/unit['CHARGE']
	#V
	unit['POTENTIAL']=unit['e']*unit['l']
	#W/m**3
	unit['POWER_DENSITY']=unit['MASS']/unit['l']/unit['t']**3
	#J
	unit['ENERGY']=unit['MASS']*unit['VELOCITY']**2
	#...
	unit['PERMITTIVITY']=(unit['t']**2*unit['CHARGE']**2)/(unit['MASS']*unit['l']**3)
	#constant
	const={}
	const['mu0']=1.
	const['epsi0']=1.
	const['c']=1.
	return unit,const
#not yet done
def gen_cgs(REAL_DX):
	#base units 4
	unit={};
	#m
	unit['l']=REAL_DX        
	#H/m
	unit['PERMEABILITY']=mu0
	#m/s
	unit['VELOCITY']=c
	#kg
	unit['MASS']=me
	#derived units
	#s
	unit['t']=unit['l']/unit['VELOCITY']
	#c
	unit['CHARGE']=np.sqrt(unit['MASS']*unit['l']/unit['PERMEABILITY'])
	#A/m**2
	unit['j']=unit['CHARGE']/unit['l']**2/unit['t']
	#T
	unit['b']=unit['MASS']/unit['t']/unit['CHARGE']
	#V/m
	unit['e']=unit['MASS']*unit['VELOCITY']/unit['t']/unit['CHARGE']
	#V
	unit['POTENTIAL']=unit['e']*unit['l']
	#W/m**3
	unit['POWER_DENSITY']=unit['MASS']/unit['l']/unit['t']**3
	#J
	unit['ENERGY']=unit['MASS']*unit['VELOCITY']**2
	#...
	unit['PERMITTIVITY']=(unit['t']**2*unit['CHARGE']**2)/(unit['MASS']*unit['l']**3)
	#constant
	const={}
	const['mu0']=1.
	const['epsi0']=1.
	const['c']=1.
	return unit,const
