import numpy as np
def readchor(file,length):
    chor=np.loadtxt(file)
    steps = int(chor.shape[0]/length)
    dim = chor.shape[1]
    newc = chor.reshape([steps,length,dim])
    return newc
#old 
def loadvar(chorus,**kargs):
    for key, value in kargs.items():
        chorus[key]=value

def loadvec(chorus,name,fmt='.out',folder=''):
    if fmt == '.out':
        for var in name:
            chorus[var] = np.loadtxt(folder+vardict[var]+fmt)
    elif fmt == '.npy':
        for var in name:
            chorus[var]=np.load(folder+vardict[var]+fmt)
    else:
        print('not loaded')
        return 1
def loadwave(chorus,name,length,fmt='.out',folder=''):
    if fmt == '.out':
        for var in name:
            dat = readchor(folder+vardict[var]+fmt,length)
            #chorus[var]=dat
            chorus[var]=(dat[:,:-1,:]+dat[:,1:,:])/2
    elif fmt == '.npy':
        for var in name:
            dat = np.load(folder+vardict[var]+fmt)
            #chorus[var]=dat
            chorus[var]=(dat[:,:-1,:]+dat[:,1:,:])/2
    else:
        print('not loaded')
        return 1

vardict = {
'z':'zpos',
'k':'kmode',
'g':'gyro',
'd':'wp2',
'j':'Jact',
'gm':'gamma',
'vr':'vr',
'vg':'vg',
's':'source',
'c':'chorus'
}
def loadall(fmt='.out',f=''):
    conf={}
    wave={}
    loadvec(conf,['z','k','d','g','j','gm','vr','vg'],fmt=fmt,folder=f)
    N = len(conf['k'])+1
    loadvar(conf,**{'ns':N})
    loadwave(wave,['c'],N,fmt=fmt,folder=f)
    loadwave(wave,['s'],N,fmt=fmt,folder=f)
    loadvar(conf,**{'nt':wave['c'].shape[0]})
    wl = (conf['k']*conf['vr']+conf['g'])[int(conf['ns']/2)]
    loadvar(conf,**{'w':wl})
    conf['dt']=np.average(abs(np.gradient(conf['z'])/conf['vr']))
    return conf,wave
#cio.loadvec(res,['zpos','kmode','gyro','Jact','vr'],fmt='.npy')

#phase space
def loaddist(file,xi_o_s,ts,n):
    cnum=np.prod(xi_o_s)
    offset=ts*cnum*8
    dist = np.reshape(np.fromfile(file,offset=offset,count=cnum*n),np.append(xi_o_s,n),'F')
    return dist
def todist(file,d):
    with open("output.bin", "ab") as f:
        np.ravel(d,'F').tofile(f)
    return 0


def pull_phase(compA,dim=0):
    return np.unwrap(np.angle(compA),axis=dim)

# ------- hdf5
import h5py
#write hdf5 
from chorus.convert import *
def dumph5(name,grp,varname,var):
    with h5py.File(name, 'a') as f:
        if grp in f:
            group = f[grp]
        else:
            #group like dict
            group = f.create_group(grp)
        #dateset is array
        for vn,v in zip(varname,var):
            if vn in f[grp]:
                del f[grp][vn]
                group.create_dataset(vn, data=v)
            else:
                group.create_dataset(vn, data=v)

#only load one eveal return
def loadh5(name,grp,varname):
    with h5py.File(name, 'r') as f:
        dat = f[grp][varname]
        dim = len(dat.shape)
        if dim == 0:
            return dat[()]
        else:
            return dat[...]

def load_desc(ori='ori.h5',post='post.h5'):
    dori = " z:zpos, k:kmode, g:gyro, d:wp2, j:Jact, p:PI, gm:gamma, vr:vr, vg:vg, s:source, c:chorus, mu."
    dpost= 'contains: absc:|a|; abss:|s|; phasew:|ph ase wave|; phases:source phase; omega:frequency chirping; wavenumber:wavenumber chirping; vg:modified vg; E:Electric field slowly varying; B:Magnetic field slowly varying; Jh:real slowly varying hot current; Jc:slowly varying cold current from dispersion relation.'
        

def init_h5(ori='ori.h5',post='post.h5',fmt='.out',folder='',step=1):
    #check members
    #load des
    load_desc(ori,post)
    oridict = {
    'z':'zpos',
    'k':'kmode',
    'g':'gyro',
    'd':'wp2',
    'j':'Jact',
    'gm':'gamma',
    'vr':'vr',
    'vg':'vg',
    's':'source',
    'c':'chorus'
    }
#load cnf
    profiles = ['z','k','d','g','j','gm','vr','vg']
    derived = ['ns','w','dz','dt','p','m']
    z = np.loadtxt(folder+oridict['z']+fmt)[::step]
    k = np.loadtxt(folder+oridict['k']+fmt)[::step]
    d = np.loadtxt(folder+oridict['d']+fmt)[::step]
    g = np.loadtxt(folder+oridict['g']+fmt)[::step]
    j = np.loadtxt(folder+oridict['j']+fmt)[::step]
    gm= np.loadtxt(folder+oridict['gm']+fmt)[::step]
    vr= np.loadtxt(folder+oridict['vr']+fmt)[::step]
    vg= np.loadtxt(folder+oridict['vg']+fmt)[::step]

    #calculate w, ns, nt
    ns = len(k)+1
    wl = (k*vr+g)[int(ns/2)]
    dz = np.gradient(z)
    dt = np.average(abs(dz/vr)[1:-1])
    PI = (wl - g)/k**2
    mu = j + PI
    datas = [z,k,d,g,j,gm,vr,vg] + [ns,wl,dz,dt,PI,mu]
    dumph5(ori,'cnf',profiles+derived,datas)
    for varname in ['c','s']:
        dat = comp(readchor(folder+oridict[varname]+fmt,(ns-2)*step+2))
        dat=(dat[::step,:-1:step]+dat[::step,1::step])/2
        dumph5(ori,'wave',[varname],[dat])
    nt = dat.shape[0]
    dumph5(ori,'cnf',['nt'],[nt])


#post script
def g_wave(input='ori.h5',output='post.h5'):
    from chorus import wave
    #load required
    with h5py.File(input, 'r') as f:
        compchor =  f['wave']['c'][...]
        compsour =  f['wave']['s'][...]
        dt = f['cnf']['dt'][()]
        dz = f['cnf']['dz'][:]
        wl = f['cnf']['w'][()]
        kl = f['cnf']['k'][:]
        wpe = np.sqrt(f['cnf']['d'][:])
        wce = f['cnf']['g'][:]
    #complex
    #amplitude a and s
    absc= abs(compchor)
    abss= abs(compsour)
    #phase
    phasew = g_phase(compchor)
    phases = g_phase(compsour)
    #w and k
    dw,dk = wave.calcWK(phasew,dt,dz)
    omega=wl+dw
    wavenumber=kl+dk
    #B = {}
    vg = wave.calcVg(dw+wl,dk+kl,dt,dz)
    E = wave.E_evlp(compchor,dt,omega)
    B = wave.B_evlp(compchor,dz,wavenumber)
    Jh = wave.Jh_evlp(compsour)
    Jc = wave.Jc_evlp(E,omega,wpe,wce)
    varnames = ['absc','abss','phasew','phases','omega','wavenumber','vg','E','B','Jh','Jc']
    varvalues = [absc,abss,phasew,phases,omega,wavenumber,vg,E,B,Jh,Jc]
    dumph5(output,'wave',varnames,varvalues)

def g_ptc(ori='ori.h5',post='post.h5'):
    from chorus import ptc
    #load required
    with h5py.File(ori, 'r') as f:
        dz = f['cnf']['dz'][:]
        gyro = f['cnf']['g'][:]
        zpos = f['cnf']['z'][:]
        dt = f['cnf']['dt'][()]
        wl = f['cnf']['w'][()]
        kl = f['cnf']['k'][:]
        vr = f['cnf']['vr'][:]
        jact = f['cnf']['j'][:]
        compchor = f['wave']['c'][...]
    with h5py.File(post, 'r') as f:
        omg = f['wave']['omega'][()]
        knum = f['wave']['wavenumber'][:]
    #particle
    dgyro = np.gradient(gyro)/np.gradient(zpos)
    dkmode = np.gradient(kl)/np.gradient(zpos)
    dwdt =np.gradient(omg,axis=0)/dt
    PI = (wl - gyro)/kl**2
    #wb2 amp and vect
    dO = ptc.f_do(omg-wl,knum-kl,vr,kl)
    wb2_cplx = ptc.f_wb2(gyro,jact,compchor,omg,knum)
    #alpha and drg
    alpha = ptc.f_alpha_w(dwdt,dgyro,abs(wb2_cplx),knum,kl,jact,PI)
    #alpha1 = ptc.f_alpha_w(0,dgyro,wb2,wavenumber,res['j'],res['j'],PI)
    drg=jact/kl * dgyro - (wl-gyro)**2/kl**4 * dkmode
    #ret   0      1    2      3      4      5        6  7  8     9      10   11 
    #return absc,abss,phasew,phases,omega,wavenumber, E, B ,vg,wb2_cplx,alpha,drg
    varnames = ['wb2_cplx','alpha','drg','dO']
    varvalues = [wb2_cplx,alpha,drg,dO]
    dumph5(post,'ptc',varnames,varvalues)
    return # type: ignore
