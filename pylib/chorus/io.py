import numpy as np
def readchor(file,length):
    chor=np.loadtxt(file)
    steps = int(chor.shape[0]/length)
    dim = chor.shape[1]
    newc = chor.reshape([steps,length,dim])
    return newc

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

def pull_phase(compA,dim=0):
    return np.unwrap(np.angle(compA),axis=dim)


