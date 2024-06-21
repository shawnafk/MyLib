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
            chorus[var]=readchor(folder+vardict[var]+fmt,length)
    elif fmt == '.npy':
        for var in name:
            chorus[var]=np.load(folder+vardict[var]+fmt)
    else:
        print('not loaded')
        return 1
    loadvar(chorus,**{'t':chorus[var].shape[0]})

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
    res={}
    loadvec(res,['z','k','d','g','j','gm','vr','vg'],fmt=fmt,folder=f)
    wl = res['k']*res['vr']+res['g']
    N = len(res['k'])+1
    loadvar(res,**{'w':wl})
    loadvar(res,**{'n':N})
    loadwave(res,['c'],N,fmt=fmt,folder=f)
    loadwave(res,['s'],N,fmt=fmt,folder=f)
    res['dt']=np.average(abs(np.gradient(res['z'])/res['vr']))
    return res
#cio.loadvec(res,['zpos','kmode','gyro','Jact','vr'],fmt='.npy')

#phase space
def loaddist(file,xi_o_s,ts,n):
    cnum=np.prod(xi_o_s)
    offset=ts*cnum*8
    dist = np.reshape(np.fromfile(file,offset=offset,count=cnum*n),np.append(xi_o_s,n),'F')
    return dist

def pull_phase(compA,dim=0):
    return np.unwrap(np.angle(compA),axis=dim)


