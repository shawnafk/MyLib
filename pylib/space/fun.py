import numpy as np
import matplotlib.pyplot  as plt
from datetime import datetime, timedelta
from astropy.time import Time
import pytz
import pyspedas
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.dates as mdates

RE = 6371.393
# ---------------------------- rotation ------------------------------
def q2rot(q1,q2,q3,N):
    mat = np.zeros((N,3,3))
    q0 = np.sqrt(1-q1**2-q2**2-q3**2)
    mat[:,0,0] = 1 - 2*q2**2 - 2* q3**2
    mat[:,0,1] = 2*q1*q2 - 2* q0 * q3
    mat[:,0,2] = 2 * q1 * q3 + 2 * q0 * q2
         
    mat[:,1,0] = 2*q1*q2 + 2* q0 * q3
    mat[:,1,1] = 1 - 2*q1**2 - 2* q3**2
    mat[:,1,2] = 2 * q2 * q3 - 2 * q0 * q1
         
    mat[:,2,0] = 2*q1*q3 - 2* q0 * q2
    mat[:,2,1] = 2*q2*q3 + 2* q0 * q1
    mat[:,2,2] = 1 - 2*q1**2 - 2*q2**2
    return mat

#only one
def rot_x(deg):
    r11 = 1
    r12 = 0
    r13 = 0
    
    r21 = 0
    r22 = np.cos(deg)
    r23 = -np.sin(deg)
    
    r31 = 0
    r32 = np.sin(deg)
    r33 = np.cos(deg)
    
    mat = np.array([[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]])
    return mat

def rot_y(deg):
    r11 = np.cos(deg)
    r12 = 0
    r13 = np.sin(deg)
    
    r21 = 0
    r22 = 1
    r23 = 0
    
    r31 = -np.sin(deg)
    r32 = 0
    r33 = np.cos(deg)
    
    mat = np.array([[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]])
    return mat

def rot_z(deg):
    r11 = np.cos(deg)
    r12 = -np.sin(deg)
    r13 = 0
    
    r21 = np.sin(deg)
    r22 = np.cos(deg)
    r23 = 0
    
    r31 = 0
    r32 = 0
    r33 = 1
    
    mat = np.array([[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]])
    return mat

# ---------------------------- time manipulation ------------------------------

def strtime2timestamp(time_str,format = '%Y-%m-%d %H:%M:%S'):
    N  = time_str.shape[0]
    dt = np.zeros(N)
    for i in range(N):
        dt[i]=datetime.strptime(time_str[i][:19], format).timestamp()
    return dt

def timestamp2isot(ts):
    N  = ts.shape[0]
    isot=np.zeros(N,dtype='object')
    for i in range(N):
        dt_object = datetime.fromtimestamp(ts[i], pytz.timezone('UTC'))
        #iso_format_string = dt_object.isoformat()
        isot[i] = Time(dt_object, scale='utc')
    return isot

#to juliandate
def isot2jd(isot):
    N  = isot.shape[0]
    t_jd=np.zeros(N)
    for i in range(N):
        t_jd[i] = isot[i].jd   
    return t_jd

def str_to_timestamp(time_str):
    # 将字符串按照'年月日时分秒'的格式解析为datetime对象
    dt = datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S')
    # 将datetime对象转换为Unix时间戳（单位为秒）
    timestamp = int(dt.timestamp())
    return timestamp

def timestamp_to_str(timestamp):
    # 将Unix时间戳转换为datetime对象
    dt = datetime.fromtimestamp(timestamp)
    # 将datetime对象格式化为'年月日时分秒'的字符串
    time_str = dt.strftime('%Y-%m-%d %H:%M:%S')
    return time_str

def find_closest_index(timestamp_list, t0):
    left, right = 0, len(timestamp_list) - 1
    closest_index = -1
    min_diff = float('inf')  # 正无穷大

    while left <= right:
        mid = (left + right) // 2
        current_diff = abs(timestamp_list[mid] - t0)

        # 如果当前的差值小于之前的最小差值，则更新最接近的索引和最小差值
        if current_diff < min_diff:
            min_diff = current_diff
            closest_index = mid

        # 如果当前时间戳小于t0，则在右半部分继续查找
        if timestamp_list[mid] < t0:
            left = mid + 1
        # 如果当前时间戳大于t0，则在左半部分继续查找
        else:
            right = mid - 1

    return closest_index

def get_idx_of_times(t0,dt,n,ts):
    t0 = str_to_timestamp(t0)
    target_time = t0 + np.arange(0,n)*dt
    idx = []
    for i in range(n):
        idx.append(find_closest_index(ts,target_time[i]))
    return idx

def convert_to_datetime(year, day_of_year, hour_fraction):
    start_date = datetime(year, 1, 1)
    
    # 将天数转换为timedelta对象并添加到开始日期
    delta = timedelta(days=day_of_year - 1)
    date = start_date + delta
    
    # 将小时分数转换为小时、分钟和秒
    all_sec = hour_fraction * 3600
    hours = all_sec // 3600
    minutes = (all_sec % 3600) // 60
    seconds = all_sec % 60
    res_microsec = all_sec % 1 * 1e6
    
    # 创建datetime对象
    datetime_obj = datetime(year, date.month, date.day, int(hours), int(minutes), int(seconds), int(res_microsec))
    # 返回格式化的字符串
    #return datetime_obj.strftime('%Y-%m-%d-%H-%M-%S.%f')
    return datetime_obj


def load_ibex(full_name):
    data = np.loadtxt(full_name,usecols=np.arange(2,76))
    #data has correspondings
    # time: 0 1 2
    # location(gse x y z): 3 4 5
    # location(R): 6
    # zaxis RA and DEC: 7 8
    # moon x y z: 9 10 11
    # moon R: 12
    # spin counts: 13
    # 0 -- 354 bins: 14 : 74
    epoch = []
    for row in data:
        epoch.append(convert_to_datetime(int(row[0]),int(row[1]),row[2]))
    orbit = data[:,3:6]
    ax_pointing = data[:,7:9]
    counts = data[:,14:]
    return np.array(epoch),orbit,ax_pointing,counts

#need bow shock, magnetospause boundary
def BS_args(B_z,beta,M_ms,D_p):
    a_1=11.1266
    a_2=0.0010
    a_3=-0.0005
    a_4=2.5966
    a_5=0.8182
    a_6=-0.0170
    a_7=-0.0122
    a_8=1.3007
    a_9=-0.0049
    a_10=-0.0328
    a_11=6.047
    a_12=1.029
    a_13=0.0231
    a_14=-0.002
    if B_z >0:
        r_0=a_1 * (1+a_2*B_z)*(1+a_9 * beta ) * (1+a_4 * ((a_8-1)*M_ms**2+2)/ (a_8+1)/M_ms**2) * D_p**(-1 / a_11)
        alpha = a_5 * (1+a_13 * B_z)*(1+ a_7 * D_p)*(1+a_10*np.log(1+beta)) * (1+a_14 * M_ms)
    else:
        r_0=a_1 * (1+a_3*B_z)*(1+a_9 * beta ) * (1+a_4 * ((a_8-1)*M_ms**2+2)/ (a_8+1)/M_ms**2) * D_p**(-1/a_11)
        alpha=a_5*(1+a_6 * B_z)*(1+ a_7 * D_p)*(1+a_10*np.log(1+beta)) * (1+a_14 * M_ms)
    epsi = a_12
    return r_0, alpha,epsi

def BS(B_z,beta,M_ms,D_p):
    r0,alp,epsi =  BS_args(B_z,beta,M_ms,D_p)
    theta = np.linspace(0,2*np.pi,128)
    r = r0 * ( (1+epsi)/(1+epsi * np.cos(theta) ) ) ** alp
    x = r * np.cos(theta)
    R = r * np.sin(theta)
    return x,R

def MP_args(B_z,D_p):
    a_1 = 11.646
    a_2 = 0.216
    a_3 = 0.122
    a_4 = 6.215
    a_5 = 0.578
    a_6 = -0.009
    a_7 = 0.012
    if B_z > 0:
        r0  = a_1 * D_p**(-1/a_4)
    elif B_z < 0 and B_z > -8:
        r0 = (a_1 + a_2* B_z) * D_p**(-1/a_4)
    else:
        r0 =( a_1 + 8*a_3 - 8 * a_2 + a_3 * B_z ) * D_p**(-1/a_4)
    alpha = (a_5 + a_6 * B_z) * (1+ a_7 * D_p)
    return r0, alpha

def MP(B_z,D_p):
    r0,alp = MP_args(B_z,D_p)
    theta = np.linspace(0,2*np.pi,128)
    r = r0 * ( 2/(1+np.cos(theta) ) ) ** alp
    x = r * np.cos(theta)
    R = r * np.sin(theta)
    return x,R

def rot_ra(ra):
    r11 = 1
    r12 = 0
    r13 = 0

    r21 = 0
    r22 = np.cos(ra)
    r23 = -np.sin(ra)

    r31 = 0
    r32 = np.sin(ra)
    r33 = np.cos(ra)

    mat = np.array([[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]])
    return mat

def rot_dec(dec):
    r11 = np.cos(dec)
    r12 = -np.sin(dec)
    r13 = 0

    r21 = np.sin(dec)
    r22 = np.cos(dec)
    r23 = 0

    r31 = 0
    r32 = 0
    r33 = 1

    mat = np.array([[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]])
    return mat


#accumulate counts
#if look from dawn,  Y > 0, Earth lies in 0 - 180 view angle  , flag = 0
#if look from dusk,  Y < 0, Earth lies in 180 - 360 view angle , flag = 1
#def xz_counts(orbit,pointing,counts,flag):
#    rows = counts.shape[0] 
#    nbins = 30
#    coords = np.zeros(((nbins)*rows,2))
#    counts = np.zeros(((nbins)*rows))
#    
#    ang_60 = (np.linspace(0,354,nbins*2) + 3)/180*np.pi
#    vec_60 = np.zeros((nbins*2,3))
#    vec_60[:,1] = np.sin(ang_60)
#    vec_60[:,2] = np.cos(ang_60)
#    #0 - 180
#    #if flag == 0:
#    #   ang = ang_360[:nbins]
#    #180 - 360
#    #elif flag == 1:
#    #   ang = ang_360[nbins:] 
#    #ang = np.linspace(180,360,nbins)
#    #map_coef = np.zeros(nbins)
#    #map_coef = np.tan((ang -3 - 90)/180*np.pi)
#    #map_coef[0] = np.inf
#    #map_coef[-1] = -np.inf
#    for r in range(rows):
#        #constructure rotation matrix
#        RA = pointing[r,0]/180*np.pi
#        DEC = pointing[r,1]/180*np.pi
#        Rx = rot_ra(RA)
#        Rz = rot_dec(DEC)
#        Rot = Rx@Rz
#        for i in range(nbins):
#            v = Rot@vec[i,:]
#        
#        #x localtions at select time intervers
#        #angle intervel
#        #z/y = tan(t) 
#        #z: t_th * 60 intervals
#        y = data[r,4]
#        #treat z as view center
#        z = y * map_coef
#        
#        Ra = data[r,7]
#        x_ra = data[r,3] + y*np.tan(Ra/180*np.pi)
#        #15 is the first 15 args
#        if flag == 0:
#            bin_s = 14
#        else:
#            bin_s = 14 + nbins
#        for i in range(nbins):
#            #if i == 0 or i == nbins - 1:
#            #    continue
#            #else:
#            #x
#            coords[i+(nbins)*r,0] = x_ra
#            #z
#            coords[i+(nbins)*r,1] = z[i]
#            #value
#            counts[i+(nbins)*r] = data[r,bin_s+i]
#    return coords,counts

def load_omni(t,str):
    date1 = t[0].strftime('%Y-%m-%d-%H-%M')
    date2 = t[-1].strftime('%Y-%m-%d-%H-%M')
    omni_vars = pyspedas.omni.data(trange=[date1, date2],notplot=True,varnames = [str])
    return omni_vars
#get X-R(Z/Y, assuming axisymmectric) location of MP and BS
def loc_bs_mp(t):
    date1 = t[0].strftime('%Y-%m-%d-%H-%M')
    date2 = t[-1].strftime('%Y-%m-%d-%H-%M')
    omni_vars = pyspedas.omni.data(trange=[date1, date2],notplot=True,varnames = ['BZ_GSM','Mgs_mach_num','Beta','Pressure'])
    def ave_array(array_with_nans):
        clean_array = array_with_nans[~np.isnan(array_with_nans)]
        return np.average(clean_array) 

    IMF_Bz = ave_array(omni_vars['BZ_GSM']['y'])
    D_p = ave_array(omni_vars['Pressure']['y'])
    Beta = ave_array(omni_vars['Beta']['y'])
    M_ms = ave_array(omni_vars['Mgs_mach_num']['y'])

    x_mp,R_mp = MP(IMF_Bz,D_p)
    x_bs,R_bs = BS(IMF_Bz,Beta,M_ms,D_p)
    return x_mp,R_mp, x_bs,R_bs

def plot_orbit(ax,t,x,y):
    tdays =  mdates.date2num(t)
    dt = tdays-tdays[0]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # 创建一个颜色映射
    cmap = plt.get_cmap('viridis')
    norm = plt.Normalize(dt.min(), dt.max())

    # 创建LineCollection对象
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    # 将时间值赋给LineCollection对象
    lc.set_array(dt)
    lc.set_linewidth(2)
    ax.add_collection(lc)
    ax.autoscale()

    # 添加颜色条
    #cbar = plt.colorbar(lc, ax=ax)
    #cbar.set_label('Time')
    return 0

def xz_counts(orbit,pointing,counts,flag):
    rows = counts.shape[0] 
    nbins = 30
    coords = np.zeros(((nbins)*rows,2))
    
    ang_60 = (np.linspace(0,354,nbins*2) + 3)/180*np.pi
    vec_60 = np.zeros((nbins*2,3))
    vec_60[:,1] = np.sin(ang_60)
    vec_60[:,2] = np.cos(ang_60)
    #0 - 180
    #if flag == 0:
    #   ang = ang_360[:nbins]
    #180 - 360
    #elif flag == 1:
    #   ang = ang_360[nbins:] 
    #ang = np.linspace(180,360,nbins)
    #map_coef = np.zeros(nbins)
    #map_coef = np.tan((ang -3 - 90)/180*np.pi)
    #map_coef[0] = np.inf
    #map_coef[-1] = -np.inf
    
    recounts = np.zeros(((nbins)*rows))
    #15 is the first 15 args
    if flag == 0:
        bin_start = 0
    else:
        bin_start = nbins
    for r in range(rows):
        #constructure rotation matrix
        RA = pointing[r,0]/180*np.pi
        DEC = pointing[r,1]/180*np.pi
        Rx = rot_ra(RA)
        Rz = rot_dec(DEC)
        Rot = Rx@Rz
        for i in range(nbins):
            v = Rot@vec_60[i+bin_start,:]
            #start
            #find cross section
            # y0 + dt * vec_y = 0
            t = (0 - orbit[r,1])/v[1]
            x_proj = (orbit[r,0] + t * v[0])/RE
            z_proj = (orbit[r,2] + t * v[2])/RE
            
            #x localtions at select time intervers
            #angle intervel
            #z/y = tan(t) 
            #z: t_th * 60 intervals
            #y = data[r,4]
            #treat z as view center
            #z = y * map_coef
            #
            #Ra = data[r,7]
            #x_ra = data[r,3] + y*np.tan(Ra/180*np.pi)
        
            #if i == 0 or i == nbins - 1:
            #    continue
            #else:
            #x
            coords[i+(nbins)*r,0] = x_proj
            #z
            coords[i+(nbins)*r,1] = z_proj
            #value
            recounts[i+(nbins)*r] = counts[r,bin_start+i]
    return coords,recounts

import os
from space import fun
import PyGeopack as gp
def plot_fieldline(time,ax,c):
    f = 'BF'+time+'.npy'
    if os.path.exists(f):
        Tall = np.load(f,allow_pickle=True)
        for i in range(30):
            for j in range(20):
                ax.plot(Tall[j].xgse[i,:],Tall[j].zgse[i,:],c,lw=0.15)
    else:
        Tall = []
        for i in np.arange(-10,10):
            #x, y, z
            Tall.append(gp.TraceField(np.linspace(-20,20,40),0,i,int(time),2))
        for i in range(30):
            for j in range(20):
                ax.plot(Tall[j].xgse[i,:],Tall[j].zgse[i,:],c,lw=0.15)
        np.save(f,Tall)
    fun.plot_earth(ax)
    return 0

def show_FOV(ax,vertex):
    x=[]
    y=[]
    for p in vertex:
        x.append(p[0])
        y.append(p[1])
        ax.scatter(p[0],p[1],marker='+',color='k')
    x.append(vertex[0][0])
    y.append(vertex[0][1])
    ax.plot(x,y)
    return 0

def plot_earth(ax):
    t = np.linspace(0,2*np.pi,128)
    x = np.cos(t)
    y = np.sin(t)
    ax.plot(x,y,'k')
    return 0



def nearest_value(t0,epoch,value):
    epoch2ts = [dt.timestamp() for dt in epoch]
    idx = find_closest_index(epoch2ts,t0)
    return value[idx]
    

