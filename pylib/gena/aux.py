#dst time plot
import numpy as np
import pandas as pd
import datetime
def convert_to_timestamp(y, d, h):
    dt = datetime.datetime(int(y), 1, 1) + datetime.timedelta(days= int(d) - 1, hours=int(h))
    return dt.timestamp()

#indicator_data = np.loadtxt('./omni2_OCAaNiQuGH.lst')
def show_dst(ax, dst_file):
    indicator_data = np.loadtxt(dst_file)
    year = indicator_data[:,0].astype(int)
    day_of_year = indicator_data[:,1].astype(int)
    hour = indicator_data[:,2].astype(int)
    dst = indicator_data[:,3]
    #nt
    #W/m^2
    #lym = indicator_data[:,4]
    vec_convert = np.vectorize(convert_to_timestamp)
    timestamps = vec_convert(year, day_of_year, hour)
    ax0_date = pd.to_datetime(timestamps, unit='s', utc=True)
    ax.plot(ax0_date, dst, label='DST')
    ax.set_title('DST')
    ax.set_ylabel('nT')
    ax.legend()
    
