import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import timedelta
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker
#rom matplotlib.dates import DateFormatter

from myplot import aux

from . import physical as phy
from . import filters as flt
from . import df as gdf

from lunardate import LunarDate

def show_channel(df,cnl):
    series = df[cnl]
    zero_count = (series == -94*0.081).sum()
    bins = [1, 100, 200, 300, 400, 500, 600, 1000]
    labels = ['1 - 100', '100 - 200', '200 - 300', '300 - 400', '400 - 500', '500 - 600', '>600']
    filtered_series = series[series > 0]
    cut_data = pd.cut(filtered_series, bins=bins, labels=labels, right=False)
    interval_counts = cut_data.value_counts()
    all_counts = pd.concat([pd.Series({'0': zero_count}), interval_counts])
    all_percentages = all_counts / all_counts.sum() * 100
    f,ax=plt.subplots()
    ax.pie(all_percentages, labels=all_percentages.index, autopct='%1.1f%%')
    return 0

def show_channel_bar(df, channel, cnl, folder):
    series = df[cnl]
    zero_count = (series == -94*0.081).sum()
    if channel == 'x':
        bins = [1, 100, 200, 300, 400, 500, 600]
        labels = ['1 - 100', '100 - 200', '200 - 300', '300 - 400', '400 - 500', '500 - 600']
    elif channel == 'y':
        bins = [1, 50, 100, 150, 200, 250, 300]
        labels = ['1 - 50', '50 - 100', '100 - 150', '150 - 200', '200 - 250', '250 - 300']
    filtered_series = series[series > 0]
    cut_data = pd.cut(filtered_series, bins=bins, labels=labels, right=False)
    interval_counts = cut_data.value_counts()
    # 重新索引，确保按照分箱顺序排列
    interval_counts = interval_counts.reindex(labels, fill_value=0)
    all_counts = pd.concat([pd.Series({'0': zero_count}), interval_counts])
    total_count = all_counts.sum()
    all_percentages = all_counts / total_count * 100


    f, ax = plt.subplots()
    # 绘制条形图
    bars = ax.bar(all_counts.index, all_counts.values)
    # 在每个条形上方添加百分比和计数标签
    for bar, (index, value), locx in zip(bars, all_counts.items(), np.linspace(0, 1, len(all_counts))):
        percentage = all_percentages[index]
        ax.text(locx, 0.2, f'{value}', transform=ax.transAxes, ha='center')
        ax.text(locx, 0.1, f'{percentage:.5f}%', transform=ax.transAxes, ha='center')


    ax.set_xlabel('Value Intervals')
    ax.set_ylabel('Count')
    ax.set_yticks([])
    ax.set_title(cnl)
    fig_w = 3.375
    f.set_size_inches(fig_w * 4, fig_w / 4 * 3)
    f.savefig(folder + cnl)



#only Date and ns col is enough
#construct it first
#return time_edges,window_histed

def df_resample(df, window='1min'):
    resampled = df.set_index('Date', drop=False).resample(window)
    return resampled

#resample hists location in windows
#level2 df
def hists_df_loc_in_window(df,vx=np.arange(-75,75),vy=np.arange(-30,30),window='10min'):
    resampled = df_resample(df,window)
    all_hist_matrix = []
    time_edges = []
    for time, group in resampled:
        time_edges.append(time)
    ht=[]
    for _, group in resampled:
        H, _ , _ = np.histogram2d(group['X1 (mm)'],group['Y1 (mm)'],[vx,vy])  if not group.empty else (np.zeros((len(vx)-1, len(vy)-1)),0,0)
        ht.append(H)
    return np.array(time_edges),np.array(ht)



def hist_df_ns_in_window(df,value_bins= np.arange(0, 1400, 2),window='1min'):
    reampled = df_resample(df,window)
    time_edges = []
    for time, _ in resampled:
        time_edges.append(time)
    
    window_histed = []
    for _, group in resampled:
        # 为每个特征列计算直方图
        #return 69? by 70(bins)
        hist, _ = np.histogram(group['values'], bins=value_bins) if not group.empty else (np.zeros(len(value_bins)-1),0)
        window_histed.append(hist)
    
    time_edges.append(time_edges[-1] + pd.Timedelta(window))
    time_edges = np.array(time_edges)
    return time_edges,window_histed
 
def add_bar(mesh,ax):
    ax_pos = ax.get_position()
    ax_width = ax_pos.x1 - ax_pos.x0
    ax_height = ax_pos.y1 - ax_pos.y0
    fraction = 0.25
    aspect = ax_width / (ax_height * fraction)

    cbar = plt.colorbar(mesh, ax=ax, orientation='horizontal', location='top',
                    fraction=fraction, pad=0.01, aspect=aspect)
    cbar.ax.set_ylabel('counts', rotation=0, labelpad=40)
    return 0
   
    
def show_startsum(axs,df):
    all_hist_matrix = []
    
    df1 = gdf.ENA_DataFrame({'Date': df['Date'], 'values': df['X1_start (ns)'] + df['X2_start (ns)']})
    df2 = gdf.ENA_DataFrame({'Date': df['Date'], 'values': df['Y1_start (ns)'] + df['Y2_start (ns)']})
    vb = np.arange(0,1400,2) 
    for df in [df1,df2]:
        time_edges, hist = hist_df_ns_in_window(df, value_bins= vb)
        all_hist_matrix.append(np.array(hist))
    
    vmax = np.max([np.max(hist) for hist in all_hist_matrix])
    vmin = np.min([np.min(hist) for hist in all_hist_matrix])
    
    for ax, hist_matrix in zip(axs, all_hist_matrix):
        mesh = ax.pcolormesh(
            time_edges,
            vb,
            hist_matrix.T,
            cmap='viridis',
            shading='flat',
            norm=LogNorm(vmin=1e-1, vmax=vmax)
        )
    add_bar(mesh,axs[0])
    
def show_startsum_old(axs,df,window='1min'):
    df.set_index('Date', inplace=True)
    resampled = df.resample(window)
    
    all_hist_matrix = []
    time_edges = []
    for time, _ in resampled:
        time_edges.append(time)
    
    value_bins = np.arange(0, 1400, 1)
    
    for (value,value_conj) in [['X1_start (ns)','X2_start (ns)'], ['Y1_start (ns)','Y2_start (ns)']]:
        window_hists = []
        for _, group in resampled:
            # 为每个特征列计算直方图
            #return 69? by 70(bins)
            hist, _ = np.histogram(group[value] + group[value_conj], bins=value_bins) if not group.empty else (np.zeros(len(value_bins)-1),0)
            window_hists.append(hist)
        all_hist_matrix.append(np.array(window_hists))
    
    # 转换矩阵结构 [时间窗口数, 特征数] -> [特征数, 时间窗口数]
    
    # 添加最后一个时间边缘
    time_edges.append(time_edges[-1] + pd.Timedelta(window))
    time_edges = np.array(time_edges)

    # 计算全局极值（保持颜色映射统一）
    vmax = np.max([np.max(hist) for hist in all_hist_matrix])
    #vmin = np.min([np.min(hist) for hist in all_hist_matrix])
    #print(vmax,vmin)
    # 绘制每个子图（修改以下部分）
    # 新增坐标变换函数（关键修改部分）
    #def forward(x):
    #    """将物理坐标映射到显示坐标"""
    #    return np.where(x <= 80, x, (x - 80)/10 + 80)
    
    #def inverse(x):
    #    """将显示坐标映射回物理坐标"""
    #    return np.where(x <= 80, x, (x - 80)*10 + 80)

    
    for ax, hist_matrix in zip(axs, all_hist_matrix):
        # 绘制时使用物理坐标分箱
        mesh = ax.pcolormesh(
            time_edges,
            value_bins,  # 使用原始分箱边界
            hist_matrix.T,
            cmap='viridis',
            shading='flat',
            norm=LogNorm(vmin=1e-1, vmax=vmax)
        )
        
#        # 设置双刻度（关键视觉效果）
#        #ax.yaxis.set_major_locator(ticker.FixedLocator([0, 20, 40, 60, 80, 200, 400, 600]))
#        #ax.yaxis.set_major_formatter(ticker.FuncFormatter(
#        #    lambda x, pos: f'{inverse(x):.0f}' if x > 80 else f'{x:.0f}'
#        #))
#    
#        # 设置自定义缩放
#        #scale = FuncScale(ax, functions=(forward, inverse))
#        #ax.set_yscale(scale)
    add_bar(mesh,axs[0])
#    #cbar.ax.set_xticks([vmin,vmax])  # 替换原来的LogFormatter
#    #cbar.formatter = ticker.LogFormatterSciNotation()  # 替换原来的LogFormatter

def show_hists_8(axs,df):
    all_hist_matrix = []
    vb = np.arange(0, 700, 1)
    for name in ['X1_start (ns)','X2_start (ns)', 'Y1_start (ns)','Y2_start (ns)', 'X1_end (ns)','X2_end (ns)', 'Y1_end (ns)','Y2_end (ns)']:
        subdf = gdf.ENA_DataFrame({'Date': df['Date'], 'values': df[name]})
        time_edges, hist = hist_df_ns_in_window(subdf, value_bins= vb)
        all_hist_matrix.append(np.array(hist))
    
    vmax = np.max([np.max(hist) for hist in all_hist_matrix])
    #vmin = np.min([np.min(hist) for hist in all_hist_matrix])
    #print(vmax,vmin)
    # 绘制每个子图（修改以下部分）
    # 新增坐标变换函数（关键修改部分）
    #def forward(x):
    #    """将物理坐标映射到显示坐标"""
    #    return np.where(x <= 80, x, (x - 80)/10 + 80)
    
    #def inverse(x):
    #    """将显示坐标映射回物理坐标"""
    #    return np.where(x <= 80, x, (x - 80)*10 + 80)

    
    for ax, hist_matrix in zip(axs, all_hist_matrix):
        # 绘制时使用物理坐标分箱
        mesh = ax.pcolormesh(
            time_edges,
            vb,  # 使用原始分箱边界
            hist_matrix.T,
            cmap='viridis',
            shading='flat',
            norm=LogNorm(vmin=1e-1, vmax=vmax)
        )
        
        # 设置双刻度（关键视觉效果）
        #ax.yaxis.set_major_locator(ticker.FixedLocator([0, 20, 40, 60, 80, 200, 400, 600]))
        #ax.yaxis.set_major_formatter(ticker.FuncFormatter(
        #    lambda x, pos: f'{inverse(x):.0f}' if x > 80 else f'{x:.0f}'
        #))
    
        # 设置自定义缩放
        #scale = FuncScale(ax, functions=(forward, inverse))
        #ax.set_yscale(scale)
    # 设置 colorbar 位置和 aspect
    add_bar(mesh,axs[0])
    #cbar.ax.set_xticks([vmin,vmax])  # 替换原来的LogFormatter
    #cbar.formatter = ticker.LogFormatterSciNotation()  # 替换原来的LogFormatter
    axs[0].set_ylabel('X1S (ns)')
    axs[1].set_ylabel('X2S (ns)')
    axs[2].set_ylabel('Y1S (ns)')
    axs[3].set_ylabel('Y2S (ns)')
    
    axs[4].set_ylabel('X1E (ns)')
    axs[5].set_ylabel('X2E (ns)')
    axs[6].set_ylabel('Y1E (ns)')
    axs[7].set_ylabel('Y2E (ns)')

    #X start
    for i in (0,1):
        axs[i].axhline(50,color='r',ls='--')
        axs[i].axhline(25,color='k',ls='-')
    #Y start
    for i in (2,3):
        axs[i].axhline(40,color='r',ls='--')
        axs[i].axhline(20,color='k',ls='-')
    
    #stop
    for i in (4,5,6,7):
        axs[i].axhline(50,color='r',ls='--')
        axs[i].axhline(100,color='r',ls='--')
    
    for i in (4,5,6,7,0,1,2,3):
        axs[i].set_yscale('log')
        axs[i].set_ylim(0.9,700)

    return 0

def show_paras(ax,r,flag,gap=10,mk='',msize=0):
    #rdate = pd.to_datetime(r.iloc[:, 0],format="%Y-%m-%d %H:%M:%S.%f")
    #ddate = pd.to_datetime(d.iloc[:, 0],format="%Y-%m-%d %H:%M:%S.%f")
    rdate = pd.to_datetime(r.iloc[:, 0])
    if flag == 'a':
        ax[0].plot(rdate[::gap], abs(r['HV_MCP1'][::gap]),marker=mk,markersize=msize)
        ax[1].plot(rdate[::gap], r['MCP1'][::gap],marker=mk,markersize=msize)
        ax[2].plot(rdate[::gap], r['AX'][::gap],marker=mk,markersize=msize,label = 'AY')
        ax[2].plot(rdate[::gap], r['AY'][::gap],marker=mk,markersize=msize,label = 'AX')
        ax[3].plot(rdate[::gap], r['TA'][::gap],marker=mk)
        ax[2].legend(frameon=False, ncol=3)
    elif flag == 'b':
        ax[0].plot(rdate[::gap], abs(r['HV_MCP2'][::gap]),marker=mk,markersize=msize)
        ax[1].plot(rdate[::gap], r['MCP2'][::gap],marker=mk,markersize=msize)
        ax[2].plot(rdate[::gap], (r['BX1'][::gap]*2+r['BG'][::gap]),marker=mk,markersize=msize,label = 'BX1')
        ax[2].plot(rdate[::gap], (r['BY1'][::gap]*2+r['BG'][::gap]),marker=mk,markersize=msize,label = 'BY1')
        ax[2].plot(rdate[::gap], (r['BX2'][::gap]*2+r['BG'][::gap]),marker=mk,markersize=msize,label = 'BX2')
        ax[2].plot(rdate[::gap], (r['BY2'][::gap]*2+r['BG'][::gap]),marker=mk,markersize=msize,label = 'BY2')
        ax[3].plot(rdate[::gap], r['TB'][::gap],marker=mk,markersize=msize)
        ax[2].legend(frameon=False, ncol=3)
    else:
        print('name not correct')
        exit(1)
    ax[0].set_ylabel('Volt')
    ax[0].set_ylim(2600,3000)
    
    ax[1].set_ylabel('MCP' )
    ax[1].set_ylim(150,1000)
    
    ax[2].set_ylabel('Anode')
    ax[2].set_ylim(20,700)
    
    #temp
    ax[3].set_ylabel('Temp')
    ax[3].set_ylim(5,50)
    minor_locator = ticker.MultipleLocator(1)
    ax[3].yaxis.set_minor_locator(minor_locator)
    return 0

#We can use timestamp and hist it to get counts with in chosen bins
#H should be in sorted order
def get_flux(H,t1,t2):
    if len(H) != 0:
        #ht = np.histogram(H,np.arange(H[0],np.array(H)[-1],60))
        ht = np.histogram(H-H[0],np.arange(t1,t2,60))
        return ht[1][1:],ht[0]
    else:
        return np.arange(t1,t2,60)[1:],np.arange(t1,t2,60)[1:]*0

def counts(ax,flux_t,ddate,tsecs,**keyargs):
    t,c = get_flux(flux_t,0,tsecs)
    timedelta_array = np.array([timedelta(seconds=int(val)) for val in t])
    ax.plot(ddate.iloc[0]+timedelta_array, c, '-', **keyargs)
    return ddate.iloc[0]+timedelta_array,c

def show_counts(ax,d,d4,d8,ns,match,onboard,msize=2,vmax=16*62*70*10):
    ddate = pd.to_datetime(d.iloc[:, 0])
    tsecs = ddate.iloc[-1].timestamp()-ddate.iloc[0].timestamp()
    
    t,call = counts(ax[0],np.array(d['Timestamp']),ddate,tsecs,markersize=msize, label = 'Trigger')
    
    counts(ax[0],np.array(d4['Timestamp']),ddate,tsecs,markersize=msize, label = 'Start')
    
    _,c8 = counts(ax[0],np.array(d8['Timestamp']),ddate,tsecs,markersize=msize, label = 'StartEnd')
    
    counts(ax[0],np.array(ns['Timestamp']),ddate,tsecs,markersize=msize, label = 'StartLim')
    
    _,cmatch = counts(ax[0],np.array(match['Timestamp']),ddate,tsecs,markersize=msize, label = 'SumMatch')
    
    counts(ax[0],np.array(onboard['Timestamp']),ddate,tsecs,markersize=msize, label = 'Onboard')
    
    ax[0].legend(frameon=False, ncol=3)
    ax[0].set_ylim(1,vmax)
    ax[0].axhline(16*62*60,color='r')
    
    #question mark
    ax[1].plot(t,c8/(call+0.0001),label = 'StartEnd/Trigger')
    ax[1].plot(t,cmatch/(c8+0.0001), label = 'SumMatch/StartEnd')
    ax[1].legend(frameon=False)
    ax[0].set_yscale('log')
    ax[1].set_yscale('log')
    ax[0].set_ylabel('Counts/Min')
    ax[1].set_ylabel('Ratio')
    ax[1].set_ylim(0.001,1)
    return 0



def show_paras_counts_hists(fig,ax,p_df,l1_df,df4,df8,ns,match,onboard,flag):
    #just plot all dst
    t1 = pd.Timestamp(pd.to_datetime(p_df['Timestamp'].iloc[0], unit='s', utc=True))
    t2 = pd.Timestamp(pd.to_datetime(p_df['Timestamp'].iloc[-1], unit='s', utc=True))
    
    show_paras(ax[0:4],p_df,flag)
    
    show_hists_8(ax[4:12],l1_df)
    #show_hists(fig,ax[5],l1_df,'X1_start (ns)',np.concatenate(([0,1],np.arange(10,700,10))))
    
    show_counts(ax[12:14],l1_df,df4,df8,ns,match,onboard)
    #ax[0].set_xlim(t1,t2)
        

    
    xdates =  pd.to_datetime(l1_df.iloc[:, 0])
    lunar_dates = []
    for date in xdates:
        lunar = LunarDate.fromSolarDate(date.year, date.month, date.day)
        lunar_str = f"{lunar.year}-{lunar.month}-{lunar.day}"
        lunar_dates.append(lunar_str)
    def format_func(x, pos=None):
        # 将Matplotlib数值转为公历日期
        gregorian = mdates.num2date(x).strftime('%Y-%m-%d  %H:%M')
        # 查找对应的农历日期
        idx = int(round(x - mdates.date2num(xdates.iloc[0])))  # 获取索引
        if 0 <= idx < len(lunar_dates):
            return f"{gregorian}\n{lunar_dates[idx]}"
        else:
            return gregorian
    ax[10].xaxis.set_major_formatter(format_func)
    #ax.xaxis.set_major_locator(mdates.MonthLocator())  # 每月一个刻度

    #date_format = mdates.DateFormatter('%Y-%m-%d %H:%M')
    #ax.xaxis.set_major_locator(mdates.MonthLocator())  # 每月一个刻度
    
    for _ in ax:
        _.grid(which='minor', linestyle='--', alpha=0.5)
        _.grid(which='major', linestyle='-', alpha=0.8)
    minor_locator = ticker.MultipleLocator(25)
    #HV, MCP, THE
    for ax_obj in [ax[1], ax[2], ax[3]]:
        minor_locator = ticker.MultipleLocator(25)
        ax_obj.yaxis.set_minor_locator(minor_locator)
    # 显示次要刻度的网格线
    fig.autofmt_xdate()
    
    #t1_str = t1.strftime("%Y-%m-%d %H:%M:%S")
    #t2_str = t2.strftime("%Y-%m-%d %H:%M:%S")
    #fig.suptitle(flag + f'{t1_str} to {t2_str}')


def show_file(prefix,fig,ax,flag,snid=1):
    phy.init_param(snid)
    p_file =  prefix + "_parameter.csv"
    l1_file = prefix + "_L1A_"+ flag + ".csv"
    p_df = gdf.read_csv(p_file)
    l1_df = gdf.read_csv(l1_file)
    _, df4 = flt.start_eff(l1_df)
    _, df8 = flt.start_stop_eff(df4)
    if flag == 'a':
        xlim = phy.A_NS_XLIM
        ylim = phy.A_NS_YLIM
    elif flag == 'b':
        xlim = phy.B_NS_XLIM
        ylim = phy.B_NS_YLIM
    _, ns = flt.start_ns_eff(df8,xlim,ylim)
    _, match = flt.match_eff(ns)
    level2_df = flt.l1_2_l2(match,flag,snid)
    _, onboard = flt.onboard_eff(level2_df)
    show_paras_counts_hists(fig,ax, p_df, l1_df, df4, df8, ns, match, onboard,  flag)
    fig.savefig(prefix +  '.jpg',dpi=300,bbox_inches='tight')



#quick
#plot from level 1
from scipy.stats import gaussian_kde
def show_start_end_level2(df,msize=1,s_or_hist=1):
    f,(ax1,ax2)=plt.subplots(2,1)
    if s_or_hist == 0:
        ax1.scatter(df['X1 (mm)'], df['Y1 (mm)'], label = 'start',s=msize)
        ax2.scatter(df['X2 (mm)'], df['Y2 (mm)'], label = 'end',s=msize)
        ax1.legend()
        ax2.legend()
    elif s_or_hist == 1:
        hh = ax1.hist2d(df['X1 (mm)'], df['Y1 (mm)'],np.arange(-100,100,10))
        f.colorbar(hh[3], ax=ax1)
        hh = ax2.hist2d(df['X2 (mm)'], df['Y2 (mm)'],np.arange(-100,100,10))
        f.colorbar(hh[3], ax=ax2)
    elif s_or_hist == 2:
        xy = np.vstack([df['X1 (mm)'], df['Y1 (mm)']])
        z = 1-gaussian_kde(xy)(xy)
        ax1.scatter(df['X1 (mm)'], df['Y1 (mm)'], c=z, s=1, cmap='Spectral')
        xy = np.vstack([df['X2 (mm)'], df['Y2 (mm)']])
        z = 1-gaussian_kde(xy)(xy)
        ax2.scatter(df['X2 (mm)'], df['Y2 (mm)'], c=z, s=1, cmap='Spectral')

    aux.draw_rectangle(ax1,90,30)
    aux.draw_rectangle(ax2,90,30)
    aux.draw_rectangle(ax1,150,60)
    aux.draw_rectangle(ax2,150,60)
    ax1.set_xlabel('X (mm)')
    ax1.set_ylabel('Y (mm)')
    ax2.set_xlabel('X (mm)')
    ax2.set_ylabel('Y (mm)')
    return f,(ax1,ax2)

def show_start_end_level1(df,flag,snid=1,msize=1,s_or_hist=1):
    df2 = flt.l1_2_l2(df,flag,snid)
    return show_start_end_level2(df2,msize,s_or_hist)
    
def mid_ary(ary):
    return (ary[:-1]+ary[1:])/2

def show_angle_level2(df,dx=2.5):
    x1,y1 = (df['X1 (mm)'], df['Y1 (mm)'])
    x2,y2 = (df['X2 (mm)'], df['Y2 (mm)'])
    px = np.arctan((x1-x2)/50)*180/np.pi
    py = np.arctan((y1-y2)/50)*180/np.pi
    f,ax=plt.subplots()
    grid_x = np.arange(-60,60,dx)
    grid_y = np.arange(-30,30,dx)
    H = np.histogram2d(px,py,[grid_x,grid_y])
    Z,X = np.meshgrid(mid_ary(grid_y),mid_ary(grid_x))
    ax.pcolormesh(X,Z,H[0])
    aux.draw_rectangle(ax,120,60)
    return f,ax

#cartoon

from matplotlib import animation
def cartoon2d_start(xyextent,datetimes,array,*args,**kwargs):
	fig=plt.figure()
	ax=plt.axes()
	time = ax.annotate(0,xy=(0.2, 0.9),xycoords='figure fraction')
	im = ax.imshow(array[0,...],aspect='auto',extent=xyextent,*args,**kwargs)
	cbar = fig.colorbar(im)
	def animate(i):
		Z=array[:,:,i]
		im.set_array(Z)
		im_min, im_max = np.min(Z), np.max(Z)
		im.set_clim(im_min, im_max)
		time.set_text(datetimes[i])
	anim = animation.FuncAnimation(fig, animate, *args,**kwargs)
	return anim

if __name__ == '__main__':
    pass
