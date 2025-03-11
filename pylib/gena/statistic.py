import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import timedelta
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker
from myplot import aux
import numpy as np
from . import csv0
import pandas as pd
from . import physical as phy
from . import filters as flt
#from level 1A
#return total
def ratio_levels(df,prob='A',snid=1):
    phy.init_param(snid)
    if prob == 'A':
        xlim = phy.A_NS_XLIM
        ylim = phy.A_NS_YLIM
    if prob == 'B':
        xlim = phy.B_NS_XLIM
        ylim = phy.B_NS_YLIM
    len1,df1  = flt.start_end_eff(df)
    len2,df2  = flt.start_end_eff(df1)
    len3,df3  = flt.start_ns_eff(df2,xlim,ylim)
    len4,df4  = flt.match_eff(df3)
    len5,_  = flt.onboard_eff(df4)
    return len1,len2,len3,len4,len5



def show_hists_8(fig,axs,df,window='5min'):
    # 对于 0 - 1 区间，保持线性关系
    # 对于大于 1 的值，进行线性变换使得 10 - 20 区间和 0 - 1 区间等宽
    # Step 1: 转换时间为datetime类型并设置索引
    df['timestamp'] = pd.to_datetime(df['Date'])
    df.set_index('timestamp', inplace=True)
    
    # 预先进行 resample 操作
    resampled = df.resample(window)
    
    # 初始化存储矩阵和边缘时间
    all_hist_matrix = []
    time_edges = []
    for time, _ in resampled:
        time_edges.append(time)
    
    # 生成组合分箱：0-80用1ns间隔，80-600用10ns间隔
    #value_bins = np.concatenate([
    #    np.arange(0, 80, 1),
    #    np.arange(80, 600, 10)
    #])
    value_bins = np.arange(0, 700, 1)
    # 遍历每个时间窗口（只循环一次）
    for value in ['X1_start (ns)','X2_start (ns)', 'Y1_start (ns)','Y2_start (ns)', 'X1_end (ns)','X2_end (ns)', 'Y1_end (ns)','Y2_end (ns)']:
        window_hists = []
        for _, group in resampled:
            # 为每个特征列计算直方图
            #return 69? by 70(bins)
            hist, _ = np.histogram(group[value], bins=value_bins) if not group.empty else (np.zeros(len(value_bins)-1),0)
            window_hists.append(hist)
        all_hist_matrix.append(np.array(window_hists))
    
    # 转换矩阵结构 [时间窗口数, 特征数] -> [特征数, 时间窗口数]
    
    # 添加最后一个时间边缘
    time_edges.append(time_edges[-1] + pd.Timedelta(window))
    time_edges = np.array(time_edges)
    value_edges = value_bins
    
    # 计算全局极值（保持颜色映射统一）
    vmax = np.max([np.max(hist) for hist in all_hist_matrix])
    vmin = np.min([np.min(hist) for hist in all_hist_matrix])
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
        
        # 设置双刻度（关键视觉效果）
        #ax.yaxis.set_major_locator(ticker.FixedLocator([0, 20, 40, 60, 80, 200, 400, 600]))
        #ax.yaxis.set_major_formatter(ticker.FuncFormatter(
        #    lambda x, pos: f'{inverse(x):.0f}' if x > 80 else f'{x:.0f}'
        #))
    
        # 设置自定义缩放
        #scale = FuncScale(ax, functions=(forward, inverse))
        #ax.set_yscale(scale)
    # 设置 colorbar 位置和 aspect
    ax_pos = axs[0].get_position()
    ax_width = ax_pos.x1 - ax_pos.x0
    ax_height = ax_pos.y1 - ax_pos.y0
    fraction = 0.25
    aspect = ax_width / (ax_height * fraction)

    # 创建 colorbar 并设置 aspect
    cbar = fig.colorbar(mesh, ax=axs[0], orientation='horizontal', location='top',
                    fraction=fraction, pad=0.01, aspect=aspect)
    cbar.ax.set_ylabel('counts', rotation=0, labelpad=40)
    cbar.formatter = ticker.LogFormatter()  # 明确指定对数格式


def show_paras(ax,r,flag,gap=10,mk='',msize=0):
    #rdate = pd.to_datetime(r.iloc[:, 0],format="%Y-%m-%d %H:%M:%S.%f")
    #ddate = pd.to_datetime(d.iloc[:, 0],format="%Y-%m-%d %H:%M:%S.%f")
    rdate = pd.to_datetime(r.iloc[:, 0])
    if flag == 'A':
        ax[0].plot(rdate[::gap], abs(r['HV_MCP1'][::gap]),marker=mk,markersize=msize)
        ax[1].plot(rdate[::gap], r['MCP1'][::gap],marker=mk,markersize=msize)
        ax[2].plot(rdate[::gap], r['AX'][::gap],marker=mk,markersize=msize,label = 'AY')
        ax[2].plot(rdate[::gap], r['AY'][::gap],marker=mk,markersize=msize,label = 'AX')
        ax[3].plot(rdate[::gap], r['TA'][::gap],marker=mk)
        ax[2].legend(frameon=False, ncol=3)
    elif flag == 'B':
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
    return 0

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
def show_counts(ax,axr,d,d4,d8,ns,match,onboard,msize=2):
    ddate = pd.to_datetime(d.iloc[:, 0])
    tsecs = ddate.iloc[-1].timestamp()-ddate.iloc[0].timestamp()
    
    t,call = counts(ax,np.array(d['Timestamp']),ddate,tsecs,markersize=msize, label = 'Trigger')
    counts(ax,np.array(d4['Timestamp']),ddate,tsecs,markersize=msize, label = 'Start')
    _,c8 = counts(ax,np.array(d8['Timestamp']),ddate,tsecs,markersize=msize, label = 'StartEnd')
    counts(ax,np.array(ns['Timestamp']),ddate,tsecs,markersize=msize, label = 'StartLim')
    _,cmatch = counts(ax,np.array(match['Timestamp']),ddate,tsecs,markersize=msize, label = 'SumMatch')
    counts(ax,np.array(match['Timestamp']),ddate,tsecs,markersize=msize, label = 'Onboard')
    ax.legend(frameon=False, ncol=3)
    ax.set_ylim(1,16*62*70)
    ax.axhline(16*62*60,color='r')
    
    #question mark
    axr.plot(t,c8/(call+0.0001),label = 'StartEnd/Trigger')
    axr.plot(t,cmatch/(c8+0.0001), label = 'SumMatch/StartEnd')
    axr.legend(frameon=False)
    return 0

# ================================ statistic ================================
def channel_count_ratio(df,c):
    if c == 'X1':
        count = ((df['X1_start (ns)'] != 0) | (df['X1_end (ns)'])).sum()
    if c == 'X2':
        count = ((df['X2_start (ns)'] != 0) | (df['X2_end (ns)'])).sum()
    if c == 'Y1':
        count = ((df['Y1_start (ns)'] != 0) | (df['Y1_end (ns)'])).sum()
    if c == 'Y2':
        count = ((df['Y2_start (ns)'] != 0) | (df['Y2_end (ns)'])).sum()
    return count     

def channel_count_ratio4(df):
    count_0 = ((df['X1_start (ns)'] == 0) & (df['X2_start (ns)'] == 0) & (df['Y1_start (ns)'] == 0) & (df['Y2_start (ns)'] == 0) & (df['X1_end (ns)'] == 0) & (df['X2_end (ns)'] == 0) & (df['Y1_end (ns)'] == 0) & (df['Y2_end (ns)'] == 0)).sum()
    count = ((df['X1_start (ns)'] != 0) & (df['X2_start (ns)'] != 0) & (df['Y1_start (ns)'] != 0) & (df['Y2_start (ns)'] != 0)).sum()
    print(count_0,len(df))
    return count
def channel_count_ratio8(df):
    count = ((df['X1_start (ns)'] != 0) & (df['X2_start (ns)'] != 0) & (df['Y1_start (ns)'] != 0) & (df['Y2_start (ns)'] != 0) & (df['X1_end (ns)'] != 0) & (df['X2_end (ns)'] != 0) & (df['Y1_end (ns)'] != 0) & (df['Y2_end (ns)'] != 0)).sum()
    return count 

def pspe(a,b):
    pe = np.roots([1,-(2*a+b),a])
    pe = pe[pe<1]
    ps = a/pe
    return ps,pe

def show_start_bar(df, prob, folder,ssid=1):
    phy.init_param(ssid)
    if prob == 'a':
        xlim = phy.A_NS_XLIM
        ylim = phy.A_NS_YLIM
    if prob == 'b':
        xlim = phy.B_NS_XLIM
        ylim = phy.B_NS_YLIM
    ybins = [1, ylim, ylim + 30, ylim + 30 *2, ylim + 30 * 3, ylim + 30 * 4, 300]
    ylabels = [str(int(y1)) + '-' + str(int(y2)) for y1,y2 in zip(ybins[:-1],ybins[1:])]
    xbins = [1, xlim, xlim + 50, xlim + 50 *2, xlim + 50 * 3, xlim + 50 * 4, xlim + 50 * 5, xlim + 50 * 6, 700]
    xlabels = [str(int(x1)) + '-' + str(int(x2)) for x1,x2 in zip(xbins[:-1],xbins[1:])]
    for i,cnl in enumerate(['X1_start (ns)', 'X2_start (ns)', 'Y1_start (ns)', 'Y2_start (ns)']):
        series = df[cnl]
        #zero_count = (series == -94*0.081).sum()
        zero_count = (series == 0).sum()
        filtered_series = series[series > 0]
        if i < 4:
            bins = xbins
            labels = xlabels
        else:
            bins = ybins
            labels = ylabels
        cut_data = pd.cut(filtered_series, bins=bins, labels=labels, right=False)
        interval_counts = cut_data.value_counts()
        # 重新索引，确保按照分箱顺序排列
        interval_counts = interval_counts.reindex(labels, fill_value=0)
        all_counts = pd.concat([pd.Series({'0': zero_count}), interval_counts])
        total_count = all_counts.sum()
        all_percentages = all_counts / total_count * 100

        f, ax = plt.subplots()
        # 绘制条形图
        #ars = ax.barh(all_counts.values,all_counts.index)
        bars = ax.barh(all_counts.index,all_counts.values)
        
        ymin, ymax = ax.get_ylim()
        y_range = ymax - ymin

        # 在每个条形的横向上添加百分比和计数标签
        for bar, (index, value) in zip(bars, all_counts.items()):
            # 获取条形的 y 坐标
            y_coord = bar.get_y() + bar.get_height() / 2
            # 计算归一化的 y 坐标
            normalized_y = (y_coord - ymin) / y_range
            # 添加计数标签
            ax.text(0.3, normalized_y, f'{value}', ha='left', va='center', transform=ax.transAxes)
            # 获取对应的百分比
            percentage = all_percentages[index]
            # 添加百分比标签
            ax.text(0.5, normalized_y, f'{percentage:.2f}%', ha='left', va='center', transform=ax.transAxes)
            ax.set_ylabel('Value Intervals')
            ax.set_xlabel('Count')
            ax.set_title(cnl)
            fig_w = 3.375
            f.set_size_inches(fig_w * 4, fig_w / 4 * 3)
            f.savefig(folder + cnl)


def show_8_dist(df, prob, folder,ssid):
    phy.init_param(ssid)
    if prob == 'a':
        xlim = phy.A_NS_XLIM
        ylim = phy.A_NS_YLIM
    if prob == 'b':
        xlim = phy.B_NS_XLIM
        ylim = phy.B_NS_YLIM
    ybins = [1, ylim, ylim + 40, ylim + 40 *2, ylim + 40 * 3, ylim + 40 * 4, ylim + 40 * 5, ylim + 40 * 6, ylim + 40 * 7, ylim + 40 * 8, 9999]
    ylabels = [str(int(y1)) + '-' + str(int(y2)) for y1,y2 in zip(ybins[:-1],ybins[1:])]
    xbins = [1, xlim, xlim + 50, xlim + 50 *2, xlim + 50 * 3, xlim + 50 * 4, xlim + 50 * 5, xlim + 50 * 6, xlim + 50 * 7, xlim + 50 * 8, 9999]
    xlabels = [str(int(x1)) + '-' + str(int(x2)) for x1,x2 in zip(xbins[:-1],xbins[1:])]
    for i,cnl in enumerate(['X1_start (ns)', 'X2_start (ns)', 'X1_end (ns)', 'X2_end (ns)', 'Y1_start (ns)', 'Y2_start (ns)','Y1_end (ns)', 'Y2_end (ns)']):
        series = df[cnl]
        #zero_count = (series == -94*0.081).sum()
        #zero_count = (series == 0).sum()
        #filtered_series = series[series > 0]
        if i < 4:
            bins = xbins
            labels = xlabels
        else:
            bins = ybins
            labels = ylabels
        cut_data = pd.cut(series, bins=bins, labels=labels, right=False)
        interval_counts = cut_data.value_counts()
        # 重新索引，确保按照分箱顺序排列
        interval_counts = interval_counts.reindex(labels, fill_value=0)
        #all_counts = pd.concat([pd.Series({'0': zero_count}), interval_counts])
        all_counts = interval_counts
        total_count = all_counts.sum()
        all_percentages = all_counts / total_count * 100

        f, ax = plt.subplots()
        # 绘制条形图
        #ars = ax.barh(all_counts.values,all_counts.index)
        bars = ax.barh(all_counts.index,all_counts.values)
        
        ymin, ymax = ax.get_ylim()
        y_range = ymax - ymin

        # 在每个条形的横向上添加百分比和计数标签
        for bar, (index, value) in zip(bars, all_counts.items()):
            # 获取条形的 y 坐标
            y_coord = bar.get_y() + bar.get_height() / 2
            # 计算归一化的 y 坐标
            normalized_y = (y_coord - ymin) / y_range
            # 添加计数标签
            ax.text(0.3, normalized_y, f'{value}', ha='left', va='center', transform=ax.transAxes)
            # 获取对应的百分比
            percentage = all_percentages[index]
            # 添加百分比标签
            ax.text(0.5, normalized_y, f'{percentage:.2f}%', ha='left', va='center', transform=ax.transAxes)
            ax.set_ylabel('Value Intervals')
            ax.set_xlabel('Count')
            ax.set_title(cnl)
            fig_w = 3.375
            f.set_size_inches(fig_w * 4, fig_w)
            f.savefig(folder + cnl)

def show_match_dist(df,folder='./',gap=10):
    sumx = df['X1_start (ns)'] +  df['X2_start (ns)']
    sumy = df['Y1_start (ns)'] +  df['Y2_start (ns)']
    ratiox, barx = np.histogram(sumx,np.arange(0,np.max(sumx),gap),density=True)
    f, ax = plt.subplots(2,1)
    ax[0].bar(barx[:-1],ratiox,gap)
    ax[0].axvline(45,color='r')
    ax[0].axvline(55,color='r')
    ax[0].set_ylabel('Sum X dist')
    ax[0].set_xlim(0,600)
    ax[0].set_ylabel('Sum X dist')
    ratioy, bary = np.histogram(sumx,np.arange(0,np.max(sumy),gap),density=True)
    ax[1].axvline(35,color='r')
    ax[1].axvline(45,color='r')
    ax[1].bar(bary[:-1],ratioy,gap)
    ax[1].set_xlabel('Sum (ns)')
    ax[1].set_ylabel('Sum Y dist')
    ax[1].set_xlim(0,600)
    fig_w = 3.375
    f.set_size_inches(fig_w * 4, fig_w*2)
    f.savefig(folder + 'match')

def show_4_dist(df, prob, folder,ssid=1):
    if prob == 'a':
        xlim = phy.A_NS_XLIM
        ylim = phy.A_NS_YLIM
    if prob == 'b':
        xlim = phy.B_NS_XLIM
        ylim = phy.B_NS_YLIM
    ybins = [1, ylim, ylim + 30, ylim + 30 *2, ylim + 30 * 3, ylim + 30 * 4, 300]
    ylabels = [str(int(y1)) + '-' + str(int(y2)) for y1,y2 in zip(ybins[:-1],ybins[1:])]
    xbins = [1, xlim, xlim + 50, xlim + 50 *2, xlim + 50 * 3, xlim + 50 * 4, xlim + 50 * 5, xlim + 50 * 6, 700]
    xlabels = [str(int(x1)) + '-' + str(int(x2)) for x1,x2 in zip(xbins[:-1],xbins[1:])]
    for i,cnl in enumerate(['X1_start (ns)', 'X2_start (ns)', 'Y1_start (ns)', 'Y2_start (ns)']):
        series = df[cnl]
        #zero_count = (series == -94*0.081).sum()
        zero_count = (series == 0).sum()
        filtered_series = series[series > 0]
        if i < 4:
            bins = xbins
            labels = xlabels
        else:
            bins = ybins
            labels = ylabels
        cut_data = pd.cut(filtered_series, bins=bins, labels=labels, right=False)
        interval_counts = cut_data.value_counts()
        # 重新索引，确保按照分箱顺序排列
        interval_counts = interval_counts.reindex(labels, fill_value=0)
        all_counts = pd.concat([pd.Series({'0': zero_count}), interval_counts])
        total_count = all_counts.sum()
        all_percentages = all_counts / total_count * 100

        f, ax = plt.subplots()
        # 绘制条形图
        #ars = ax.barh(all_counts.values,all_counts.index)
        bars = ax.barh(all_counts.index,all_counts.values)
        
        ymin, ymax = ax.get_ylim()
        y_range = ymax - ymin

        # 在每个条形的横向上添加百分比和计数标签
        for bar, (index, value) in zip(bars, all_counts.items()):
            # 获取条形的 y 坐标
            y_coord = bar.get_y() + bar.get_height() / 2
            # 计算归一化的 y 坐标
            normalized_y = (y_coord - ymin) / y_range
            # 添加计数标签
            ax.text(0.3, normalized_y, f'{value}', ha='left', va='center', transform=ax.transAxes)
            # 获取对应的百分比
            percentage = all_percentages[index]
            # 添加百分比标签
            ax.text(0.5, normalized_y, f'{percentage:.2f}%', ha='left', va='center', transform=ax.transAxes)
            ax.set_ylabel('Value Intervals')
            ax.set_xlabel('Count')
            ax.set_title(cnl)
            fig_w = 3.375
            f.set_size_inches(fig_w * 4, fig_w / 4 * 3)
            f.savefig(folder + cnl)


#quick
#plot from level 1
def show_start_end_level1(df,flag,snid=1,msize=1,s_or_hist=1):
    phy.init_param(snid)
    if flag == 'a':
        xc = phy.A_X_CENTER_NS
        yc = phy.A_Y_CENTER_NS
        xcoef = phy.A_X_COF
        ycoef = phy.A_Y_COF
    if flag == 'b':
        xc = phy.B_X_CENTER_NS
        yc = phy.B_Y_CENTER_NS
        xcoef = phy.B_X_COF
        ycoef = phy.B_Y_COF
        
    f,(ax1,ax2)=plt.subplots(2,1)
    x1 = (df['X2_start (ns)'] - df['X1_start (ns)'])
    y1 = (df['Y2_start (ns)'] - df['Y1_start (ns)'])
    x2 = (df['X2_end (ns)'] - df['X1_end (ns)'])
    y2 = (df['Y2_end (ns)'] - df['Y1_end (ns)'])

    x1 =  (x1 - xc) * xcoef
    y1 =  (y1 - yc) * ycoef
    x2 =  (x2 - xc) * xcoef
    y2 =  (y2 - yc) * ycoef
    if s_or_hist == 0:
        ax1.scatter(x1,y1, label = 'start',s=msize)
        ax2.scatter(x2,y2, label = 'end',s=msize)
        ax1.legend()
        ax2.legend()
    elif s_or_hist == 1:
        hh = ax1.hist2d(x1,y1,np.arange(-100,100,10))
        f.colorbar(hh[3], ax=ax1)
        hh = ax2.hist2d(x2,y2,np.arange(-100,100,10))
        f.colorbar(hh[3], ax=ax2)

    aux.draw_rectangle(ax1,90,30)
    aux.draw_rectangle(ax2,90,30)
    aux.draw_rectangle(ax1,150,60)
    aux.draw_rectangle(ax2,150,60)
    ax1.set_xlabel('X (mm)')
    ax1.set_ylabel('Y (mm)')
    ax2.set_xlabel('X (mm)')
    ax2.set_ylabel('Y (mm)')
    return f,(ax1,ax2)

#level2
def show_start_end_level2(df,msize=1):
    f,(ax1,ax2)=plt.subplots(2,1)
    ax1.scatter(df['X1 (mm)'], df['Y1 (mm)'], label = 'start',markersize=msize)
    ax2.scatter(df['X2 (mm)'], df['Y2 (mm)'], label = 'end',markersize=msize)
    ax1.legend()
    ax2.legend()
    aux.draw_rectangle(ax1,90,30)
    aux.draw_rectangle(ax1,150,60)
    aux.draw_rectangle(ax2,90,30)
    aux.draw_rectangle(ax2,150,60)
    ax1.set_xlabel('X (mm)')
    ax1.set_ylabel('Y (mm)')
    ax2.set_xlabel('X (mm)')
    ax2.set_ylabel('Y (mm)')
    return f,(ax1,ax2)


if __name__ == '__main__':
    #Level 1 A
    df_1A = csv0.wrap_df("./L1_a.csv")
    df_1A_1 = csv0.df_slices_time(df_1A,'2024-10-01 00:00:00','2024-10-15 00:00:00')

    #phase 1
    savepath='phase1_a/'
    csv0.folder(savepath)
    csv0.show_channel_bar(df_1A_1,'x','X1_start (ns)',savepath)
    csv0.show_channel_bar(df_1A_1,'x','X1_end (ns)',savepath)
    csv0.show_channel_bar(df_1A_1,'x','X2_start (ns)',savepath)
    csv0.show_channel_bar(df_1A_1,'x','X2_end (ns)',savepath)
    csv0.show_channel_bar(df_1A_1,'y','Y1_start (ns)',savepath)
    csv0.show_channel_bar(df_1A_1,'y','Y1_end (ns)',savepath)
    csv0.show_channel_bar(df_1A_1,'y','Y2_start (ns)',savepath)
    csv0.show_channel_bar(df_1A_1,'y','Y2_end (ns)',savepath)


    #df_1A_2 = df_slices_time(df_1A,'2024-10-15 00:00:00','2024-10-30 00:00:00')
#
    ##phase 2
    #savepath='phase2_a/'
    #folder(savepath)
    #show_channel_bar(df_1A_2,'x','X1_start (ns)',savepath)
    #show_channel_bar(df_1A_2,'x','X1_end (ns)',savepath)
    #show_channel_bar(df_1A_2,'x','X2_start (ns)',savepath)
    #show_channel_bar(df_1A_2,'x','X2_end (ns)',savepath)
#
    #show_channel_bar(df_1A_2,'y','Y1_start (ns)',savepath)
    #show_channel_bar(df_1A_2,'y','Y1_end (ns)',savepath)
    #show_channel_bar(df_1A_2,'y','Y2_start (ns)',savepath)
    #show_channel_bar(df_1A_2,'y','Y2_end (ns)',savepath)
#
    #df_1A_3 = df_slices_time(df_1A,'2024-11-01 00:00:00','2024-11-15 00:00:00')
#
    ##phase 1
    #savepath='phase3_a/'
    #folder(savepath)
    #show_channel_bar(df_1A_3,'x','X1_start (ns)',savepath)
    #show_channel_bar(df_1A_3,'x','X1_end (ns)',savepath)
    #show_channel_bar(df_1A_3,'x','X2_start (ns)',savepath)
    #show_channel_bar(df_1A_3,'x','X2_end (ns)',savepath)
#
    #show_channel_bar(df_1A_3,'y','Y1_start (ns)',savepath)
    #show_channel_bar(df_1A_3,'y','Y1_end (ns)',savepath)
    #show_channel_bar(df_1A_3,'y','Y2_start (ns)',savepath)
    #show_channel_bar(df_1A_3,'y','Y2_end (ns)',savepath)
#
#
    #df_1A_4_work = df_slices_time(df_1A,'2024-12-03 15:47:00','2024-12-03 15:47:30')
#
    #savepath='phase4_a/'
    #folder(savepath)
    #show_channel_bar(df_1A_4_work,'x','X1_start (ns)',savepath)
    #show_channel_bar(df_1A_4_work,'x','X1_end (ns)',savepath)
    #show_channel_bar(df_1A_4_work,'x','X2_start (ns)',savepath)
    #show_channel_bar(df_1A_4_work,'x','X2_end (ns)',savepath)
    #show_channel_bar(df_1A_4_work,'y','Y1_start (ns)',savepath)
    #show_channel_bar(df_1A_4_work,'y','Y1_end (ns)',savepath)
    #show_channel_bar(df_1A_4_work,'y','Y2_start (ns)',savepath)
    #show_channel_bar(df_1A_4_work,'y','Y2_end (ns)',savepath)
#
#
#
    ##Level 1 A
    #df_1B = wrap_df("./L1_b.csv")
    #df_1B_1 = df_slices_time(df_1B,'2024-10-01 00:00:00','2024-10-15 00:00:00')
#
    ##phase 1
    #savepath='phase1_b/'
    #folder(savepath)
    #show_channel_bar(df_1B_1,'x','X1_start (ns)',savepath)
    #show_channel_bar(df_1B_1,'x','X1_end (ns)',savepath)
    #show_channel_bar(df_1B_1,'x','X2_start (ns)',savepath)
    #show_channel_bar(df_1B_1,'x','X2_end (ns)',savepath)
#
    #show_channel_bar(df_1B_1,'y','Y1_start (ns)',savepath)
    #show_channel_bar(df_1B_1,'y','Y1_end (ns)',savepath)
    #show_channel_bar(df_1B_1,'y','Y2_start (ns)',savepath)
    #show_channel_bar(df_1B_1,'y','Y2_end (ns)',savepath)
#
#
    #df_1B_2 = df_slices_time(df_1B,'2024-10-15 00:00:00','2024-10-30 00:00:00')
#
    ##phase 2
    #savepath='phase2_b/'
    #folder(savepath)
    #show_channel_bar(df_1B_2,'x','X1_start (ns)',savepath)
    #show_channel_bar(df_1B_2,'x','X1_end (ns)',savepath)
    #show_channel_bar(df_1B_2,'x','X2_start (ns)',savepath)
    #show_channel_bar(df_1B_2,'x','X2_end (ns)',savepath)
#
    #show_channel_bar(df_1B_2,'y','Y1_start (ns)',savepath)
    #show_channel_bar(df_1B_2,'y','Y1_end (ns)',savepath)
    #show_channel_bar(df_1B_2,'y','Y2_start (ns)',savepath)
    #show_channel_bar(df_1B_2,'y','Y2_end (ns)',savepath)
#
    #df_1B_3 = df_slices_time(df_1B,'2024-11-01 00:00:00','2024-11-15 00:00:00')
#
    ##phase 1
    #savepath='phase3_b/'
    #folder(savepath)
    #show_channel_bar(df_1B_3,'x','X1_start (ns)',savepath)
    #show_channel_bar(df_1B_3,'x','X1_end (ns)',savepath)
    #show_channel_bar(df_1B_3,'x','X2_start (ns)',savepath)
    #show_channel_bar(df_1B_3,'x','X2_end (ns)',savepath)
#
    #show_channel_bar(df_1B_3,'y','Y1_start (ns)',savepath)
    #show_channel_bar(df_1B_3,'y','Y1_end (ns)',savepath)
    #show_channel_bar(df_1B_3,'y','Y2_start (ns)',savepath)
    #show_channel_bar(df_1B_3,'y','Y2_end (ns)',savepath)
#
#
    #df_2A = wrap_df("./L2_a.csv")
    #df_2A_1 = df_slices_time(df_2A,'2024-10-01 00:00:00','2024-10-15 00:00:00')
    #df_2A_2 = df_slices_time(df_2A,'2024-10-15 00:00:00','2024-11-01 00:00:00')
    #df_2A_3 = df_slices_time(df_2A,'2024-11-01 00:00:00','2024-11-15 00:00:00')
#
    #df_2B = wrap_df("./L2_b.csv")
    #df_2B_1 = df_slices_time(df_2B,'2024-10-01 00:00:00','2024-10-15 00:00:00')
    #df_2B_2 = df_slices_time(df_2B,'2024-10-15 00:00:00','2024-11-01 00:00:00')
    #df_2B_3 = df_slices_time(df_2B,'2024-11-01 00:00:00','2024-11-15 00:00:00')
#
