import matplotlib.pyplot as plt
from myplot import aux
import numpy as np
from gena import csv
import pandas as pd
X_MCP = 145
Y_MCP = 55
#ns limit need divided by 2
A_X_LEFT_NS  = -50.44742
A_X_RIGHT_NS = 54.91726

A_Y_LEFT_NS  = -34.71570
A_Y_RIGHT_NS = 23.59747

B_X_LEFT_NS  = -51.18360
B_X_RIGHT_NS = 53.40479

#SN_ID 1 orbit 2 ground
SN_ID = 1
if SN_ID == 1:
    #B_Y_LEFT_NS  = -16.68598
    #B_Y_RIGHT_NS = 41.62719
    B_Y_LEFT_NS  = -34.71570
    B_Y_RIGHT_NS = 23.59747

if SN_ID == 2:
    B_Y_LEFT_NS  = -34.71570
    B_Y_RIGHT_NS = 23.59747

A_X_MAX_NS = (abs(A_X_RIGHT_NS) + abs(A_X_LEFT_NS));
A_Y_MAX_NS = (abs(A_Y_RIGHT_NS) + abs(A_Y_LEFT_NS));
#need divided by 2 and yes
A_NS_XLIM = A_X_MAX_NS/2;
A_NS_YLIM = A_Y_MAX_NS/2;
A_X_CENTER_NS = (A_X_LEFT_NS + A_X_RIGHT_NS);
A_Y_CENTER_NS = (A_Y_LEFT_NS + A_Y_RIGHT_NS);
A_X_COF = X_MCP/A_X_MAX_NS;
A_Y_COF = Y_MCP/A_Y_MAX_NS;

B_X_MAX_NS = abs(B_X_RIGHT_NS) + abs(B_X_LEFT_NS);
B_Y_MAX_NS = abs(B_Y_RIGHT_NS) + abs(B_Y_LEFT_NS);
B_NS_XLIM = B_X_MAX_NS/2;
B_NS_YLIM = B_Y_MAX_NS/2;
#should divided by 2 but not yet
B_X_CENTER_NS = (B_X_LEFT_NS + B_X_RIGHT_NS);
B_Y_CENTER_NS = (B_Y_LEFT_NS + B_Y_RIGHT_NS);
B_X_COF = X_MCP/B_X_MAX_NS;
B_Y_COF = Y_MCP/B_Y_MAX_NS;

def show_start_end1(df,flag):
    if flag == 'a':
        xc = A_X_CENTER_NS
        yc = A_Y_CENTER_NS
        xcoef = A_X_COF
        ycoef = A_Y_COF
    if flag == 'b':
        xc = B_X_CENTER_NS
        yc = B_Y_CENTER_NS
        xcoef = B_X_COF
        ycoef = B_Y_COF
    f,(ax1,ax2)=plt.subplots(2,1)
    x1 = (df['X2_start (ns)'] - df['X1_start (ns)'])
    y1 = (df['Y2_start (ns)'] - df['Y1_start (ns)'])
    x2 = (df['X2_end (ns)'] - df['X1_end (ns)'])
    y2 = (df['Y2_end (ns)'] - df['Y1_end (ns)'])

    x1 =  (x1 - xc) * xcoef
    y1 =  (y1 - yc) * ycoef
    x2 =  (x2 - xc) * xcoef
    y2 =  (y2 - yc) * ycoef
    ax1.scatter(x1,y1, label = 'start',s=1)
    ax2.scatter(x2,y2, label = 'end',s=1)
    ax1.legend()
    ax2.legend()
    aux.draw_rectangle(ax1,90,30)
    aux.draw_rectangle(ax2,150,60)
    ax1.set_xlabel('X (mm)')
    ax1.set_ylabel('Y (mm)')
    ax2.set_xlabel('X (mm)')
    ax2.set_ylabel('Y (mm)')
    return f

#level2
def show_start_end(df):
    f,(ax1,ax2)=plt.subplots(2,1)
    ax1.scatter(df['X1 (mm)'], df['Y1 (mm)'], label = 'start',s=1)
    ax2.scatter(df['X2 (mm)'], df['Y2 (mm)'], label = 'end',s=1)
    ax1.legend()
    ax2.legend()
    aux.draw_rectangle(ax1,90,30)
    aux.draw_rectangle(ax2,150,60)
    ax1.set_xlabel('X (mm)')
    ax1.set_ylabel('Y (mm)')
    ax2.set_xlabel('X (mm)')
    ax2.set_ylabel('Y (mm)')
    return f

#level 1A
def ratio_levels(df):
    sum = np.array((df['X1_start (ns)'] * df['X2_start (ns)'] * df['Y1_start (ns)'] * df['Y2_start (ns)'] * df['X1_end (ns)'] * df['X2_end (ns)'] * df['Y1_end (ns)'] * df['Y2_end (ns)']))
    df1B = df[sum!=0]
    all = sum.shape[0]
    nonzero = np.where(sum!=0)[0].shape[0]
    
    df1C = df1B[(df1B['X1_start (ns)']<53) & (df1B['X2_start (ns)']<53) & (df1B['Y1_start (ns)']<40) & (df1B['Y2_start (ns)']<40)]
    nseff = df1C.shape[0]
    
    sumx =  np.array((df1C['X1_start (ns)'] + df1C['X2_start (ns)']))
    sumy =  np.array((df1C['Y1_start (ns)'] + df1C['Y2_start (ns)']))
    
    df2 = df1C[(sumx<55) & (sumx>45) & (sumy<45) & (sumy>35)]
    match_eff = df2.shape[0]
    return [nonzero, all, nonzero/all],[nseff,nonzero,nseff/nonzero],[match_eff,nseff,match_eff/nseff]


def ratio_new(df,prob):
    if prob == 'A':
        xlim = A_NS_XLIM
        ylim = A_NS_YLIM
    if prob == 'B':
        xlim = B_NS_XLIM
        ylim = B_NS_YLIM
    sum = np.array((df['X1_start (ns)'] * df['X2_start (ns)'] * df['Y1_start (ns)'] * df['Y2_start (ns)'] ))
    df1B = df[sum!=0]
    df1B = df1B[(df1B['X1_start (ns)']<xlim) & (df1B['X2_start (ns)']<xlim) & (df1B['Y1_start (ns)']<ylim) & (df1B['Y2_start (ns)']<ylim)]
    all = len(df)
    len1B = len(df1B)
    sum = np.array((df1B['X1_end (ns)'] * df1B['X2_end (ns)'] * df1B['Y1_end (ns)'] * df1B['Y2_end (ns)'] ))
    df1C = df1B[sum!=0]
    len1C = len(df1C)
    
    sumx =  np.array((df1C['X1_start (ns)'] + df1C['X2_start (ns)']))
    sumy =  np.array((df1C['Y1_start (ns)'] + df1C['Y2_start (ns)']))
    
    df1D = df1C[(sumx<55) & (sumx>45) & (sumy<45) & (sumy>35)]
    len1D = len(df1D)
    return [all, len1B, len1C,len1D]

def ratio_l2(df):
    df_rational= df[(df['X1 (mm)']<90) & (df['Y1 (mm)']<30) & (df['X2 (mm)']<150) & (df['Y2 (mm)']<60)]
    num_rational = df_rational.shape[0]
    all=df.shape[0]
    return [num_rational,all,num_rational/all]

def channel_count_rate(df,c):
    if c == 'X1':
        count = ((df['X1_start (ns)'] != 0) | (df['X1_end (ns)'])).sum()
    if c == 'X2':
        count = ((df['X2_start (ns)'] != 0) | (df['X2_end (ns)'])).sum()
    if c == 'Y1':
        count = ((df['Y1_start (ns)'] != 0) | (df['Y1_end (ns)'])).sum()
    if c == 'Y2':
        count = ((df['Y2_start (ns)'] != 0) | (df['Y2_end (ns)'])).sum()
    return count     

def channel_count_rate4(df):
    count_0 = ((df['X1_start (ns)'] == 0) & (df['X2_start (ns)'] == 0) & (df['Y1_start (ns)'] == 0) & (df['Y2_start (ns)'] == 0) & (df['X1_end (ns)'] == 0) & (df['X2_end (ns)'] == 0) & (df['Y1_end (ns)'] == 0) & (df['Y2_end (ns)'] == 0)).sum()
    count = ((df['X1_start (ns)'] != 0) & (df['X2_start (ns)'] != 0) & (df['Y1_start (ns)'] != 0) & (df['Y2_start (ns)'] != 0)).sum()
    print(count_0,len(df))
    return count
def channel_count_rate8(df):
    count = ((df['X1_start (ns)'] != 0) & (df['X2_start (ns)'] != 0) & (df['Y1_start (ns)'] != 0) & (df['Y2_start (ns)'] != 0) & (df['X1_end (ns)'] != 0) & (df['X2_end (ns)'] != 0) & (df['Y1_end (ns)'] != 0) & (df['Y2_end (ns)'] != 0)).sum()
    return count 

def pspe(a,b):
    pe = np.roots([1,-(2*a+b),a])
    pe = pe[pe<1]
    ps = a/pe
    return ps,pe

def show_start_bar(df, prob, folder):
    if prob == 'a':
        xlim = A_NS_XLIM
        ylim = A_NS_YLIM
    if prob == 'b':
        xlim = B_NS_XLIM
        ylim = B_NS_YLIM
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


def show_1c_bar(df, prob, folder):
    if prob == 'a':
        xlim = A_NS_XLIM
        ylim = A_NS_YLIM
    if prob == 'b':
        xlim = B_NS_XLIM
        ylim = B_NS_YLIM
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



if __name__ == '__main__':
    #Level 1 A
    df_1A = csv.wrap_df("./L1_a.csv")
    df_1A_1 = csv.df_slices_time(df_1A,'2024-10-01 00:00:00','2024-10-15 00:00:00')

    #phase 1
    savepath='phase1_a/'
    csv.folder(savepath)
    csv.show_channel_bar(df_1A_1,'x','X1_start (ns)',savepath)
    csv.show_channel_bar(df_1A_1,'x','X1_end (ns)',savepath)
    csv.show_channel_bar(df_1A_1,'x','X2_start (ns)',savepath)
    csv.show_channel_bar(df_1A_1,'x','X2_end (ns)',savepath)
    csv.show_channel_bar(df_1A_1,'y','Y1_start (ns)',savepath)
    csv.show_channel_bar(df_1A_1,'y','Y1_end (ns)',savepath)
    csv.show_channel_bar(df_1A_1,'y','Y2_start (ns)',savepath)
    csv.show_channel_bar(df_1A_1,'y','Y2_end (ns)',savepath)


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