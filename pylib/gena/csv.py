import pandas as pd
import numpy as np
from datetime import datetime
def dt_from_ts(ts):
    dt = datetime.strptime("21-01-01 00:00:00", "%y-%m-%d %H:%M:%S").timestamp()
    return  [datetime.fromtimestamp(timestamp+dt) for timestamp in ts]
    
def seconds_to_reference_time(target_time_str):
    # 将字符串格式的时间转换为datetime对象
    target_time = datetime.strptime(target_time_str, '%Y-%m-%d %H:%M:%S')
    reference_time = datetime(2021, 1, 1, 0, 0, 0)
    # 计算时间差
    time_diff = target_time - reference_time
    # 返回时间差的秒数
    return time_diff.total_seconds()

def df_slices_time_by(df,t1,t2,by='Satellite seconds'):
    s1 = seconds_to_reference_time(t1)
    s2 = seconds_to_reference_time(t2)
    # is it correct to use & here?
    #return df[(df['TimeStamp'] > s1) & (df['TimeStamp'] < s2)]
    return df[(df[by] > s1) & (df[by] < s2)]

def wrap_df(fn):
    return  pd.read_csv(fn,encoding='gbk',header=0).iloc[:,:]


from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 18})
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

def folder(fn):
    import os
    folder_path = fn
    try:
        os.makedirs(folder_path)
        print(f"文件夹 {folder_path} 创建成功")
    except FileExistsError:
        print(f"文件夹 {folder_path} 已经存在")
    except Exception as e:
        print(f"创建文件夹时出错: {e}")
