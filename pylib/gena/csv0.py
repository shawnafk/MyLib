import pandas as pd
import numpy as np
# df has shape 
def df_t1_t2(df,t1,t2):
    s1 = pd.to_datetime(t1)
    s2 = pd.to_datetime(t2)
    df_date = pd.to_datetime(df['Date'])
    return df[(df_date > s1) & (df_date < s2)]

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
