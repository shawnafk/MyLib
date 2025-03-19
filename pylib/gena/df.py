import pandas as pd
import numpy as np
# df has shape

class ENA_DataFrame(pd.DataFrame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #change Date string to datetime, when init
        self['Date'] = pd.to_datetime(self['Date'])
    
    def split_time(self,t1,t2):
        s1 = pd.to_datetime(t1)
        s2 = pd.to_datetime(t2)
        df_date = pd.to_datetime(self['Date'])
        return ENA_DataFrame(self[(df_date > s1) & (df_date < s2)])

    def separate_df(self,tgap=600):
        time_diff = self['Date'].diff()
        split_indices = time_diff[time_diff > pd.Timedelta(seconds=tgap)].index
        experiments = []
        time_interval = []
        start_index = 0
        for split_index in split_indices:
            time_interval.append([self.iloc[start_index,0],self.iloc[split_index,0]])
            experiment = ENA_DataFrame(self[start_index:split_index])
            experiments.append(experiment)
            start_index = split_index
        
        last_experiment = self[start_index:]
        experiments.append(last_experiment)
        return experiments, time_interval

def read_csv(csv_path,*args, **kwargs):
    df = pd.read_csv(csv_path,*args, **kwargs).drop_duplicates().sort_values(by='Timestamp')
    edf = ENA_DataFrame(df)
    return edf
    
def separate_dfs(another_df, date_ranges):
    another_df.iloc[:, 0] = pd.to_datetime(another_df.iloc[:, 0])
    sub_df = []
    for i, (start_date, end_date) in enumerate(date_ranges):
        print(start_date, end_date)
        # 筛选出在当前日期区间内的数据
        sub_df.append(another_df[(another_df.iloc[:, 0] >= start_date) & (another_df.iloc[:, 0] <= end_date)])
    return sub_df

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
