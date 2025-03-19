import pandas as pd
def separate_df(df,tgap=600):
    # 计算相邻行的时间间隔
    df.iloc[:, 0] = pd.to_datetime(df.iloc[:, 0])
    time_diff = df.iloc[:, 0].diff()
    # 找出时间间隔大于10秒的行的索引
    split_indices = time_diff[time_diff > pd.Timedelta(seconds=tgap)].index
    # 初始化实验列表
    experiments = []
    time_interval = []
    # 起始索引
    start_index = 0
    # 遍历分割索引
    for split_index in split_indices:
        # 提取当前实验的数据
        time_interval.append([df.iloc[start_index,0],df.iloc[split_index,0]])
        experiment = df[start_index:split_index]
        experiments.append(experiment)
        # 更新起始索引
        start_index = split_index
    # 提取最后一个实验的数据
    last_experiment = df[start_index:]
    experiments.append(last_experiment)
    return experiments, time_interval

def separate_dfs(another_df, date_ranges):
    another_df.iloc[:, 0] = pd.to_datetime(another_df.iloc[:, 0])
    sub_df = []
    for i, (start_date, end_date) in enumerate(date_ranges):
        print(start_date, end_date)
        # 筛选出在当前日期区间内的数据
        sub_df.append(another_df[(another_df.iloc[:, 0] >= start_date) & (another_df.iloc[:, 0] <= end_date)])
    return sub_df


