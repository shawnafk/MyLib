import numpy as np
import pandas as pd

#select data with some conditions
#return len and sub dataframe
def select_data(df,cond):
    sub_df = df[cond]
    len_sub = len(sub_df)
    return len_sub, sub_df

#condition
def count_eff(df): 
    sum = np.array((df['X1_start (ns)'] + df['X2_start (ns)'] + df['Y1_start (ns)'] + df['Y2_start (ns)'] + df['X1_end (ns)'] + df['X2_end (ns)'] + df['Y1_end (ns)'] + df['Y2_end (ns)']))
    cond = (sum!=0)
    len,sub_df = select_data(df,cond)
    return len,sub_df

def start_eff(df): 
    multi = np.array((df['X1_start (ns)'] * df['X2_start (ns)'] * df['Y1_start (ns)'] * df['Y2_start (ns)']))
    cond = (multi!=0)
    len,sub_df = select_data(df,cond)
    return len,sub_df

def start_stop_eff(df): 
    multi = np.array((df['X1_start (ns)'] * df['X2_start (ns)'] * df['Y1_start (ns)'] * df['Y2_start (ns)'] * df['X1_end (ns)'] * df['X2_end (ns)'] * df['Y1_end (ns)'] * df['Y2_end (ns)']))
    cond = (multi!=0)
    len,sub_df = select_data(df,cond)
    return len,sub_df

def start_ns_eff(df,xlim,ylim):
    cond = (df['X1_start (ns)']<xlim) & (df['X2_start (ns)']) & (df['Y1_start (ns)']) & (df['Y2_start (ns)']<ylim)
    len,sub_df = select_data(df,cond)
    return len,sub_df

def match_eff(df):
    sumx =  np.array((df['X1_start (ns)'] + df['X2_start (ns)']))
    sumy =  np.array((df['Y1_start (ns)'] + df['Y2_start (ns)']))
    cond = (sumx<55) & (sumx>45) & (sumy<45) & (sumy>35)
    len,sub_df = select_data(df,cond)
    return len,sub_df
    
#level 2
def onboard_eff(df):
    cond = (df['X1 (mm)']<90) & (df['Y1 (mm)']<30) & (df['X2 (mm)']<150) & (df['Y2 (mm)']<60)
    len,sub_df = select_data(df,cond)
    return len,sub_df
    