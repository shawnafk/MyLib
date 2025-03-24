import numpy as np
from . import physical as phy
from . import df as gdf

#select data with some conditions
#return len and sub dataframe
def select_data(df,cond):
    sub_df = df[cond]
    len_sub = len(sub_df)
    return len_sub, sub_df


#bad condition

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

def match_eff(df,xl=45,xu=55,yl=35,yu=45):
    sumx =  np.array((df['X1_start (ns)'] + df['X2_start (ns)']))
    sumy =  np.array((df['Y1_start (ns)'] + df['Y2_start (ns)']))
    cond = (sumx<xu) & (sumx>xl) & (sumy<yu) & (sumy>yl)
    len,sub_df = select_data(df,cond)
    return len,sub_df


    
#level 2
def l1_2_l2(df,flag,snid):
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
    x1 = (df['X2_start (ns)'] - df['X1_start (ns)'])
    y1 = (df['Y2_start (ns)'] - df['Y1_start (ns)'])
    x2 = (df['X2_end (ns)'] - df['X1_end (ns)'])
    y2 = (df['Y2_end (ns)'] - df['Y1_end (ns)'])

    x1 =  (x1 - xc) * xcoef
    y1 =  (y1 - yc) * ycoef
    x2 =  (x2 - xc) * xcoef
    y2 =  (y2 - yc) * ycoef
    
    return gdf.ENA_DataFrame({'Date': df['Date'], 'X1 (mm)': x1, 'Y1 (mm)': y1, 'X2 (mm)': x2, 'Y2 (mm)': y2})

def onboard_eff(df):
    cond = (df['X1 (mm)']<90) & (df['Y1 (mm)']<30) & (df['X2 (mm)']<150) & (df['Y2 (mm)']<60)
    len,sub_df = select_data(df,cond)
    return len,sub_df
    

