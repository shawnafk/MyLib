import matplotlib.pyplot as plt
from myplot import aux
import numpy as np
from gena import csv
def show_start_end(df):
    f,(ax1,ax2)=plt.subplots(2,1)
    ax1.scatter(df['X1 (mm)']/2, df['Y1 (mm)']/2, label = 'start',s=1)
    ax2.scatter(df['X2 (mm)']/2, df['Y2 (mm)']/2, label = 'end',s=1)
    ax1.legend()
    ax2.legend()
    aux.draw_rectangle(ax1,90,30)
    aux.draw_rectangle(ax2,150,60)
    ax1.set_xlabel('X (mm)')
    ax1.set_ylabel('Y (mm)')
    ax2.set_xlabel('X (mm)')
    ax2.set_ylabel('Y (mm)')

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


def ratio_l2(df):
    df_rational= df[(df['X1 (mm)']<90) & (df['Y1 (mm)']<30) & (df['X2 (mm)']<150) & (df['Y2 (mm)']<60)]
    num_rational = df_rational.shape[0]
    all=df.shape[0]
    return [num_rational,all,num_rational/all]

def pspe(a,b):
    pe = np.roots([1,-(2*a+b),a])
    pe = pe[pe<1]
    ps = a/pe
    return ps,pe
#Level 1 A
df_1A = gena_csv.wrap_df("./L1_a.csv")
df_1A_1 = df_slices_time(df_1A,'2024-10-01 00:00:00','2024-10-15 00:00:00')

#phase 1
savepath='phase1_a/'
folder(savepath)
show_channel_bar(df_1A_1,'x','X1_start (ns)',savepath)
show_channel_bar(df_1A_1,'x','X1_end (ns)',savepath)
show_channel_bar(df_1A_1,'x','X2_start (ns)',savepath)
show_channel_bar(df_1A_1,'x','X2_end (ns)',savepath)
                         
show_channel_bar(df_1A_1,'y','Y1_start (ns)',savepath)
show_channel_bar(df_1A_1,'y','Y1_end (ns)',savepath)
show_channel_bar(df_1A_1,'y','Y2_start (ns)',savepath)
show_channel_bar(df_1A_1,'y','Y2_end (ns)',savepath)


df_1A_2 = df_slices_time(df_1A,'2024-10-15 00:00:00','2024-10-30 00:00:00')

#phase 2
savepath='phase2_a/'
folder(savepath)
show_channel_bar(df_1A_2,'x','X1_start (ns)',savepath)
show_channel_bar(df_1A_2,'x','X1_end (ns)',savepath)
show_channel_bar(df_1A_2,'x','X2_start (ns)',savepath)
show_channel_bar(df_1A_2,'x','X2_end (ns)',savepath)
                     
show_channel_bar(df_1A_2,'y','Y1_start (ns)',savepath)
show_channel_bar(df_1A_2,'y','Y1_end (ns)',savepath)
show_channel_bar(df_1A_2,'y','Y2_start (ns)',savepath)
show_channel_bar(df_1A_2,'y','Y2_end (ns)',savepath)

df_1A_3 = df_slices_time(df_1A,'2024-11-01 00:00:00','2024-11-15 00:00:00')

#phase 1
savepath='phase3_a/'
folder(savepath)
show_channel_bar(df_1A_3,'x','X1_start (ns)',savepath)
show_channel_bar(df_1A_3,'x','X1_end (ns)',savepath)
show_channel_bar(df_1A_3,'x','X2_start (ns)',savepath)
show_channel_bar(df_1A_3,'x','X2_end (ns)',savepath)
                     
show_channel_bar(df_1A_3,'y','Y1_start (ns)',savepath)
show_channel_bar(df_1A_3,'y','Y1_end (ns)',savepath)
show_channel_bar(df_1A_3,'y','Y2_start (ns)',savepath)
show_channel_bar(df_1A_3,'y','Y2_end (ns)',savepath)


df_1A_4_work = df_slices_time(df_1A,'2024-12-03 15:47:00','2024-12-03 15:47:30')

savepath='phase4_a/'
folder(savepath)
show_channel_bar(df_1A_4_work,'x','X1_start (ns)',savepath)
show_channel_bar(df_1A_4_work,'x','X1_end (ns)',savepath)
show_channel_bar(df_1A_4_work,'x','X2_start (ns)',savepath)
show_channel_bar(df_1A_4_work,'x','X2_end (ns)',savepath)
show_channel_bar(df_1A_4_work,'y','Y1_start (ns)',savepath)
show_channel_bar(df_1A_4_work,'y','Y1_end (ns)',savepath)
show_channel_bar(df_1A_4_work,'y','Y2_start (ns)',savepath)
show_channel_bar(df_1A_4_work,'y','Y2_end (ns)',savepath)



#Level 1 A
df_1B = wrap_df("./L1_b.csv")
df_1B_1 = df_slices_time(df_1B,'2024-10-01 00:00:00','2024-10-15 00:00:00')

#phase 1
savepath='phase1_b/'
folder(savepath)
show_channel_bar(df_1B_1,'x','X1_start (ns)',savepath)
show_channel_bar(df_1B_1,'x','X1_end (ns)',savepath)
show_channel_bar(df_1B_1,'x','X2_start (ns)',savepath)
show_channel_bar(df_1B_1,'x','X2_end (ns)',savepath)
                         
show_channel_bar(df_1B_1,'y','Y1_start (ns)',savepath)
show_channel_bar(df_1B_1,'y','Y1_end (ns)',savepath)
show_channel_bar(df_1B_1,'y','Y2_start (ns)',savepath)
show_channel_bar(df_1B_1,'y','Y2_end (ns)',savepath)


df_1B_2 = df_slices_time(df_1B,'2024-10-15 00:00:00','2024-10-30 00:00:00')

#phase 2
savepath='phase2_b/'
folder(savepath)
show_channel_bar(df_1B_2,'x','X1_start (ns)',savepath)
show_channel_bar(df_1B_2,'x','X1_end (ns)',savepath)
show_channel_bar(df_1B_2,'x','X2_start (ns)',savepath)
show_channel_bar(df_1B_2,'x','X2_end (ns)',savepath)
                     
show_channel_bar(df_1B_2,'y','Y1_start (ns)',savepath)
show_channel_bar(df_1B_2,'y','Y1_end (ns)',savepath)
show_channel_bar(df_1B_2,'y','Y2_start (ns)',savepath)
show_channel_bar(df_1B_2,'y','Y2_end (ns)',savepath)

df_1B_3 = df_slices_time(df_1B,'2024-11-01 00:00:00','2024-11-15 00:00:00')

#phase 1
savepath='phase3_b/'
folder(savepath)
show_channel_bar(df_1B_3,'x','X1_start (ns)',savepath)
show_channel_bar(df_1B_3,'x','X1_end (ns)',savepath)
show_channel_bar(df_1B_3,'x','X2_start (ns)',savepath)
show_channel_bar(df_1B_3,'x','X2_end (ns)',savepath)
                     
show_channel_bar(df_1B_3,'y','Y1_start (ns)',savepath)
show_channel_bar(df_1B_3,'y','Y1_end (ns)',savepath)
show_channel_bar(df_1B_3,'y','Y2_start (ns)',savepath)
show_channel_bar(df_1B_3,'y','Y2_end (ns)',savepath)


df_2A = wrap_df("./L2_a.csv")
df_2A_1 = df_slices_time(df_2A,'2024-10-01 00:00:00','2024-10-15 00:00:00')
df_2A_2 = df_slices_time(df_2A,'2024-10-15 00:00:00','2024-11-01 00:00:00')
df_2A_3 = df_slices_time(df_2A,'2024-11-01 00:00:00','2024-11-15 00:00:00')

df_2B = wrap_df("./L2_b.csv")
df_2B_1 = df_slices_time(df_2B,'2024-10-01 00:00:00','2024-10-15 00:00:00')
df_2B_2 = df_slices_time(df_2B,'2024-10-15 00:00:00','2024-11-01 00:00:00')
df_2B_3 = df_slices_time(df_2B,'2024-11-01 00:00:00','2024-11-15 00:00:00')
