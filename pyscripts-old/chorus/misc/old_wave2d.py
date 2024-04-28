import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

dT=48.8988646329616
zpos=np.loadtxt("../data/ori/zpos.out")
sour  = np.load("../data/nj1_beta0_adv_sour.npy")
tstep = np.arange(0,sour.shape[0])*dT
Z,T = np.meshgrid(zpos,tstep)

fig,ax=plt.subplots()
ax.pcolormesh(Z,T,sour[:,:-1,1],cmap='bwr',shading='gouraud',rasterized=True,vmin=-1e-7,vmax=1e-7)
#axins = inset_axes(ax, width="40%", height="30%", loc='lower left', bbox_to_anchor=(0.1, 0.1, 1, 1), bbox_transform=ax.transAxes) 
axins = ax.inset_axes((0.15, 0.2, 0.4, 0.3))
axins.pcolormesh(Z,T,sour[:,:-1,1],cmap='bwr',shading='gouraud',rasterized=True,vmin=-1e-7,vmax=1e-7)
Zvl=[-1.0e3,-1.98e3]
Tvl=[3.73e4,4.18e4]
vleft = - abs( (Zvl[0]-Zvl[1])/(Tvl[0]-Tvl[1]))
axins.plot(Zvl,Tvl,'k--',lw=2.5)

Zvr=[-1.4e3,-9.7e2]
Tvr=[4.05e4,4.278e4]
vright =  abs( (Zvr[0]-Zvr[1])/(Tvr[0]-Tvr[1]))
axins.plot(Zvr,Tvr,'k-',lw=2.5)
# 设置放大区间
zone_left = -2000
zone_right = -900
zone_top = 4.4e4
zone_bottom = 3.6e4

#axins.set_xlim(zone_left, zone_right)
#axins.set_ylim(zone_bottom, zone_top)
# 坐标轴的扩展比例（根据实际数据调整）
x_ratio = 0. # x轴显示范围的扩展比例
y_ratio = 0. # y轴显示范围的扩展比例

# X轴的显示范围
zone_width  = zone_right - zone_left
xlim0 = zone_left-(zone_width)*x_ratio
xlim1 = zone_right+(zone_width)*x_ratio

# Y轴的显示范围
zone_height  = zone_top - zone_bottom
ylim0 = zone_bottom-(zone_height)*y_ratio
ylim1 = zone_top+(zone_height)*y_ratio

# 调整子坐标系的显示范围
axins.set_xlim(xlim0, xlim1)
axins.set_ylim(ylim0, ylim1)
  
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec='k', lw=1)
