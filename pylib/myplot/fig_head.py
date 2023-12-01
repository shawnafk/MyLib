import numpy as np
import matplotlib.pyplot as plt
import os
plt.style.use('/Users/shawn/Nextcloud/MyLib/pylib/myplot/prl.mplstyle')
shadeopts = {'cmap': 'jet', 'shading': 'gouraud',"rasterized":True}
#for p in os.environ.get('PYTHONPATH').split(':'):
#	fname=p+'./prl.mplstyle'
#	if os.path.exists(fname):
#		plt.style.use(fname)
#		break
import matplotlib as mpl
import matplotlib.font_manager as font_manager
font = font_manager.FontProperties(family='Times New Roman',
                                   weight='bold',
                                   style='normal')
import matplotlib.ticker as ticker
from matplotlib.pyplot import MultipleLocator
axlbfs=10
def multiple_formatter(denominator=2, number=np.pi, latex='\pi'):
    def gcd(a, b):
        while b:
            a, b = b, a%b
        return a
    def _multiple_formatter(x, pos):
        den = denominator
        num = int(np.rint(den*x/number))
        com = gcd(num,den)
        (num,den) = (int(num/com),int(den/com))
        if den==1:
            if num==0:
                return r'$0$'
            if num==1:
                return r'$%s$'%latex
            elif num==-1:
                return r'$-%s$'%latex
            else:
                return r'$%s%s$'%(num,latex)
        else:
            if num==1:
                return r'$\frac{%s}{%s}$'%(latex,den)
            elif num==-1:
                return r'$\frac{-%s}{%s}$'%(latex,den)
            else:
                return r'$\frac{%s%s}{%s}$'%(num,latex,den)
    return _multiple_formatter
class Multiple:
    def __init__(self, denominator=2, number=np.pi, latex='\pi'):
        self.denominator = denominator
        self.number = number
        self.latex = latex
    def locator(self):
        return plt.MultipleLocator(self.number / self.denominator)
    def formatter(self):
        return plt.FuncFormatter(multiple_formatter(self.denominator, self.number, self.latex))
def add_top_cax(ax,pad,width,ratio):
    axpos=ax.get_position()
    #x0 y0 x1 y1
    caxpos=mpl.transforms.Bbox.from_extents(
        axpos.x0 + (1-ratio)/2,
        axpos.y1 + pad,
        axpos.x1 - (1-ratio)/2,
        axpos.y1 + pad + width,
        ) 
    cax = ax.figure.add_axes(caxpos)
    return cax
def add_right_cax(ax,pad,width):
    axpos=ax.get_position()
    caxpos=mpl.transforms.Bbox.from_extents(
        axpos.x1+pad,
        axpos.y0,
        axpos.x1+pad+width,
        axpos.y1,
        )
    cax = ax.figure.add_axes(caxpos)
    return cax
'''
x0=ax1.get_position().get_points()[0][0]
x1=ax1.get_position().get_points()[1][0]
y0=ax1.get_position().get_points()[0][1]
y1=ax1.get_position().get_points()[1][1]
pos=fig.add_axes([x0,y1+5/100.,x1-x0,0.05])
plt.colorbar(im,cax=pos,orientation='horizontal')
'''

from mpl_toolkits.mplot3d import axes3d
class MyAxes3D(axes3d.Axes3D):

    def __init__(self, baseObject, sides_to_draw):
        self.__class__ = type(baseObject.__class__.__name__,
                              (self.__class__, baseObject.__class__),
                              {})
        self.__dict__ = baseObject.__dict__
        self.sides_to_draw = list(sides_to_draw)
        self.mouse_init()

    def set_some_features_visibility(self, visible):
        for t in self.w_zaxis.get_ticklines() + self.w_zaxis.get_ticklabels():
            t.set_visible(visible)
        self.w_zaxis.line.set_visible(visible)
        self.w_zaxis.pane.set_visible(visible)
        self.w_zaxis.label.set_visible(visible)

    def draw(self, renderer):
        # set visibility of some features False 
        self.set_some_features_visibility(False)
        # draw the axes
        super(MyAxes3D, self).draw(renderer)
        # set visibility of some features True. 
        # This could be adapted to set your features to desired visibility, 
        # e.g. storing the previous values and restoring the values
        self.set_some_features_visibility(True)

        zaxis = self.zaxis
        draw_grid_old = zaxis.axes._draw_grid
        # disable draw grid
        zaxis.axes._draw_grid = False

        tmp_planes = zaxis._PLANES

        if 'l' in self.sides_to_draw :
            # draw zaxis on the left side
            zaxis._PLANES = (tmp_planes[2], tmp_planes[3],
                             tmp_planes[0], tmp_planes[1],
                             tmp_planes[4], tmp_planes[5])
            zaxis.draw(renderer)
        if 'r' in self.sides_to_draw :
            # draw zaxis on the right side
            zaxis._PLANES = (tmp_planes[3], tmp_planes[2], 
                             tmp_planes[1], tmp_planes[0], 
                             tmp_planes[4], tmp_planes[5])
            zaxis.draw(renderer)

        zaxis._PLANES = tmp_planes

        # disable draw grid
        zaxis.axes._draw_grid = draw_grid_old
