import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def draw_rectangle(ax,length, width):
    # 计算矩形左下角的坐标
    left_bottom_x = -length / 2
    left_bottom_y = -width / 2

    rect = Rectangle((left_bottom_x, left_bottom_y), length, width, linewidth=1, edgecolor='r', facecolor='none')
    ax.add_patch(rect)
    return 0

import pickle
def openfig(fname):
    with open(fname, 'rb') as f:
        loaded_fig = pickle.load(f)
    return loaded_fig

def savefig(fname,fig):
    with open(name, 'wb') as f:
        pickle.dump(fig, f)
    return 0
