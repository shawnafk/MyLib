import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def draw_rectangle(ax,length, width):
    # 计算矩形左下角的坐标
    left_bottom_x = -length / 2
    left_bottom_y = -width / 2

    rect = Rectangle((left_bottom_x, left_bottom_y), length, width, linewidth=1, edgecolor='r', facecolor='none')
    ax.add_patch(rect)
    return 0
