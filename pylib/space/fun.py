import numpy as np
import matplotlib.pyplot  as plt
def plot_earth(ax):
    t = np.linspace(0,2*np.pi,128)
    x = np.cos(t)
    y = np.sin(t)
    ax.plot(x,y,'k')
    return 0
