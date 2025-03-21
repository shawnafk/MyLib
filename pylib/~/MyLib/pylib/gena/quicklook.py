# ... existing code ...
def cartoon2d_start(xyextent,datetimes,array,*args,**kwargs):
    fig=plt.figure()
    ax=plt.axes()
    # Bug fix: Adjust indentation
    #aux.draw_rectangle(ax,150,60)
    #aux.draw_rectangle(ax,90,30)
    time = ax.annotate(0,xy=(0.2, 0.9),xycoords='figure fraction')
    im = ax.imshow(array[0,...],aspect='auto',extent=xyextent,*args,**kwargs)
    cbar = fig.colorbar(im)
    def animate(i):
        Z=array[i,:,:]
        im.set_array(Z)
        im_min, im_max = np.min(Z), np.max(Z)
        im.set_clim(im_min, im_max)
        time.set_text(datetimes[i])
    anim = animation.FuncAnimation(fig, animate, *args,**kwargs)
    return anim
# ... existing code ...