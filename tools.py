import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import Planetary_data_file as pd

d2r=np.pi/180

def plot_n_orbits(rs,labels,cb=pd.earth,show_plot=False,save_plot=False):
        
    fig = plt.figure(figsize=(50,50))          # projection - '3d' essential import
    ax = fig.add_subplot(111,projection='3d')  # add subplot 111 - 1st row,1st column 1st value

    n = 0
    for r in rs:                                                                # plor trajectory and starting point
        ax.plot(r[:,0],r[:,1],r[:,2],'k', label=labels[n])               # satallite trajectory plot
        ax.plot([r[0,0]],[r[0,1]],[r[0,2]]) # satellites initial position plot
        n+=1

                                             # plot earth
    _u,_v = np.mgrid [0:2*np.pi:20j,0:np.pi:20j] # define sphere  (VIDEO L2 EXPLANATION)
    _x = cb['radius'] * np.cos(_u) * np.sin(_v)        # trig
    _y = cb['radius'] * np.sin(_u) * np.sin(_v)        # trig
    _z = cb['radius'] * np.cos(_v)                     # trig
    ax.plot_surface(_x,_y,_z, cmap='Blues')      # surface plot (x,y,z variables cmap=colour plot)


                                     # plot the x, y, z vectors
    l=cb['radius']*2.0
    x,y,z = [[0,0,0],[0,0,0],[0,0,0]]    # origin of arrow plot
    u,v,w = [[50,0,0],[0,50,0],[0,0,50]] # finish of arrow plot

    ax.quiver(x,y,z,u,v,w,color='k')  # quiver is the arrow function with the above arguements and k=colour
    max_val=np.max(np.abs(rs))         # this helps normalise the axis and displays equal magnitudes i.e cubic looking


                                  # set labels and titles
    ax.set_xlim([-max_val,max_val])
    ax.set_ylim([-max_val,max_val])
    ax.set_zlim([-max_val,max_val])

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    ax.set_title('Satellite orbit in GSO') # title

    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(title+'.png',dpi=300)
