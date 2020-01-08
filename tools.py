import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math as m 

import planetary_data_file as pd

d2r = np.pi/180.0

def plot_n_orbits(rs,labels,cb=pd.earth, show_plot=False,save_plot=False, title='Multiple Orbits'):
                    
    fig = plt.figure(figsize=(16,8))          # projection - '3d' essential import
    ax = fig.add_subplot(111,projection='3d')  # add subplot 111 - 1st row,1st column 1st value

    n=0
    for r in rs:     
        ax.plot(r[:,0],r[:,1],r[:,2], label=labels[n])               # satallite trajectory plot
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

    ax.quiver(x,y,z,u,v,w,color='k')                            # quiver is the arrow function with the above arguements and k=colour
    max_val=np.max(np.abs(rs))          # this helps normalise the axis and displays equal magnitudes i.e cubic looking


                                                  # set labels and titles
    ax.set_xlim([-max_val,max_val])
    ax.set_ylim([-max_val,max_val])
    ax.set_zlim([-max_val,max_val])

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    ax.set_title('Title') # title
    plt.legend()

    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(title+'.png',dpi=300)


def coes2rv(coes,deg=False,mu=pd.earth['mu']):
    if deg:
        a,e,i,ta,aop,raan=coes
        i*=d2r
        ta*=d2r
        aop*=d2r
        raan*=d2r
    else:
        a,e,i,ta,aop,raan = coes

    E=ecc_anomaly([ta,e], 'tae')

    r_norm=a*(1-e**2)/(1+e*np.cos(ta))

    #calculate r and v vectors in perifocal frame
    r_perif = r_norm*np.array([m.cos(ta),m.sin(ta),0])
    v_perif=m.sqrt(mu*a)/r_norm*np.array([-m.sin(E),m.cos(E)*m.sqrt(1-e**2),0])

    #rotation matrix from the perifocal to ECI
    perif2eci = np.transpose(eci2perif(raan,aop,i))

    #calculate r and v vectors in inertial frame
    r = np.dot(perif2eci,r_perif)
    v = np.dot(perif2eci,v_perif)

    return r,v

def eci2perif(raan,aop,i):
    
    row0 =[-m.sin(raan)*m.cos(i)*m.sin(aop)+ m.cos(raan)*m.cos(aop), m.cos(raan)*m.cos(i)*m.sin(aop)+ m.sin(raan)*m.cos(aop), m.sin(i)*m.sin(aop)]
    row1 =[-m.sin(raan)*m.cos(i)*m.cos(aop)- m.cos(raan)*m.sin(aop), m.cos(raan)*m.cos(i)*m.cos(aop)- m.sin(raan)*m.sin(aop), m.sin(i)*m.cos(aop)]
    row2 =[m.sin(raan)*m.sin(i), -m.cos(raan)*m.sin(i), m.cos(i)]

    return np.array([row0,row1,row2])

def ecc_anomaly(arr,method,tol=1e-8):
    if method=='newton':

        #newtons method for iteratively finding E
        Me,e=arr
        if Me<np.pi/2.0: E0=Me+e/2.0
        else: E0 = Me- e
        for n in range(200):
            ratio=(E0-e*np.sin(E0)-Me)/(1-e*np.cos(E0));
            if abs(ratio)<tol:
                if n==0: return E0
                else: return E1
            else:
                E1=E0-ratio
                E0=E1

        #failure to converge
        return False
    elif method == 'tae':
        ta,e=arr
        return 2*m.atan(m.sqrt((1-e)/(1+e))*m.tan(ta/2.0))
    else:
        print('Invalid method for eccentric anomaly')
            

    
    

