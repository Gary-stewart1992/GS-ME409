import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math as m
import datetime

import Planetary_data_file as pd

d2r = np.pi/180.0

def plot_n_orbits(rs,labels,cb=pd.earth, show_plot=False,save_plot=False, title='Multiple Orbits'):
                    
    fig = plt.figure(figsize=(16,8))          # projection - '3d' essential import
    ax = fig.add_subplot(111,projection='3d')  # add subplot 111 - 1st row,1st column 1st value

    n=0
    for r in rs:     
        ax.plot(r[:,0],r[:,1],r[:,2],label=labels[n])               # satallite trajectory plot
        ax.plot([r[0,0]],[r[0,1]],[r[0,2]],'ko')                    # satellites initial position plot
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

    ax.set_title('Multiple Orbits') # title
    plt.legend()

    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(title+'.png',dpi=300)


def coes2rv(coes,deg=False,mu=pd.earth['mu']):
    if deg:
        a,e,i,ta,aop,raan,date=coes
        i*=d2r
        ta*=d2r
        aop*=d2r
        raan*=d2r
    else:
        a,e,i,ta,aop,raan,date = coes

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

    return r,v,date

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
            

def tle2coes(tle_filename,mu=pd.earth['mu']):
    
    #read the tle text file from celetrak
    with open(tle_filename, 'r') as f:
        lines=f.readlines()

    #break the text file into three seperate lines
    line0 = lines[0].strip() #satellite name
    line1 = lines[1].strip().split
    line2 = lines[2].strip().split

    #epoch (year and day)
    epoch=line1()[3]
    year,month,day,hour=calc_epoch(epoch)

    #gather COE's from .txt

    #inclination
    i=float(line2()[2])*d2r #in radians

    #RAAN
    raan=float(line2()[3])*d2r #in radians

    #eccentricity
    e_string=line2()[4]
    e=float('0.'+e_string)

    #arguement of perigee
    aop=float(line2()[5])*d2r #in radians

    #mean_anomaly
    Me=float(line2()[6])*d2r #in radians

    #mean motion
    mean_motion=float(line2()[7]) #in revolutions per day

    #period
    T=1/mean_motion*24*3600 #in units of seconds

    #semi-major axis
    a=(T**2*mu/4.0/np.pi**2)**(1/3.0)

    #calculate eccentric anomaly
    E=ecc_anomaly([Me,e],'newton')

    #calculate true anomally
    ta=true_anomaly([E,e])

    #magnitude of radius vector
    r_mag=a*(1-e*np.cos(E))    #####potential issue

    return a,e,i,ta,aop,raan,[year,month,day,hour]

def calc_epoch(epoch):

    #epoch year
    year=int('20'+epoch[:2])

    epoch=epoch[2:].split('.')

    #day of year
    day_of_year=int(epoch[0])-1

    #decimal hour of day
    hour=float('0.'+epoch[1])*24

    #year/month/day
    date=datetime.date(year,1,1)+datetime.timedelta(day_of_year)

    #extract month and day from tle
    month=float(date.month)
    day=float(date.day)

    return year,month,day,hour

def true_anomaly(arr):
    E,e=arr
    return 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2.0))

def tle2rv(tle_filename):
    return coes2rv(tle2coes(tle_filename))



    
    
    

