#Gary Stewart
#Tool functions

#add in packages
import math as m
import datetime
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#user defined
import Planetary_data_file as pd

#radians/degree conversion
d2r = np.pi/180.0
r2d = 180.0/np.pi

#normalised velocity and unit velocity vector 
def norm(v):
    return np.linalg.norm(v)

def normed(v):
    return np.array(v)/norm(v)


#same as plot_3d in orbit_propagator but was used to plot multiple orbits on one plot.
def plot_n_orbits(rs,label,cb=pd.earth, show_plot=False,save_plot=False, title='Multiple Orbits',dpi=500):
    
    #projection - '3d' essential import                
    fig = plt.figure(figsize=(16,8))

    #add subplot 111 - 1st row,1st column 1st value
    ax = fig.add_subplot(111,projection='3d')  

    #counter
    n=0
    for r in rs:
        
        # satallite trajectory plot
        ax.plot(r[:,0],r[:,1],r[:,2],label='Trajectory')
        
        # satellites initial position plot
        ax.plot([r[0,0]],[r[0,1]],[r[0,2]], 'ko')                    
        n+=1

    #lot earth
    _u,_v = np.mgrid [0:2*np.pi:20j,0:np.pi:20j]       # define sphere
    _x = cb['radius'] * np.cos(_u) * np.sin(_v)        # trig
    _y = cb['radius'] * np.sin(_u) * np.sin(_v)        # trig
    _z = cb['radius'] * np.cos(_v)                     # trig
    ax.plot_surface(_x,_y,_z, cmap='Blues')            # surface plot (x,y,z variables cmap=colour plot)


    #plot the x, y, z vectors
    l=cb['radius']*2.0
    x,y,z = [[0,0,0],[0,0,0],[0,0,0]]               # origin of arrow plot
    u,v,w = [[500,0,0],[0,500,0],[0,0,500], 'ko']   # finish of arrow plot

    #quiver is the arrow function with the above arguements and k=colour
    ax.quiver(x,y,z,u,v,w,color='k')

    #this helps normalise the axis and displays equal magnitudes i.e cubic looking
    max_val=np.max(np.abs(rs))          


    #set labels and titles
    ax.set_xlim([-max_val,max_val])
    ax.set_ylim([-max_val,max_val])
    ax.set_zlim([-max_val,max_val])

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    ax.set_title('Deorbiting Manoeuvre') # title
    plt.legend()

    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(title+'.png',dpi=500)


#COES2RV algorithm at Rene-Schwarz m003
def coes2rv(coes,deg=True,mu=pd.earth['mu']):

    #if imput is in degrees
    if deg:
        a,e,i,ta,aop,raan,date=coes
        i*=d2r
        ta*=d2r
        aop*=d2r
        raan*=d2r

    #if input is radians    
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

#RV2COES Algorithm at Rene-Schwarz m002
def rv2coes(r,v,mu=pd.earth['mu'],degrees=False,print_results=False):
    
    #norm of position vector
    r_norm=norm(r)


    #specific angular momentum vector
    h=np.cross(r,v)
    h_norm=norm(h)

    #inclination
    i=m.acos(h[2]/h_norm)

    #eccentricity vector
    e=((norm(v)**2-mu/r_norm)*r-np.dot(r,v)*v)/mu

    #eccentricity scalar
    e_norm=norm(e)

    #node-line
    N=np.cross([0,0,1],h)
    N_norm=norm(N)

    #RAAN
    raan=m.acos(N[0]/N_norm)
    if N[1]<0: raan=2*np.pi-raan #quad check

    #argument of perigee
    aop=m.acos(np.dot(N,e)/N_norm/e_norm)
    if e[2]<0: aop=2*np.pi-aop #another quad check

    #true anomaly
    ta=m.acos(np.dot(e,r)/e_norm/r_norm)
    if np.dot(r,v)<0: ta=2*np.pi-ta #another quad check

    #semi-major axis
    a=r_norm*(1+e_norm*m.cos(ta))/(1-e_norm**2)

    if print_results:
        print('a',a)
        print('e',e_norm)
        print('i',i*r2d)
        print('RAAN',raan*r2d)
        print('AOP',aop*r2d)
        print('TA',ta*r2d)


    #convert to degrees if it has been specificed 
    if degrees: return [a,e_norm,i*r2d,ta*r2d,aop*r2d,raan*r2d]
    else: return [a,e_norm,i,ta,aop,raan]

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
            
#used to convert two line elements to orbital elements
def tle2coes(tle_filename,mu=pd.earth['mu'],deg=True,print_results=False):
    
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


#calculate atmospheric density from given altitude
def calc_atmospheric_density(z):
    rhos,zs=find_rho_z(z)
    if rhos[0]==0: return 0.0

    Hi=-(zs[1]-zs[0])/np.log(rhos[1]/rhos[0])

    return rhos[0]*np.exp(-(z-zs[0])/Hi)

#find endpoints of altitude and density surrounding input altitude
def find_rho_z(z,zs=pd.earth['zs'],rhos=pd.earth['rhos']):
    if not 1.0<z<1000.0:
        return [[0.0,0.0],[0.0,0.0]]

    #find the two poiunts surrounding the given input
    for n in range(len(rhos)-1):
        if zs[n]<z<zs[n+1]:
            return[[rhos[n],rhos[n+1]],[zs[n],zs[n+1]]]

    #if out of range return zeros
    return [[0.0,0.0],[0.0,0.0]]
                

    
    
    

