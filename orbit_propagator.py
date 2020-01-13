import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D

import Planetary_data_file as pd
import tools as t

hours=3600.0
days=hours*24

def null_perts():
    return {
        'aero':False,
        'Cd':0,
        'A':0,
        'mu':0,
        'rho':False,
        'J2':False,
        'Aerodrag':False,
        'thrust':0,
        'thrust_direction':0,
        'isp':0
    }


def norm(v):
    return np.linalg.norm(v)

def normed(v):
    return np.array(v)/norm(v)

class orbit_propagator:
    
    def __init__(self,state0,tspan,dt,coes=False,deg=True,mass0=0,perts=null_perts(),cb=pd.earth,propagator='lsoda'):
        
        if coes:
            self.r0,self.v0,_=t.coes2rv(state0,deg=deg,mu=cb['mu'])
        else:
            self.r0 = state0[:3]
            self.v0 = state0[3:]
        self.cb=cb
        self.dt=dt
        self.mass0=mass0                                                
        self.tspan=tspan

 
    
        self.n_steps = int(np.ceil(self.tspan/self.dt))+1                                                   # ceil. function rounds float up to nearest whole number and int. transforms the float to a interger
        self.ts=np.zeros((self.n_steps,1))                                                                                 # initialise arrays
        self.y=np.zeros((self.n_steps,7))
        self.propagator=propagator                                                                      #6 states (vx,vy,vz,ax,ay,az) preallocating memory (instead of creating a new list it allows memory to overwrite existing list
        self.step = 1
        

        self.y[0,:] = self.r0.tolist() + self.v0.tolist()+[self.mass0]                    #initial condition at first step

        self.solver = ode(self.diffy_q)                                                  # initiate solver (lsoda)fast, high order
        self.solver.set_integrator(self.propagator)                                     # Adam-Bashford multistep
        self.solver.set_initial_value(self.y[0,:],0)                                        # initial state

        self.perts=perts
        
        self.propagate_orbit()

    def propagate_orbit(self):

        print('Propagating orbit...')

        while self.solver.successful() and self.step<self.n_steps:                   # propogate orbit, solver does its work, timestep to small
            self.solver.integrate(self.solver.t+self.dt)                             # while its successful the solver can have a number of errors

            self.ts[self.step] = self.solver.t                                               # i.e time step can be to small or too rigid
            self.y[self.step] = self.solver.y                                        # step<n_step means that after time steps done we exit while loop
            self.step += 1


        self.ts=self.ts[:self.step]
        self.rs = self.y[:self.step,:3]                                                      # extract the position array(60x6) we want all rows and all steps up to upto coloum 0,1,2
        self.vs = self.y[:self.step,3:6]
        self.masses=self.y[:self.step,-1]
        self.alts=(np.linalg.norm(self.rs,axis=1)-self.cb['radius']).reshape((self.step,1))


    def diffy_q(self,t_,y,):                                                           # first imput into the differential equation solver
        rx,ry,rz,vx,vy,vz,mass= y                                                             # unpack state: the ode is a function solver and needs time, state and mu
        r = np.array([rx,ry,rz])
        v = np.array([vx,vy,vz])# distance/positional array to be a vector to be used in the law of gravitation
        

                                                                              # norm of the radius vector because because perbubations require the norm of the input - this lowers computational cost
        norm_r = np.linalg.norm(r)
                                                                              # linalg is a sub library of numpy for equations and methods
        a = -r * self.cb['mu'] / norm_r**3                                    # law of gravitation, as r is vector a has output as a vector

        ##orbit propagator J2 
        if self.perts['J2']:
            z2=r[2]**2
            r2=norm_r**2
            tx=r[0]/norm_r*(5*z2/r2-1)
            ty=r[1]/norm_r*(5*z2/r2-1)
            tz=r[2]/norm_r*(5*z2/r2-3)
            a+=1.5*self.cb['J2']*self.cb['mu']*self.cb['radius']**2.0/norm_r**4.0*np.array([tx,ty,tz])


#        #aero drag calculations
#        if self.perts['aero']:
#            z=norm_r-self.cb['radius']
#            rho=t.calc_atmospheric_density(z)
#
#
#            v_rel=v-np.cross(self.cb['atm_rot_vector'],r)
#            drag=-v_rel*0.5*rhp*t.norm(v_rel)*self.perts['Cd']*self.perts['A']/self.mass
#            a+=drag

        if self.perts['thrust']:
            a+self.perts['thrust_direction']*normed(v)*self.perts['thrust']/mass/1000.0
            dmdt=-self.perts['thrust']/self.perts['isp']/9.81


        return [vx,vy,vz,a[0],a[1],a[2], dmdt]


    def calculate_coes(self,degrees=True,print_results=False):
        print('Calculating COEs....')

        self.coes=np.zeros((self.n_steps,6))

        for n in range(self.n_steps):
            self.coes[n,:]=t.rv2coes(self.rs[n,:],self.vs[n,:],mu=self.cb['mu'],degrees=degrees)



    def plot_coes(self,hours=False,days=False,show_plot=False,save_plot=False,title='Change in Orbital Elements',figsize=(16,8)):
        print('Plotting COEs...')

        fig1,axs1 =plt.subplots(nrows=2,ncols=3,figsize=figsize)

        #figure titles 
        fig1.suptitle(title,fontsize=20)

        #x-axis
        if hours:
            ts=self.ts/3600.0
            xlabel='Time (hours)'

        elif self.days:
            ts=self.days/3600.0/24.0
            xlabel='Time Elapsed (days)'

        else:
            ts=self.ts
            xlabel='Time Elapsed (seconds)'


        fig1.tight_layout(pad=6.0)
        
        
        #plotting true anomaly
        axs1[0,0].plot(ts,self.coes[:,3])
        axs1[0,0].set_title('True Anomaly vs. Time')
        axs1[0,0].grid(True)
        axs1[0,0].set_ylabel('True Anomaly (deg)')
        axs1[1,1].set_xlabel(xlabel)

        #plotting semi major axis
        axs1[1,0].plot(ts,self.coes[:,0])
        axs1[1,0].set_title('Semi-Major Axis vs. Time')
        axs1[1,0].grid(True)
        axs1[1,0].set_ylabel('Semi-Major Axis (km)')
        axs1[1,0].set_xlabel(xlabel)

        #plotting eccentricity
        axs1[0,1].plot(ts,self.coes[:,1])
        axs1[0,1].set_title('Eccentricity vs. Time')
        axs1[0,1].grid(True)
        axs1[0,1].set_xlabel(xlabel)

        #plotting argument of periapse
        axs1[0,2].plot(ts,self.coes[:,4])
        axs1[0,2].set_title('Argument of Perigee vs. Time')
        axs1[0,2].grid(True)
        axs1[0,2].set_ylabel('Argument of Perigee (deg)')
        axs1[0,2].set_xlabel(xlabel)
        
        #plotting inclination
        axs1[1,1].plot(ts,self.coes[:,2])
        axs1[1,1].set_title('Inclination vs. Time')
        axs1[1,1].grid(True)
        axs1[1,1].set_ylabel('Inclination (deg)')
        axs1[1,1].set_xlabel(xlabel)

         #plotting raan
        axs1[1,2].plot(ts,self.coes[:,5])
        axs1[1,2].set_title('RAAN vs. Time')
        axs1[1,2].grid(True)
        axs1[1,2].set_ylabel('RAAN (deg)')
        axs1[1,2].set_xlabel(xlabel)

        if show_plot:
            plt.show()

        if save_plot:
            plt.savefig(title+'.png',dpi=500)


    def plot_alts(self,show_plot=False,save_plot=False,hours=False,days=False,title='Radial Distance vs. Time',figsize=(16,8),dpi=500):

        if hours:
            ts=self.ts/3600.0
            xunit='Time (hours)'

        elif self.days:
            ts=self.days/(3600.0/24.0)
            xunit='Time Elapsed (days)'

        else:
            ts=self.ts
            xunit='Time Elapsed (seconds)'

        plt.figure(figsize=figsize)
        plt.plot(ts,self.alts, 'k')
        plt.grid(True)
        plt.xlabel('Time (%s)' % xunit)
        plt.ylabel('altitude (km)')
        plt.title(title)
        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=dpi)

            
    


    def plot_3d(self,show_plot=False,save_plot=False, title='Deorbiting Manoeuvre Trajectory',dpi=500):
        
        fig0 = plt.figure(figsize=(16,8))          # projection - '3d' essential import
        ax0 = fig0.add_subplot(111,projection='3d')  # add subplot 111 - 1st row,1st column 1st value


                                                                        # plor trajectory and starting point
        ax0.plot(self.rs[:,0],self.rs[:,1],self.rs[:,2],'k', label='Trajectory')               # satallite trajectory plot
        ax0.plot([self.rs[0,0]],[self.rs[0,1]],[self.rs[0,2]],'ko', label ='Initial Position') # satellites initial position plot


                                                 # plot earth
        _u,_v = np.mgrid [0:2*np.pi:20j,0:np.pi:20j] # define sphere  (VIDEO L2 EXPLANATION)
        _x = self.cb['radius'] * np.cos(_u) * np.sin(_v)        # trig
        _y = self.cb['radius'] * np.sin(_u) * np.sin(_v)        # trig
        _z = self.cb['radius'] * np.cos(_v)                     # trig
        ax0.plot_surface(_x,_y,_z, cmap='Blues')      # surface plot (x,y,z variables cmap=colour plot)


                                         # plot the x, y, z vectors
        l=self.cb['radius']*2.0
        x,y,z = [[0,0,0],[0,0,0],[0,0,0]]    # origin of arrow plot
        u,v,w = [[50,0,0],[0,50,0],[0,0,50]] # finish of arrow plot

        ax0.quiver(x,y,z,u,v,w,color='k')  # quiver is the arrow function with the above arguements and k=colour
        max_val=np.max(np.abs(self.rs))         # this helps normalise the axis and displays equal magnitudes i.e cubic looking



                                      # set labels and titles
        ax0.set_xlim([-max_val,max_val])
        ax0.set_ylim([-max_val,max_val])
        ax0.set_zlim([-max_val,max_val])

        ax0.set_xlabel('X (km)')
        ax0.set_ylabel('Y (km)')
        ax0.set_zlabel('Z (km)')

        ax0.set_title('Electric Propulsion Manoeuvre Trajectory') # title

        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=300)


  





