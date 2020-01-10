import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D

import Planetary_data_file as pd
import tools as t

def null_perts():
    return {
        'aero':False,
        'Cd':0,
        'A':0,
        'mu':0,
        'rho':False,
        'J2':False,
        'Aerodrag':False,
        'thrust':False,
        'thrust_direction':False,
        'isp':0
    }

class orbit_propagator:
    
    def __init__(self,state0,tspan,dt,coes=False,deg=True,cb=pd.earth,perts=null_perts()):
        
        if coes:
            self.r0,self.v0,_=t.coes2rv(state0,deg=deg,mu=cb['mu'])
        else:
            self.r0 = state0[:3]
            self.v0 = state0[3:]
            
                                                     #POTENTIAL ERROR  self.y0 = self.r0.tolist()+self.v0.tolist()..>> self.r0=r0 and self.v0=v0
        self.y0=self.r0.tolist()+self.v0.tolist()  
        self.tspan=tspan
        self.dt=dt
        self.cb=cb

        self.n_steps = int(np.ceil(self.tspan/self.dt)) +1     # ceil. function rounds float up to nearest whole number and int. transforms the float to a interger


        self.ts=np.zeros((self.n_steps,1))                                                                                 # initialise arrays
        self.ys=np.zeros((self.n_steps,6))                   # (6 states (vx,vy,vz,ax,ay,az) preallocating memory (instead of creating a new list it allows memory to overwrite existing list
        self.ts[0]=0
        self.ys[0,:] = self.y0                    #initial condition at first step
        self.step = 1 



        self.solver = ode(self.diffy_q)                   # initiate solver (lsoda)fast, high order
        self.solver.set_integrator('lsoda')               # Adam-Bashford multistep
        self.solver.set_initial_value(self.y0,0)          # initial state

        self.perts=perts
        
        self.propagate_orbit()

    def propagate_orbit(self):

            while self.solver.successful() and self.step<self.n_steps:                   # propogate orbit, solver does its work, timestep to small
                self.solver.integrate(self.solver.t+self.dt)                             # while its successful the solver can have a number of errors
                self.ts[self.step] = self.solver.t                                               # i.e time step can be to small or too rigid
                self.ys[self.step] = self.solver.y                                        # step<n_step means that after time steps done we exit while loop
                self.step += 1


                self.rs = self.ys[:,:3]                                                      # extract the position array(60x6) we want all rows and all steps up to upto coloum 0,1,2
                self.vs = self.ys[:,3:]


    def diffy_q(self,t,y,):                                                           # first imput into the differential equation solver
        rx,ry,rz,vx,vy,vz = y                                                             # unpack state: the ode is a function solver and needs time, state and mu
        r = np.array([rx,ry,rz])
        v = np.array([vx,vy,vz])                                                    # distance/positional array to be a vector to be used in the law of gravitation


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

        #aero drag calculations
        if self.perts['aero']:

            #calculate altitude and air density
            z=norm_r-self.cb['radius']
            rho=t.calc_atmospheric_density(z)

            #calculate motion of s/c wrt rotating frame

            v_rel=v-np.cross(self.cb['atm_rot_vector'],r)
            
            drag=-v_rel*0.5*rhp*t.norm(v_rel)*self.perts['Cd']*self.perts['A']/self.mass
            
            a+=drag

            

        return [vx,vy,vz,a[0],a[1],a[2]]


    def calculate_coes(self,degrees=True):
        print('Calculating COEs....')

        self.coes=np.zeros((self.n_steps,6))

        for n in range(self.n_steps):
            self.coes[n,:]=t.rv2coes(self.rs[n,:],self.vs[n,:],mu=self.cb['mu'],degrees=degrees)
        

    def plot_coes(self,hours=False,days=False,show_plot=False,save_plot=False,title='Change in Classical Orbital Elements',figsize=(16,8)):
        print('Plotting COEs...')

        fig,axs =plt.subplots(nrows=2,ncols=3,figsize=figsize)

        #figure titles 
        fig.suptitle(title,fontsize=20)

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


        fig.tight_layout(pad=6.0)
        
        
        #plotting true anomaly
        axs[0,0].plot(ts,self.coes[:,3])
        axs[0,0].set_title('True Anomaly vs. Time')
        axs[0,0].grid(True)
        axs[0,0].set_ylabel('True Anomaly (deg)')
        axs[1,1].set_xlabel(xlabel)

        #plotting semi major axis
        axs[1,0].plot(ts,self.coes[:,0])
        axs[1,0].set_title('Semi-Major Axis vs. Time')
        axs[1,0].grid(True)
        axs[1,0].set_ylabel('Semi-Major Axis (km)')
        axs[1,0].set_xlabel(xlabel)

        #plotting eccentricity
        axs[0,1].plot(ts,self.coes[:,1])
        axs[0,1].set_title('Eccentricity vs. Time')
        axs[0,1].grid(True)
        axs[0,1].set_xlabel(xlabel)

        #plotting argument of periapse
        axs[0,2].plot(ts,self.coes[:,4])
        axs[0,2].set_title('Argument of Perigee vs. Time')
        axs[0,2].grid(True)
        axs[0,2].set_ylabel('Argument of Perigee (deg)')
        axs[0,2].set_xlabel(xlabel)
        
        #plotting inclination
        axs[1,1].plot(ts,self.coes[:,2])
        axs[1,1].set_title('Inclination vs. Time')
        axs[1,1].grid(True)
        axs[1,1].set_ylabel('Inclination (deg)')
        axs[1,1].set_xlabel(xlabel)

         #plotting raan
        axs[1,2].plot(ts,self.coes[:,5])
        axs[1,2].set_title('RAAN vs. Time')
        axs[1,2].grid(True)
        axs[1,2].set_ylabel('RAAN (deg)')
        axs[1,2].set_xlabel(xlabel)

        if show_plot:
            plt.show()

        if save_plot:
            plt.savefig(title+'.png',dpi=300)


    def plot_3d(self,show_plot=False,save_plot=False):
        
        fig = plt.figure(figsize=(16,8))          # projection - '3d' essential import
        ax = fig.add_subplot(111,projection='3d')  # add subplot 111 - 1st row,1st column 1st value


                                                                        # plor trajectory and starting point
        ax.plot(self.rs[:,0],self.rs[:,1],self.rs[:,2],'k', label='Trajectory')               # satallite trajectory plot
        ax.plot([self.rs[0,0]],[self.rs[0,1]],[self.rs[0,2]],'ko', label ='Initial Position') # satellites initial position plot


                                                 # plot earth
        _u,_v = np.mgrid [0:2*np.pi:20j,0:np.pi:20j] # define sphere  (VIDEO L2 EXPLANATION)
        _x = self.cb['radius'] * np.cos(_u) * np.sin(_v)        # trig
        _y = self.cb['radius'] * np.sin(_u) * np.sin(_v)        # trig
        _z = self.cb['radius'] * np.cos(_v)                     # trig
        ax.plot_surface(_x,_y,_z, cmap='Blues')      # surface plot (x,y,z variables cmap=colour plot)


                                         # plot the x, y, z vectors
        l=self.cb['radius']*2.0
        x,y,z = [[0,0,0],[0,0,0],[0,0,0]]    # origin of arrow plot
        u,v,w = [[50,0,0],[0,50,0],[0,0,50]] # finish of arrow plot

        ax.quiver(x,y,z,u,v,w,color='k')  # quiver is the arrow function with the above arguements and k=colour
        max_val=np.max(np.abs(self.rs))         # this helps normalise the axis and displays equal magnitudes i.e cubic looking


                                      # set labels and titles
        ax.set_xlim([-max_val,max_val])
        ax.set_ylim([-max_val,max_val])
        ax.set_zlim([-max_val,max_val])

        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('Z (km)')

        ax.set_title('Electric Propulsion Manoeuver Trajectory') # title

        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=300)



