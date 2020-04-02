import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D
import xlwt
import Planetary_data_file as pd
import tools as t

hours=3600.0
days=hours*24

def null_perts():
    return {
        'J2':False,
        'aerodrag':False,
        'thrust':0,
        'thrust_direction':0,
        'isp':0,
        'Cd':0,
        'rho':0,
        'A':0
    }


def norm(v):
    return np.linalg.norm(v)

def normed(v):
    return np.array(v)/norm(v)

class orbit_propagator:
    
    def __init__(self,initial_state,time_span,time_step,coes=False,deg=True,mass0=0,perts=null_perts(),cb=pd.earth,propagator='lsoda',sc={}):
        
        if coes:
            self.r0,self.v0,_=t.coes2rv(initial_state,deg=deg,mu=cb['mu'])
        else:
            self.r0 = initial_state[:3]
            self.v0 = initial_state[3:]
        self.cb=cb
        self.time_step=time_step
        self.mass0=mass0                                                
        self.time_span=time_span

 
    
        self.n_steps = int(np.ceil(self.time_span/self.time_step))+1                                                   # ceil. function rounds float up to nearest whole number and int. transforms the float to a interger
        self.ts=np.zeros((self.n_steps+1,1))                                                                                 # initialise arrays
        self.y=np.zeros((self.n_steps+1,7))
        self.alts=np.zeros((self.n_steps+1))
        self.propagator=propagator                                                                      #6 states (vx,vy,vz,ax,ay,az) preallocating memory (instead of creating a new list it allows memory to overwrite existing list
        self.step = 0
        

        self.y[0,:] = self.r0.tolist() + self.v0.tolist()+[self.mass0]                    #initial condition at first step
        self.alts[0]=t.norm(self.r0)-self.cb['radius']


        self.solver = ode(self.ODE)                                                  # initiate solver (lsoda)fast, high order
        self.solver.set_integrator(self.propagator)                                     # Adam-Bashford multistep
        self.solver.set_initial_value(self.y[0,:],0)                                        # initial state

        self.perts=perts
    

        #store stop conditions and dictionary
        self.stop_conditions_dict=sc

        #define dictionary to map internals method
        self.stop_conditions_map={'min_alt':self.check_min_alt}
        
        #create stop conditions function list with deorbit always checked
        self.stop_condition_functions=[self.check_deorbit]

        #fill in the rest of the stop conditions.
        for key in self.stop_conditions_dict.keys():
            if key in self.stop_conditions_map:
                self.stop_condition_functions.append(self.stop_conditions_map[key])

            
        #propagate the orbit
        self.propagate_orbit()


    #check if satellite has deorbited
    def check_deorbit(self):
        if self.alts[self.step]<self.cb['deorbit_altitude']:
            print('Satellite deorbited after %.1f seconds' % self.ts[self.step])

            return False
        return True


    #check if minimum altitude exceeded
    def check_min_alt(self):
        if self.alts[self.step]<self.stop_conditions_dict['min_alt']:
            self.ts_in_hours=self.ts[self.step]/3600.0
            self.ts_in_days=self.ts[self.step]/86400.0
            
            print('Satellite reached minimum altitude after %.1f seconds' % self.ts[self.step])
            print('Which amounts to %.1f hours' % self.ts_in_hours)
            print('Which is approximately %.1f days' % self.ts_in_days)
            
            return False
        return True

    #function called at each timestep to check all stop conditions
    def check_stop_conditions(self):
            
        #for each stop condition
        for sc in self.stop_condition_functions:

            #if returns False
            if not sc():

                #stop conditions reached and will return False
                return False

        #if no stop conditions reached, return true.
        return True

        
    #propagate orbit on the initial conditions defined within the init() function
    def propagate_orbit(self):
        print('Propagating orbit...')

        
        ##propagate orbit. check for max time and stop conditions at each time step 
        while self.solver.successful() and self.step<self.n_steps and self.check_stop_conditions():
            # propogate orbit integrator step
            self.solver.integrate(self.solver.t+self.time_step)                                                            # while its successful the solver can have a number of errors
            self.step += 1

            
            self.ts[self.step] = self.solver.t                                                                       # i.e time step can be to small or too rigid
            self.y[self.step] = self.solver.y

            self.alts[self.step]=t.norm(self.solver.y[:3])-self.cb['radius']
          

        self.ts=self.ts[:self.step]
        self.rs = self.y[:self.step,:3]                                                                    # extract the position array(60x6) we want all rows and all steps up to upto coloum 0,1,2
        self.vs = self.y[:self.step,3:6]
        self.masses=self.y[:self.step,6]
        self.alts=self.alts[:self.step]

    def ODE(self,t_,y,):                                                           
        rx,ry,rz,vx,vy,vz,mass = y                                                             
        r = np.array([rx,ry,rz])
        v = np.array([vx,vy,vz])
        #print(norm(v))

                                                                              # norm of the radius vector because because perbubations require the norm of the input - this lowers computational cost
        norm_r = np.linalg.norm(r)
                                                                              # linalg is a sub library of numpy for equations and methods
        a = -r * self.cb['mu'] / norm_r**3                                    # law of gravitation, as r is vector a has output as a vector


        ##orbit propagator J2 ######################################################################################################
        if self.perts['J2']:
            
            z2=r[2]**2
            r2=norm_r**2
            
            tx=r[0]/norm_r*(5*z2/r2-1)
            ty=r[1]/norm_r*(5*z2/r2-1)
            tz=r[2]/norm_r*(5*z2/r2-3)
            
            a+=1.5*self.cb['J2']*self.cb['mu']*self.cb['radius']**2.0/norm_r**4.0*np.array([tx,ty,tz])

        #calculate aerodynamic drag
        if self.perts['aerodrag']:

            #calculate altitude and air density
            z=norm_r-self.cb['radius']  # find altitude
            rho=t.calc_atmospheric_density(z) #find air density at given altitude

            #calculate motion of s/c with repsect to a rotating atmosphere
            v_rel=v-np.cross(self.cb['atm_rot_vector'],r)

            drag=-v_rel*0.5*rho*t.norm(v_rel)*self.perts['Cd']*self.perts['A']/mass

            a+=drag

        #calculate thrus
        if self.perts['thrust']:
            a+=(self.perts['thrust_direction']*t.normed(v)*self.perts['thrust']/mass)/1000.0
            mass_flow=-self.perts['thrust']/(self.perts['isp']*9.81)

        return [vx,vy,vz,a[0],a[1],a[2], mass_flow]

    def calculate_coes(self,degrees=True,print_results=False):
        print('Calculating COEs....')

        self.coes=np.zeros((self.n_steps,6))

        for n in range(self.n_steps):
            self.coes[n,:]=t.rv2coes(self.rs[n,:],self.vs[n,:],mu=self.cb['mu'],degrees=degrees)


    def plot_coes(self,hours=False,days=False,show_plot=False,save_plot=False,title='Change In COES at 35786km with BPT4000',figsize=(16,8)):
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

    def plot_alts(self,show_plot=False,save_plot=False,hours=False,days=False,title='Satellite Mass v.s Manoeuvre Time - Geostationary Orbit, XIPS-25 (35786km)',figsize=(16,8),dpi=500):


        if hours:
            ts=self.ts/3600.0
            xunit='Time Elapsed (hours'
            
        elif self.days:
            ts=self.days/(3600.0/24.0)
            xunit='Time Elapsed (days'

        else:
            ts=self.ts
            xunit='Time Elapsed (seconds)'

        plt.figure(figsize=figsize)

        SIZE = 20
        SIZE1 = 12

        
        plt.plot(ts,self.alts, Label="XIPS-25")

        plt.rc('font', size=SIZE1)          # controls default text sizes
        plt.rc('axes', titlesize=SIZE1)     # fontsize of the axes title
        plt.rc('axes', labelsize=SIZE1)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SIZE1)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SIZE1)    # fontsize of the tick labels
        plt.rc('legend', fontsize=SIZE1)    # legend fontsize
        plt.rc('figure', titlesize=SIZE1)  # fontsize of the figure title          
        plt.grid(True)
        plt.xlabel('%s)'% xunit)
        plt.ylabel('Altitude (km)')
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
            plt.savefig(title+'.png',dpi=500)


    


