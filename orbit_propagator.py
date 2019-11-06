import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

import Planetary_data_file as pd

class orbit_propagator:

    def __init__(self,r0,v0,tspan,dt,cb=pd.earth):
        self.r0=r0
        self.v0=v0
        self.tspan=tspan
        self.dt=dt
        self.cb=cb

    def propagate_orbit(self):

            self.n_steps = int(np.ceil(self.tspan/self.dt))  # ceil. function rounds float up to nearest whole number and int. transforms the float to a interger


                                                 # initialise arrays
            self.ys=np.zeros((self.n_steps,6)) # (6 states (vx,vy,vz,ax,ay,az) preallocating memory (instead of creating a new list it allows memory to overwrite existing list
            self.ts=np.zeros((self.n_steps,1)) # (1 state (time)



                                                        #initial conditions
            self.y0 = self.r0 + self.v0                 #add lists together to concatenate (not element by element)
            self.ys[0] = np.array(self.y0)              #initial condition at first step
            self.step = 1



            self.solver = ode(self.diffy_q)                   # initiate solver (lsoda)fast, high order
            self.solver.set_integrator('lsoda')               # Adam-Bashford multistep
            self.solver.set_initial_value(self.y0,0)          # initial state






            while self.solver.successful() and self.step<self.n_steps:     # propogate orbit, solver does its work, timestep to small
                self.solver.integrate(self.solver.t+self.dt)                             # while its successful the solver can have a number of errors
                self.ts[self.step] = self.solver.t                                               # i.e time step can be to small or too rigid
                self.ys[self.step] = self.solver.y                         # step<n_step means that after time steps done we exit while loop
                self.step += 1


                self.rs = self.ys[:,:3]  # extract the position array(60x6) we want all rows and all steps up to upto coloum 0,1,2
                self.vs = self.ys[:,3:]


    def diffy_q(self,t,y,):              # first imput into the differential equation solver
        rx,ry,rz,vx,vy,vz = y         # unpack state: the ode is a function solver and needs time, state and mu
        r = np.array([rx,ry,rz])      # distance/positional array to be a vector to be used in the law of gravitation


                                                          # norm of the radius vector because because perbubations require the norm of the input - this lowers computational cost
        norm_r = np.linalg.norm(r)                       # linalg is a sub library of numpy for equations and methods
        ax,ay,az = -r * self.cb['mu'] / norm_r**3        # law of gravitation, as r is vector a has output as a vector
        return [vx,vy,vz,ax,ay,az]                      #input = state(position, velocity) so we want to return derivative(velocity,accelleration)


    def plot_3d(self,show_plot=False,save_plot=False):
        
        fig = plt.figure(figsize=(50,50))          # projection - '3d' essential import
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

        ax.set_title('Satellite orbit in GSO') # title

        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=300)


