##########################################################################################################################################################################################
#       TITLE & NOTES

#   Gary Stewart, ME409,420,421 Individual Project
#   Traditional two body problem will form the basis of simulation work for project
#   numpy used as a toolbox for scientific computing and scipy as the integrator (ode solver)
#   matplotlib and mpl are standard tools for plotting graphs in 2D and 3D in matplotlib

#   all supporting documentation on functions can be found on add-in packages websites
#   equations and constants @ 'orbital mechanics for engineering students, Howard Curtis'

###########################################################################################################################################################################################
##      IMPORTS

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D


##########################################################################################################################################################################################
#       ORBITAL PLOTTING FUNCTION

def plot(r):                                   # create 3D plot
    fig = plt.figure(figsize=(50,50))          # projection - '3d' essential import
    ax = fig.add_subplot(111,projection='3d')  # add subplot 111 - 1st row,1st column 1st value


                                                                        # plor trajectory and starting point
    ax.plot(r[:,0],r[:,1],r[:,2],'k', label='Trajectory')               # satallite trajectory plot
    ax.plot([r[0,0]],[r[0,1]],[r[0,2]],'ko', label ='Initial Position') # satellites initial position plot

    r_plot = earth_radius   # define earths radius


                                                 # plot earth
    _u,_v = np.mgrid [0:2*np.pi:20j,0:np.pi:20j] # define sphere  (VIDEO L2 EXPLANATION)
    _x = r_plot * np.cos(_u) * np.sin(_v)        # trig
    _y = r_plot * np.sin(_u) * np.sin(_v)        # trig
    _z = r_plot * np.cos(_v)                     # trig
    ax.plot_surface(_x,_y,_z, cmap='Blues')      # surface plot (x,y,z variables cmap=colour plot)


                                         # plot the x, y, z vectors
    l=r_plot*2.0
    x,y,z = [[0,0,0],[0,0,0],[0,0,0]]    # origin of arrow plot
    u,v,w = [[50,0,0],[0,50,0],[0,0,50]] # finish of arrow plot

    ax.quiver(x,y,z,u,v,w,color='k')  # quiver is the arrow function with the above arguements and k=colour
    max_val=np.max(np.abs(r))         # this helps normalise the axis and displays equal magnitudes i.e cubic looking


                                      # set labels and titles
    ax.set_xlim([-max_val,max_val])
    ax.set_ylim([-max_val,max_val])
    ax.set_zlim([-max_val,max_val])

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    ax.set_title('Satellite orbit in GSO') # title

    plt.show() # show plot

#########################################################################################################################################################################################
#       IMPLIMENT ODE SOLVER


earth_radius = 6378.0        # units: km (equatorial)
earth_mu = 10000.0           # units: km3/s2 (gravitational constant for earth - see newtons law of gravitation)


def diffy_q(t,y,mu):              # first imput into the differential equation solver
    rx,ry,rz,vx,vy,vz = y         # unpack state: the ode is a function solver and needs time, state and mu
    r = np.array([rx,ry,rz])      # distance/positional array to be a vector to be used in the law of gravitation


                                          # norm of the radius vector because because perbubations require the norm of the input - this lowers computational cost
    norm_r = np.linalg.norm(r)            # linalg is a sub library of numpy for equations and methods
    ax,ay,az = -r * mu / norm_r**3
                                   # law of gravitation, as r is vector a has output as a vector
    ay = ay + 0.000000015
    ax = ax + 0.000000015
    return [vx,vy,vz,ax,ay,az]            # input = state(position, velocity) so we want to return derivative(velocity,accelleration)


if __name__ == '__main__':                # special variable which defines the code is being written in main script and not imported


                                            # initial conditions of orbit parameters
    r_mag = earth_radius + 10000.0          # magnitutde of orbit distance(km)
    v_mag = np.sqrt(earth_mu / r_mag)       # magnitude of velocity of a circular orbit follows this equation



    r0 = [r_mag,0,0]        # initial position and velcity vectors                  )
    v0 = [0,v_mag,0]        # velocity will always be perpendicular to position (centrifugal force



    tspan = 100 * 24 * 60.0 * 60.0        # timespan of simulation,  we know duration of one orbit in GSO = 24 hours
    dt = 100.0                         # in order to get our total number of steps (timespam/timestep)
    n_steps = int(np.ceil(tspan/dt))  # ceil. function rounds float up to nearest whole number and int. transforms the float to a interger


                             # initialise arrays
    ys=np.zeros((n_steps,6)) # (6 states (vx,vy,vz,ax,ay,az) preallocating memory (instead of creating a new list it allows memory to overwrite existing list
    ts=np.zeros((n_steps,1)) # (1 state (time)



                                 #initial conditions
    x0 = r0 + v0                 #add lists together to concatenate (not element by element)
    ys[0] = np.array(x0)         #initial condition at first step
    step = 1



    solver = ode(diffy_q)                   # initiate solver (lsoda)fast, high order
    solver.set_integrator('lsoda')          # Adam-Bashford multistep
    solver.set_initial_value(x0,0)          # initial state
    solver.set_f_params(earth_mu)           # define 3rd argument mu

        




    while solver.successful() and step<n_steps:     # propogate orbit, solver does its work, timestep to small
        solver.integrate(solver.t+dt)               # while its successful the solver can have a number of errors
        ts[step] = solver.t                         # i.e time step can be to small or too rigid
        ys[step] = solver.y                         # step<n_step means that after time steps done we exit while loop
        step += 1


    rs = ys[:,:3]  # extract the position array(60x6) we want all rows and all steps up to upto coloum 0,1,2

    plot(rs)

 #       END CODE
################################################################################################################################################################################################
