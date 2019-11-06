import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D

from orbit_propagator import orbit_propagator as op
import planetary_data_file as pd

cb=pd.earth


if __name__ == '__main__': # special variable which defines the code is being written in main script and not imported

                                            # initial conditions of orbit parameters
    r_mag = cb['radius'] + 35786.0          # magnitutde of orbit distance(km)
    v_mag = np.sqrt(cb['mu'] / r_mag)       # magnitude of velocity of a circular orbit follows this equation




    r0 = [r_mag,r_mag*0.1,r_mag*-0.1]        # initial position and velcity vectors                  )
    v0 = [0,v_mag,v_mag*0.3]                          # velocity will always be perpendicular to position (centrifugal force



    tspan = 24.0 * 3600              # timespan of simulation,  we know duration of one orbit in GSO = 24 hours
    dt = 20.0                         # in order to get our total number of steps (timespam/timestep)
    
    op=op(r0,v0,tspan,dt)
    op.propagate_orbit()
    op.plot_3d(show_plot=True)
