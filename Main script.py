import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D

import Planetary_data_file as pd
import n_plot_tool as npt
from orbit_propagator import orbit_propagator as op

tspan = 48 * 3600
dt = 20.0

cb=pd.earth


if __name__ == '__main__': # special variable which defines the code is being written in main script and not imported

                                            # initial conditions of orbit parameters
    r_mag = cb['radius'] + 36000         # magnitutde of orbit distance(km)
    v_mag = np.sqrt(cb['mu'] / r_mag)       # magnitude of velocity of a circular orbit follows this equation
    r0=np.array([r_mag,0,0])                        
    v0=np.array([0,v_mag,0])
    
    r_mag = cb['radius'] + 16000
    v_mag = np.sqrt(cb['mu'] / r_mag)
    r00 =np.array([r_mag,0,0])                        
    v00 = np.array([0,v_mag,0.9]) 

    
    op0=op(r0,v0,tspan,dt)
    op00=op(r00,v00,tspan,dt)
    
    op0.propagate_orbit()
    op00.propagate_orbit()
    
    npt.plot_n_orbits([op0.rs,op00.rs],labels=['LEO','GEO'],show_plot=True)
