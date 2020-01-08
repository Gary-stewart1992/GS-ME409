##Main Script
import numpy as np
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D

import tools as t
from orbit_propagator import orbit_propagator as op
import planetary_data_file as pd


tspan = 24.0 * 3600                  
dt = 20.0

cb=pd.earth


if __name__ == '__main__':                          # special variable which defines the code is being written in main script and not imported

    c0=[cb['radius']+4000.0,0.0006189,51.6393,0.0,234.1955,105.6372]  #ISS

    c1=[cb['radius']+35800.0,0.0,0.0,0.0,0.0,0.0]                    #GEO

    c2=[cb['radius']+6000.0,0.3,20.0,0.0,15.0,40.0]                  #MEO
        
        
    op0=op(c0,tspan,dt,coes=True)
    op1=op(c1,tspan,dt,coes=True)
    op2=op(c2,tspan,dt,coes=True)

    op0.propagate_orbit()
    op1.propagate_orbit()
    op2.propagate_orbit()

    t.plot_n_orbits([op0.rs,op1.rs,op2.rs],labels=['ISS','GSO','MEO'],title=['Multiple Orbits'], show_plot=True)



