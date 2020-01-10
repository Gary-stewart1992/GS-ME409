##Main Script
import numpy as np
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D
import planetary_data_file as pd
import tools as t
from orbit_propagator import orbit_propagator as op
from orbit_propagator import null_perts


tspan = 1000 * 3600                  
dt = 10.0

cb=pd.earth


if __name__ == '__main__':                          # special variable which defines the code is being written in main script and not imported
    
    perts=null_perts()
    perts['thrust']=0.327
    perts['thrust_direction']=-1
    perts['isp']=4300
    perts['J2']=True

    mass0 = 1400 #kg

    #Galileo-022
    op0=op(t.tle2coes('galileo.txt'),tspan,dt,coes=True,deg=True,perts=perts)

    #SKYNET-4C com sat (1500kg,2.2kW) 
    #op1=op(t.tle2coes('SKYNET4Ctle.txt'),tspan,dt,coes=True,deg=True,perts=perts)

    #GLONASS 
    #op2=op(t.tle2coes('GLONASS-M.txt'),tspan,dt,coes=True,deg=True,perts=perts)

    #POLAR
    #op3=op(t.tle2coes('POLAR.txt'),tspan,dt,coes=True,deg=True,perts=perts)

    #POLAR
    op4=op(t.tle2coes('ISS.txt'),tspan,dt,coes=True,deg=True,perts=perts)

    #t.plot_n_orbits([op0.rs,op1.rs,op2.rs,op3.rs,op4.rs],labels=['GALIL022','SKYNET4C','GLONASS-M','POLAR','ISS'],title=['Multiple Orbits'], show_plot=True)

    op4.plot_3d(show_plot=True)






