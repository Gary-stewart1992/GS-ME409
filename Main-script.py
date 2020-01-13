##Main Script
import numpy as np
import datetime
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D
import Planetary_data_file as pd
import tools as t
from orbit_propagator import orbit_propagator as op
from orbit_propagator import null_perts


tspan = 24*3600                   
dt = 250.0

cb=pd.earth


if __name__ == '__main__':                          # special variable which defines the code is being written in main script and not imported
    
    perts=null_perts()
    perts['aero']=False
    perts['Cd']=2.2
    perts['A']=(1e-3)**2/4.0 #km^2
    perts['thrust']=0.5
    perts['thrust_direction']=1
    perts['isp']=5000.0
    perts['J2']=True



    mass0 = 1500.0 #kg

    #SKYNET-4C com sat (1500kg,2.2kW) 
    op=op(t.tle2coes('ISS.txt'),tspan,dt,coes=True,deg=True,mass0=mass0,perts=perts)


    
    op.plot_alts(show_plot=True,hours=True)
    op.plot_3d(show_plot=True)             
    op.calculate_coes()
    op.plot_coes(show_plot=True,hours=True)
     


















    #Galileo-022
    #op0=op(t.tle2coes('galileo.txt'),tspan,dt,coes=True,deg=True,perts=perts)

    #SKYNET-4C com sat (1500kg,2.2kW) 
    #op1=op(t.tle2coes('SKYNET4Ctle.txt'),tspan,dt,coes=True,deg=True,perts=perts)

    #GLONASS 
    #op2=op(t.tle2coes('GLONASS-M.txt'),tspan,dt,coes=True,deg=True,perts=perts)

    #POLAR
    #op3=op(t.tle2coes('POLAR.txt'),tspan,dt,coes=True,deg=True,perts=perts)

    #ISS
    #op=op(t.tle2coes('ISS.txt'),tspan,dt,coes=True,deg=True,mass0=mass0,perts=perts)




