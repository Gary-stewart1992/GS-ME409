##Main Script
import numpy as np
import datetime
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D
import planetary_data_file as pd
import tools as t
from orbit_propagator import orbit_propagator as op
from orbit_propagator import null_perts


tspan = 48* 3600                   
dt = 100.0

cb=pd.earth


if __name__ == '__main__':                          # special variable which defines the code is being written in main script and not imported
    
    perts=null_perts()
    perts['aero']=False
    perts['Cd']=2.2
    perts['A']=(1e-3)**2/4.0 #km^2
    perts['thrust']=0.327
    perts['thrust_direction']=-1
    perts['isp']=4300
    perts['J2']=True
 
    Mass0= 150.0


  
   

    #Galileo-022
    op0=op(t.tle2coes('galileo.txt'),tspan,dt,coes=True,deg=True,perts=perts)

    #SKYNET-4C com sat (1500kg,2.2kW) 
    op1=op(t.tle2coes('SKYNET4Ctle.txt'),tspan,dt,coes=True,deg=True,perts=perts)

    #GLONASS 
    op2=op(t.tle2coes('GLONASS-M.txt'),tspan,dt,coes=True,deg=True,perts=perts)

    #POLAR
    op3=op(t.tle2coes('POLAR.txt'),tspan,dt,coes=True,deg=True,perts=perts)

    #ISS
    op4=op(t.tle2coes('ISS.txt'),tspan,dt,coes=True,deg=True,perts=perts)

    state0=t.tle2coes('ISS.txt',deg=True)
    



    #book example propagate
    op=op(state0,tspan,dt,perts=perts,deg=True,coes=True)
    op.plot_3d(show_plot=True)
    op.calculate_coes()
    op.plot_coes(show_plot=True,hours=True)

     

















  #book example
    #r0=np.array([-2384.46,5729.01,3050.45])
    #v0=np.array([-7.36138,-2.98997,1.64352])
    #state0=np.array(t.rv2coes(r0,v0,print_results=True))
    #t.plot_n_orbits([op0.rs,op1.rs,op2.rs,op3.rs,op4.rs],labels=['GALIL022','SKYNET4C','GLONASS-M','POLAR','ISS'],title=['Multiple Orbits'], show_plot=True)


     #ISS propogate
    #op4.plot_3d(show_plot=True)
    #op4.calculate_coes()
    #op4.plot_coes(show_plot=True,hours=True)




