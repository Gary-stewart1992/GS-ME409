##Gary Stewart
##Main Script

import numpy as np
import datetime
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D
import planetary_data_file as pd
import tools as t
from orbit_propagator import orbit_propagator as op
from orbit_propagator import null_perts



#span of time in which the simulation occurs
tspan = 350*24*3600
#time step
dt = 500.0

#loads earths data from the planatery data file
cb=pd.earth

# special variable which defines the code is being written in main script and not imported
if __name__ == '__main__':   


    #attributes values to the null dictionary defined in orbit propagator 
    perts=null_perts()
    perts['thrust']=0.28
    perts['thrust_direction']=-1
    perts['isp']=3550
    perts['J2']=True
    perts['aerodrag']=True
    perts['Cd']=2.2
    perts['A']= 5.2e-6 #km2

    #satellite dry mass 
    mass0 = 2000

    #initialises the solver
    op=op(t.tle2coes('EUTELSAT3B.txt'),tspan,dt,coes=True,deg=False,mass0=mass0,perts=perts)

    #plots altitude vs time plot 
    op.plot_alts(show_plot=True,hours=True)
    
    #plots 3D plot 
    op.plot_3d(show_plot=True)

    #calculates plots
    op.calculate_coes()

    #plots coes, requires time within bounds of r1  and r2
    op.plot_coes(show_plot=True,hours=True)


 




