##Main Script
import numpy as np
import datetime
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D
import planetary_data_file as pd
import tools as t
from orbit_propagator import orbit_propagator as op
from orbit_propagator import null_perts


tspan = 45*24*3600                   
dt = 500.0


cb=pd.earth


if __name__ == '__main__':                          # special variable which defines the code is being written in main script and not imported
    
    perts=null_perts()
    perts['thrust']=0.27
    perts['thrust_direction']=-1
    perts['isp']=1950
    perts['J2']=True
    perts['aerodrag']=True
    perts['Cd']=2.2
    perts['A']=(1e-3)**2/4.0

    mass0 = 700
    
    
    op=op(t.tle2coes('IRIDIUM.txt'),tspan,dt,coes=True,deg=False,mass0=mass0,perts=perts)

    op.plot_alts(show_plot=True,hours=True)
    op.plot_3d(show_plot=True)
    op.calculate_coes()
    op.plot_coes(show_plot=True,hours=True)
     


#Satellite TLE's 
    #EUTELSATK
    #Galileo
    #03BFM5
    #SKYNET-4Ctle
    #GLONASS-M
    #IRIDIUM

 




