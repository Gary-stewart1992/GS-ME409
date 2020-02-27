##Main Script
import numpy as np
import datetime
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D
import Planetary_data_file as pd
import tools as t
from orbit_propagator import orbit_propagator as op
from orbit_propagator import null_perts


time_span = (24*3600)*38.8                 
time_step = 50.0

cb=pd.earth


if __name__ == '__main__':                          # special variable which defines the code is being written in main script and not imported
    perts['thrust']=0.05
    perts['thrust_direction']=-1
    perts['isp']=310
    perts['J2']=True

    
    #defined stop conditions library
    sc={'min_alt':350,
        'max_alt':43000.0
    }

    mass0 = 100.0 #kg

    #SKYNET-4C com sat (1500kg,2.2kW) 
    op=op(t.tle2coes('03BFM5.txt'),time_span,time_step,coes=True,deg=False,mass0=mass0,perts=perts,sc=sc)

    op.plot_alts(show_plot=True,hours=True)
    op.plot_3d(show_plot=True)             
    op.calculate_coes()
    op.plot_coes(show_plot=True,hours=True)
     

##Satellite TLE's
#03BFM5.txt
#galileo.txt
#SKYNET-4C.txt
#GLONASS-M 
#POLAR.txt
#ISS





