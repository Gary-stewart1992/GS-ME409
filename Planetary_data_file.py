#Gary Stewart
#Planatary data file

#Imports
import math as m
import datetime
import numpy as np


#Universal constant
G_meters = 6.67408e-11
G=G_meters*10**-9 #required units


#sun dictionary definition
sun={
    'name':'Sun',
    'mass':1.989e30,
    'mu':1.32712e11,
    'radius':695700.0

}

#earth atmospheric array
atm=np.array([[63.096,2.059e-4],[251.189,5.909e-11],[1000.0,3.561e-15]]) #[altitude,density]

#earth value dictionary
earth={
    
        'name':'Earth',
        'mass':5.972e24,        
        'mu':5.972e24*G,
        'radius':6378.0,      #earth radius
        'J2':-1.082635854e-3, #j2 constant
        'deorbit_altitude':0, #deorbit altitude (km)
        'zs':atm[:,0], #altitude (km)
        'rhos':atm[:,1]*10**8, #density (kgm3)
        'atm_rot_vector':np.array([0.0,0.0,72.9211e-6]) #how earth rotates (rads/s)

}
