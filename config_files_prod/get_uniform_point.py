from random import *
import numpy as np
#import matplotlib.pyplot as plt


def get_tp_point(theta_min = 0, theta_max = 1.1, phi_min = 0, phi_max = np.pi/2):
    ''' This is a function that gives back points uniformly spread accross the dome
    This is done by ensuring that the number of points on a given ring is proportionnal 
    to the perimeter of that ring (hence the linear function in theta in the if condition)
    '''
    a = False
    while a != True:
        u = random()
        v = random()
        test = random()
        theta =  2 * np.pi * u
        phi = 2 * np.pi * v
        b = 0
        c = 0.4
#        theta_max = 1.1
        if (test) <= c * u + b  and theta <= theta_max  and phi <= phi_max and phi >= phi_min and theta >= theta_min:#
            a = True
            return theta, phi

'''
Notes from the meeting:
piece-wise polynomial response function
fake data points fitting to the data - moving them around
response function
spline of fake data points  - fitting to the data
ask Mark whether it is a problem that two of my points are actual points? - does give weight but it's the end points and also whether we are closer or not to others might also influence this so yeah...
wait the next steps  - sort out the stats and binning first before worrying about the number of hits
then pred vs data in the fit:
response function prediction (map in colors - with bins substraction to compare properly ! and geographically)
can also predict maps completely -> pull the map fully out of the simulated splines -> good
next: expand in the R direction my fits.
'''


