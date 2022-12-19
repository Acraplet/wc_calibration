from random import *
import numpy as np
import matplotlib.pyplot as plt


def get_tp_point(theta_max = 1.1, phi_min = 0, phi_max = np.pi/2):
    a = False
    while a != True:
        u = random()
        v = random()
        test = random()

        theta =  2 * np.pi * u
        phi = 2 * np.pi * v

        b = 0
        c = 0.4

        theta_max = 1.1
        if (test) <= c * u + b  and theta <= theta_max  and phi <= phi_max and phi >= phi_min:#
            a = True
            return theta, phi
        #ax2 = fig.add_subplot(122, projection='polar')
        #ax2.set_title("Q Comparision \n %s-%s"%(filename, comparision_file))

        #zr2 = (df['Q']/df['events'] - df_compa['Q']/df_compa['events']) #/df['Q']/df['events']
        #im2 = ax2.scatter(df['phi'], df['theta'], c = zr2, cmap = "nipy_spectral",s = 40)
        #col3 = plt.colorbar(im2, label=f"Difference in number of *raw* P.E. in mPMT58 per photon", orientation="vertical", format= "%.2f")

        #ax2.set_thetamin(0)
        #ax2.set_thetamax(phi_max)
        #ax2.xaxis.labelpad = 10

        #ax2.set_xlabel(f'$\Theta$(rad)')
        #plt.savefig('/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_pictures/Comparisions/Non-interpolated/%s.png'%outputfile_name)
        ##plt.show()
    #plt.plot(hist.T)
    #plt.show()


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


