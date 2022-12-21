from random import *
import numpy as np
import matplotlib.pyplot as plt

c = 0
d = 0
e = 0

polar = False
n = 100000000
if polar == True:
    n = 100000
theta_list = []
phi_list = []
test_list = []
phi_min = 0
phi_max = np.pi/2
if polar == True:
    ax = plt.axes(projection = 'polar')
for i in range(n):
    u = random()
    v = random()
    test = random()
    theta =  2 * np.pi * u
    phi = 2 * np.pi * v
    a = 28.3699314
    b = 0
    c = 0.4
    frac = 0.1
    theta_max = 1.1
    if (test) <= c * u + b  and theta <= theta_max  and phi <= phi_max and phi >= phi_min:
        if polar == True:
            ax.plot(phi, theta, 'x')
        if 0.1 * (1 - frac) <= theta and 0.1 * (1 + frac)>= theta:
            c = c +1
        if 0.6 * (1 - frac) <= theta and 0.6 * (1 + frac)>= theta:
            e = e +1
        if 1 * (1 - frac) <= theta and 1 * (1 + frac)>= theta:
            d = d +1

        theta_list.append(theta)
        phi_list.append(phi)
        test_list.append(test)

if polar == False:
    theta_list = np.array(theta_list)
    plt.figure()
    plt.hist(theta_list, bins = 50)
    plt.xlabel(f'$\Theta$ of sampled points')
    plt.ylabel('Occurences')

    fig=plt.figure()
    ax=plt.axes()
    pc = ax.hist2d(phi_list, np.sin(theta_list), bins = 20, cmap = "viridis")
    cb = fig.colorbar(pc[3])
    cb.ax.set_label("counts in bin")
    ax.set_xlabel(f'$\phi$ (rad)')
    ax.set_ylabel(f'sin($\Theta$)')
    plt.show()



print('\n', c/(2*np.pi*0.1*len(theta_list)) * 100, e/(2*np.pi * 0.6*len(theta_list)) * 100, d/(2 * np.pi *1.0 *len(theta_list)) * 100, '\n')

if polar == True:
    # define binning
    rbins = np.linspace(0, 2*np.pi, 10)
    abins = np.linspace(0, 1.1, 10)

    #calculate histogram
    hist, _, _ = np.histogram2d(phi_list, theta_list, bins=(rbins, abins))
    R, A= np.meshgrid(rbins, abins)

    # plot
    fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))

    pc = ax.pcolormesh(R, A, hist.T, cmap="magma_r")
    fig.colorbar(pc)
    ax.grid(True)
    plt.show()

    plt.plot(test_list, theta_list, 'x')
    plt.xlabel('test')
    plt.ylabel('theta')
    plt.show()

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


