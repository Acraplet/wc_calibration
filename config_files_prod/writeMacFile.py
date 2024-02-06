#This is a python file to write multiple occurences of a WCSim config file stored in WCSim_config folder

#It runs with the following inputs:
# (R) R - distance from the PMT to the source
# (t) nTheta - the number of theta points we want uniformly spaced in sin theta smaller than theta_max
# (p) nPhi - the number of phi points we want uniformly spaced in phi smaller than phi_max
# (a) absff - the abwff absorption coefficient in WCSim
# (r) rayff - the rayleigh scattering coefficient in WCSim
# (d) random_positions - to sample uniformly spaced in phi and sin(theta) but not nicely arranged as it is without this option toggled: here we input the number of points that we want to have   
# (u) uniform_positions - to sample uniformly i.e. a given number of points per area on the sphere: input the number of points we want to sample  -> This is what we need to use

#to use: python writeMacFile.py -R 10 -t 10 -p 10 -f 10 -a 0.000555 -r 10e10  -u/d 3810 

#the total number of points is therefore nTheta*nPhi
#Note that the source positions are defined with respect to the centre of mPMT 58 (in WCTE) the dimensions of the origin of the light is hard coded in Module_writeFiles.py 

#import uproot
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import getopt
import get_uniform_point as gp
import random
import Module_writeFiles as md
#Here are the default values for R, nTheta, nPhi
nTheta = 10
nPhi = 10
R = 10
nEvent = 1000 #default number of events: 1000

#Works naively  with neatly arranged (nTheta, nPhi) positions 
ordered_positions = True 
#instead of a clean phi, theta array choose random theta and phi positions - still sampled in sin(theta), phi
random_positions = False 
#Now true uniformity over the sphere (i.e. constant number of points per unit area)
uniform_positions = False

offset = 0 #1/6 * np.pi#add a bit of margin not to cut into the middle of a bin (the ones on the edge of the quarter)

#The limits for our theta and phi
#These are for a 'classic' quarter of the mPMT dome
theta_min = 0#0.59 #0 #0.59 #0
theta_max = 1.1 #1.1#0.61 #1.1 #0.61 #1.1 # #1 #now it is done in sin(theta) #np.pi/2 #up to 90 degrees for now
phi_min = 0 #np.pi/2 - offset/2 #3*np.pi/2 - have checked the symmetry: we can extrapolate a quarter to the full map 
phi_max = 2 * np.pi #np.pi + offset #2*np.pi  #up to 90 degrees for now -  second ring of bins is causing the problem 
#so we need a bit less offset

#This is for an intensive scan around a given bin (24)
#theta_min = 0.23445
#theta_max = 0.3917
#phi_min = 0.94205
#phi_max = 1.0899

#Set some default values
absff = 1.30
rayff = 0.75
FileID = 0

#Read the user inputs
argv = sys.argv[1:]
opts, args = getopt.getopt(argv, "t:R:p:f:e:a:r:d:u:")
for opt, arg in opts:
        if opt in ['-R']:
            R = float(arg)
        elif opt in ['-p']:
            nPhi = int(arg)
        elif opt in ['-t']:
            nTheta = int(arg)
        elif opt in ['-f']:
            FileID = int(arg)
        elif opt in ['-e']:
            nEvent = int(arg)
        elif opt in ['-a']:
            absff = float(arg)
        elif opt in ['-r']:
            rayff = float(arg)
        elif opt in ['-d']:
            random_positions = True
            ordered_positions = False
            nPositions = int(arg)
            print("careful - we are using random source positions instead of an array of %i pos in phi and %i pos in theta"%(nPhi, nTheta))
        elif opt in ['-u']:
            uniform_positions = True
            ordered_positions = False
            nPositions = int(arg)
            print("careful - we are using %i uniform source positions instead of an array of %i pos in phi and %i pos in theta"%(nPositions, nPhi, nTheta))

if ordered_positions == True:
    nPositions = nTheta * nPhi
#need a file name and print the config we are choosing
print("\n---------------------------------------------------------------------------------------------------------")
try:
    FileID
    print("File ID: %i" % FileID)
except NameError:
    print("Error: No Filename  detected")
print("Config:\n  R = %.2f \n  nEvent = %i \n  abwff = %.3e \n  rayff = %.3e"%(R, nEvent, absff, rayff))
print("Sampling method: \n  ordered = %s (%i x %i pos) \n  random = %s (%i pos) \n  uniform = %s (%i pos)"%(ordered_positions, nTheta, nPhi, random_positions, nPositions, uniform_positions, nPositions))
print("\n%.3f < theta < %.3f and %.3f < phi < %.3f"%(theta_min, theta_max, phi_min, phi_max))
print("---------------------------------------------------------------------------------------------------------\n")

#get the correct name for the file
alpha_mode = "Absff%.3e_Rayff%.3e"%(absff, rayff)
if absff>=10 and rayff <=10:
    alpha_mode = "Absff%.1e_Rayff%.3e"%(absff, rayff)
if rayff>=10 and absff <=10:
    alpha_mode = "Absff%.3e_Rayff%.1e"%(absff, rayff)
if rayff>=10 and absff >=10:
    alpha_mode = "Absff%.1e_Rayff%.1e"%(absff, rayff)

#First write the tuning file:
md.makeTuningConfigFile(FileID, absff, rayff)

########################## The make the config files : if ordered ######################################
if ordered_positions == True:
    range_theta = np.linspace(np.sin(theta_min), np.sin(theta_max), nTheta)
    range_phi = np.linspace(phi_min, phi_max, nPhi)
    range_theta = np.arcsin(range_theta)

    for theta in range_theta:
        w = 0
        for phi in range_phi:
            if theta!=0 or w == 0:   
              source_xpos, source_ypos, source_zpos = md.Rtp_to_xyz_source(R, theta, phi)
              md.makeConfigFile(source_xpos, source_ypos, source_zpos, alpha_mode,theta, phi, R, FileID, nEvent)

########################## Make the files : if random  ######################################
if random_positions == True:
    range_theta = []
    range_phi = []
    for k in range(int(nPositions)):
        range_theta.append(random.uniform(np.sin(theta_min), np.sin(theta_max))) #here random draw in sin theta
        range_phi.append(random.uniform(phi_min, phi_max))
    range_theta = np.arcsin(np.array(range_theta)) #convert back to theta
    range_phi = np.array(range_phi) #for now use completely random phi
    for t in range(len(range_theta)):
        phi = range_phi[t];
        theta = range_theta[t];
        if theta!=0 or w == 0: 
            source_xpos, source_ypos, source_zpos = md.Rtp_to_xyz_source(R, theta, phi)
            md.makeConfigFile(source_xpos, source_ypos, source_zpos, alpha_mode,theta, phi, R, FileID, nEvent)
            w = 1

########################## Make the files : if uniform ######################################
if uniform_positions == True:
    range_theta = []
    range_phi = []
    for k in range(nPositions):
        theta, phi = gp.get_tp_point(theta_min = theta_min, theta_max = theta_max, phi_min = phi_min, phi_max = phi_max)
        range_theta.append(float("%.3f"%theta))
        range_phi.append(float("%.3f"%phi))
    range_phi = np.array(range_phi) #for now use completely random phi
    range_theta = np.array(range_theta)
    print(range_phi)
    for t in range(len(range_theta)):
        phi = range_phi[t]
        theta = range_theta[t]
        source_xpos, source_ypos, source_zpos = md.Rtp_to_xyz_source(R, theta, phi)
        md.makeConfigFile(source_xpos, source_ypos, source_zpos, alpha_mode,theta, phi, R, FileID, nEvent)






