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
#import matplotlib.pyplot as plt
import sys
import getopt
import get_uniform_point as gp
import random
import Module_writeFiles as md
#Here are the default values for R, nTheta, nPhi
nEvent = 1000 #default number of events: 1000

#default source position: centre of the tank
x_source = 0
y_source = 0
z_source = 0

#Set some default values
absff = 1.30
rayff = 0.75
FileID = 0
coneTheta = 160

#Read the user inputs
argv = sys.argv[1:]
opts, args = getopt.getopt(argv, "t:R:p:t:f:e:a:r:d:u:x:y:z:")
for opt, arg in opts:
    if opt in ['-f']:
        FileID = int(arg)
    elif opt in ['-e']:
        nEvent = int(arg)
    elif opt in ['-a']:
        absff = float(arg)
    elif opt in ['-r']:
        rayff = float(arg)
    elif opt in ['-x']:
        x_source = float(arg)
    elif opt in ['-y']:
        y_source = float(arg)
    elif opt in ['-z']:
        z_source = float(arg)
    elif opt in ['-t']:
        coneTheta = float(arg)
    elif opt in ['-R']:
        #the number of photons - useful to keep
        R = int(arg)

#need a file name and print the config we are choosing
print("\n---------------------------------------------------------------------------------------------------------")
try:
    FileID
    print("File ID: %i" % FileID)
except NameError:
    print("Error: No Filename  detected")
print("Config: \n  number of photons: %i \n nEvent = %i \n  abwff = %.3e \n  rayff = %.3e"%(R, nEvent, absff, rayff))
print("Source position: \n [%.3f, %.3f, %.3f]"%(x_source, y_source, z_source))
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

########################## Make the files : with source positions  ######################################
flag = -9999

md.makeConfigFile_LB(x_source, y_source, z_source, alpha_mode,flag, flag, R, FileID, nEvent, coneTheta)

print("CAREFULL, the dark cone is %.2f degrees" %coneTheta)




