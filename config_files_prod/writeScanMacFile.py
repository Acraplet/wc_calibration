#This is a python file to write multiple occurences of a WCSim file
#this is similar to its mother code: write Macfile but it is only for four specific positions
#(0.00, 0.00)
#(0.49, 1.07)
#(0.33, 1.07)
#(0.92, 1.07)

#It runs with the following inputs:
# R - distance from the PMT to the source
# a - abwff the absorption coefficient in WCSim
# r - rayff the Rayleigh scattering coef in WCSim
# 

#import uproot
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import getopt
import random

#Here are the default values for R, nTheta, nPhi
R = 10
nEvent = 10000

#These are the given theta and phi positions over which we'll scan
range_theta = [0.0, 0.33, 0.49, 0.92]
range_phi = [0.0, 1.07, 1.07, 1.07]

random_positions = False 
#instead of a clean phi, theta array choose random theta and phi positions
theta_min = 0
theta_max = np.sin(1.1) # #1 #now it is done in sin(theta) #np.pi/2 #up to 90 degrees for now

phi_min = 0
phi_max = np.pi/2 #up to 90 degrees for now

argv = sys.argv[1:]

opts, args = getopt.getopt(argv, "t:R:p:f:e:a:r:d:")
for opt, arg in opts:
        if opt in ['-R']:
            R = float(arg)
        elif opt in ['-f']:
            FileID = int(arg)
        elif opt in ['-e']:
            nEvent = int(arg)
        elif opt in ['-a']:
            absff = float(arg)
        elif opt in ['-r']:
            rayff = float(arg)

#All the informations about the file are stored in its name
if absff>=10 and rayff <=10:
    alpha_mode = "Absff%.1e_Rayff%.3e"%(absff, rayff)
if rayff>=10 and abwff <=10:
    alpha_mode = "Absff%.3e_Rayff%.1e"%(absff, rayff)
if rayff>=10 and abwff >=10:
    alpha_mode = "Absff%.1e_Rayff%.1e"%(absff, rayff)

#need a file name
try:
    FileID
    print("File ID: %i" % FileID)
except NameError:
    print("Error: No Filename  detected")

print("----------------------------------------------------------------------------------------------------------")
print("Config:\n R = %.2f \n nEvent = %i \n abwff = %.3e \n rayff = %.3e"%(R, nEvent, absff, rayff))
print("Source positions: \n theta  = %s \n phi = %s"%(*range_theta, range_phi))
print("----------------------------------------------------------------------------------------------------------")


def makeConfigFile(source_xpos, source_ypos, source_zpos, alpha_mode, theta, phi, R, nEvent = 10000):

    run_beam_on = nEvent

    template_txtFile = open("WCSim_template.txt","r")
    data_saving_path = "/vols/t2k/users/ac4317/WCTE/WCSim/mPMTmapping/data"
    config_saving_path = "/home/ac4317/Laptops/Year1/WCTE/wc_calibration/WCSim_ScanConfigFiles"

    config_file_name = "%s/WCSim_config_mPMTmapping_401nm_FileID%i_%s_x%.2f_y%.2f_z%.2f_t%.2f_p%.2f_R%.2f.mac"%(config_saving_path, FileID,alpha_mode,source_xpos,source_ypos,source_zpos, theta, phi, R)
    data_file_name = "wcsim_mPMTmapping_401nm_FileID%i_%s_x%.2f_y%.2f_z%.2f_t%.2f_p%.2f_R%.2f.root"%(FileID, alpha_mode,source_xpos,source_ypos,source_zpos, theta, phi, R)
    #making sure the value we write down is the correct value we implement in the config file
    theta = float("%.2f"%theta)
    phi = float("%.2f"%phi)
    #print(theta, phi)
    source_dirx = - np.sin(theta) * np.cos(phi)
    source_diry = - np.cos(theta)
    source_dirz = - np.sin(theta) * np.sin(phi)
    orientation_string = "/gun/direction %.5f %.5f %.5f\n"%(source_dirx, source_diry, source_dirz) 
    #can be more precise once we save the source direction infp
    position_string = "/gun/position %.2f %.2f %.2f \n"%(source_xpos, source_ypos, source_zpos)

    with open("%s"%config_file_name, "w") as file:
        for line in template_txtFile:
            file.write(line)
        file.write(orientation_string)
        file.write(position_string + "\n")
        file.write("/WCSim/tuning/abwff %s\n"%absff)
        file.write("/WCSim/tuning/rayff %s\n"%rayff)
        file.write("/Tracking/fractionOpticalPhotonsToDraw 100.0 \n")
        file.write("/WCSimIO/RootFile %s \n"%(data_file_name))
        #now we are running in temp folder so no need to have the path to the file
        file.write("/WCSimIO/SaveRooTracker 0 \n")
        file.write("/run/beamOn %i"%run_beam_on)
    template_txtFile.close()



centre_offset = 27.4 #cm - distance from the PMT posiiton to the centre of the mPMT sphere

targetPMT_xpos = 0.0
targetPMT_ypos = -128.05 - centre_offset 
targetPMT_zpos = 0.0

source_xpos = 0
source_ypos = 0
source_zpos = 0
mPMT_radius = 34.2 
#need to put the source away from the centre of the mPMT sphere

for i in range(len(range_theta)):
    theta = range_theta[i]
    phi = range_phi[i]
    source_xpos = targetPMT_xpos + (R + mPMT_radius) * np.sin(theta) * np.cos(phi)
    source_ypos = targetPMT_ypos + (R + mPMT_radius) * np.cos(theta) 
    #careful, y is the "typical" vertical z coordinate
    source_zpos = targetPMT_zpos + (R + mPMT_radius) * np.sin(theta) * np.sin(phi)
    makeConfigFile(source_xpos, source_ypos, source_zpos, alpha_mode,theta, phi, R, nEvent)







