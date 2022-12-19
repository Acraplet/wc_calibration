#This is a python file to write multiple occurences of a WCSim file

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

import uproot
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import getopt
import get_uniform_point as gp
import random
#Here are the default values for R, nTheta, nPhi
nTheta = 10
nPhi = 0
R = 10
nEvent = 10000

#Works naively  with neatly arranged (nTheta, nPhi) positions 
ordered positions = True 
#instead of a clean phi, theta array choose random theta and phi positions - still sampled in sin(theta), phi
random_positions = False 
#Now true uniformity over the sphere (i.e. constant number of points per unit area)
uniform_positions = False

#The limits for our theta and phi
theta_min = 0
theta_max = np.sin(1.1) # #1 #now it is done in sin(theta) #np.pi/2 #up to 90 degrees for now
phi_min = 0
phi_max = np.pi/2 #up to 90 degrees for now

#Geometrical factors -  to position the source above and around the correct mPMT
centre_offset = 27.4 #cm - distance from the PMT posiiton to the centre of the mPMT sphere
targetPMT_xpos = 0.0
targetPMT_ypos = -128.05 - centre_offset 
targetPMT_zpos = 0.0
source_xpos = 0
source_ypos = 0
source_zpos = 0
mPMT_radius = 34.2 #need to put the source away from the centre of the mPMT sphere


#Read the user inputs
argv = sys.argv[1:]
opts, args = getopt.getopt(argv, "t:R:p:f:e:a:r:d:")
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



#need a file name and print the config we are choosing
print("----------------------------------------------------------------------------------------------------------")
try:
    FileID
    print("File ID: %i" % FileID)
except NameError:
    print("Error: No Filename  detected")
print("Config:\n R = %.2f \n nEvent = %i \n abwff = %.3e \n rayff = %.3e"%(R, nEvent, absff, rayff))
print("Sampling method: \n ordered = %s (%i x %i pos) \n random (in sin(theta), phi)  = %s (%i pos) \n uniform = %s (%i pos)"%(ordered_positions, nTheta, nPhi, random_positions, nPositions, uniform_positions, nPositions))
print("----------------------------------------------------------------------------------------------------------")

#get the correct name for the file
if absff>=10 and rayff <=10:
    alpha_mode = "Absff%.1e_Rayff%.3e"%(absff, rayff)
if rayff>=10 and abwff <=10:
    alpha_mode = "Absff%.3e_Rayff%.1e"%(absff, rayff)
if rayff>=10 and abwff >=10:
    alpha_mode = "Absff%.1e_Rayff%.1e"%(absff, rayff)


def makeConfigFile(source_xpos, source_ypos, source_zpos, alpha_mode, theta, phi, R, nEvent = 10000):
    ''' 
    This function write a config file from the template we have and with the given inputs
    '''
    run_beam_on = nEvent
    template_txtFile = open("WCSim_template.txt","r")
    data_saving_path = "/vols/t2k/users/ac4317/WCTE/WCSim/mPMTmapping/data"
    config_saving_path = "/home/ac4317/Laptops/Year1/WCTE/wc_calibration/WCSim_configFiles"
    config_file_name = "%s/WCSim_config_mPMTmapping_401nm_FileID%i_%s_x%.2f_y%.2f_z%.2f_t%.2f_p%.2f_R%.2f.mac"%(config_saving_path, FileID,alpha_mode,source_xpos,source_ypos,source_zpos, theta, phi, R)
    data_file_name = "wcsim_mPMTmapping_401nm_FileID%i_%s_x%.2f_y%.2f_z%.2f_t%.2f_p%.2f_R%.2f.root"%(FileID, alpha_mode,source_xpos,source_ypos,source_zpos, theta, phi, R)
    #making sure the value we write down is the correct value we implement in the config file
    theta = float("%.2f"%theta)
    phi = float("%.2f"%phi)
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

def Rtp_to_xyz_source(R, theta, phi):
    '''
    extract the source positions corresponding to the given R theta phi coordinates
    '''
    global mPMT_radius, targetPMT_xpos, targetPMT_ypos, target_zpos
    source_xpos = targetPMT_xpos + (R + mPMT_radius) * np.sin(theta) * np.cos(phi)
    source_ypos = targetPMT_ypos + (R + mPMT_radius) * np.cos(theta)
    source_zpos = targetPMT_zpos + (R + mPMT_radius) * np.sin(theta) * np.sin(phi)
    return source_xpos, source_ypos, source_zpos 


########################## Make the files : if ordered ######################################
if ordered_positions == True:
    range_theta = np.linspace(theta_min, theta_max, nTheta)
    range_phi = np.linspace(phi_min, phi_max, nPhi)
    range_theta = np.arcsin(range_theta)

    for theta in range_theta:
        w = 0
        for phi in range_phi:
            if theta!=0 or w == 0:   
              source_xpos, source_ypos, source_zpos = Rtp_to_xyz_source(R, theta, phi)
              makeConfigFile(source_xpos, source_ypos, source_zpos, alpha_mode,theta, phi, R, nEvent)

########################## Make the files : if random  ######################################
if random_positions == True:
    range_theta = []
    range_phi = []
    for k in range(int(nPositions)):
        range_theta.append(random.uniform(theta_min, theta_max)) #here random draw in sin theta
        range_phi.append(random.uniform(phi_min, phi_max))
    range_theta = np.arcsin(np.array(range_theta)) #convert back to theta
    range_phi = np.array(range_phi) #for now use completely random phi
    for t in range(len(range_theta)):
        phi = range_phi[t];
        theta = range_theta[t];
        if theta!=0 or w == 0: 
            source_xpos, source_ypos, source_zpos = Rtp_to_xyz_source(R, theta, phi)
            makeConfigFile(source_xpos, source_ypos, source_zpos, alpha_mode,theta, phi, R, nEvent)
            w = 1

########################## Make the files : if uniform ######################################
if uniform_positions == True:
    range_theta = []
    range_phi = []
    for k in range(nPositions):
        theta, phi = gp.get_tp_point(theta_max = theta_max, phi_min = phi_min, phi_max = phi_max)
        range_theta.append(float("%.2f"%theta))
        range_phi.append(float("%.2f"%phi))
    range_phi = np.array(range_phi) #for now use completely random phi
    range_theta = np.array(range_theta)
    for t in range(len(range_theta)):
        phi = range_phi[t]
        theta = range_theta[t]
        source_xpos, source_ypos, source_zpos = Rtp_to_xyz_source(R, theta, phi)
        makeConfigFile(source_xpos, source_ypos, source_zpos, alpha_mode,theta, phi, R, nEvent)
        wt.makeTuningConfigFile(FileID,abwff = absff, rayff = rayff)






