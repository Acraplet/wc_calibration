#This is a python file with the function to write the config and tuning WCSim file and other useful functions 

import uproot
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import getopt

#Geometrical factors -  to position the source above and around the correct mPMT
centre_offset = 27.4 #cm - distance from the PMT posiiton to the centre of the mPMT sphere
targetPMT_xpos = 0.0
targetPMT_ypos = -128.05 - centre_offset
targetPMT_zpos = 0.0
source_xpos = 0
source_ypos = 0
source_zpos = 0
mPMT_radius = 34.2 #need to put the source away from the centre of the mPMT sphere


def makeTuningConfigFile(FileID, abwff=10e10, rayff = 10e10):
    template_txtFile = open("WCSim_tuning_template.txt","r")
    config_saving_path = "/home/ac4317/Laptops/Year1/WCTE/wc_calibration/WCSim_tuningFiles"

    config_file_name = "%s/tuning_parameters_FileID%i.mac"%(config_saving_path, FileID)
    with open("%s"%config_file_name, "w") as file:
        file.write("/WCSim/tuning/abwff %s\n"%abwff)
        file.write("/WCSim/tuning/rayff %s\n"%rayff)
        for line in template_txtFile:
            file.write(line)
    #print(config_file_name)
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

def makeConfigFile(source_xpos, source_ypos, source_zpos, alpha_mode, theta, phi, R, FileID, nEvent = 10000):
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
        file.write("/Tracking/fractionOpticalPhotonsToDraw 100.0 \n")
        file.write("/WCSimIO/RootFile %s \n"%(data_file_name))
        #now we are running in temp folder so no need to have the path to the file
        file.write("/WCSimIO/SaveRooTracker 0 \n")
        file.write("/run/beamOn %i"%run_beam_on)
    template_txtFile.close()


