#This is a python file with the function to write the config and tuning WCSim file and other useful functions 

#import uproot
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

# Get environment variables
WORKDIR = os.getenv('WCCALIB')

def makeTuningConfigFile(FileID, abwff=10e10, rayff = 10e10):
    template_txtFile = open(WORKDIR+"/config_files_prod/WCSim_tuning_template.txt","r")
    config_saving_path = WORKDIR+"/WCSim_tuningFiles"
#"/vols/t2k/users/ac4317/WCTE/WCSim/mPMTmapping/tuning_mac"

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
    source_xpos = targetPMT_xpos + (R + mPMT_radius) * np.sin(theta) * np.cos(phi) #* np.sign(phi)
    source_ypos = targetPMT_ypos + (R + mPMT_radius) * np.cos(theta)
    source_zpos = targetPMT_zpos + (R + mPMT_radius) * np.sin(theta) * np.sin(phi) #* np.sin(phi-np.pi/2)
    return source_xpos, source_ypos, source_zpos

def makeConfigFile(source_xpos, source_ypos, source_zpos, alpha_mode, theta, phi, R, FileID, nEvent = 10000):
    ''' 
    This function write a config file from the template we have and with the given inputs
    Craefull - now we have a 3 s.f. precision on the source position -> this will have to be updated everywhere
    '''
    run_beam_on = nEvent
    template_txtFile = open(WORKDIR+"/config_files_prod/WCSim_template.txt","r")
    data_saving_path = WORKDIR+"/mPMTmapping/data"
    config_saving_path = WORKDIR+"/WCSim_configFiles"
	#"/vols/t2k/users/ac4317/WCTE/WCSim/mPMTmapping/config_files"
    config_file_name = "%s/WCSim_config_mPMTmapping_401nm_FileID%i_%s_x%.3f_y%.3f_z%.3f_t%.3f_p%.3f_R%.2f.mac"%(config_saving_path, FileID,alpha_mode,source_xpos,source_ypos,source_zpos, theta, phi, R)
    data_file_name = "%s/wcsim_mPMTmapping_401nm_FileID%i_%s_x%.3f_y%.3f_z%.3f_t%.3f_p%.3f_R%.2f.root"%(data_saving_path,FileID, alpha_mode,source_xpos,source_ypos,source_zpos, theta, phi, R)
    #making sure the value we write down is the correct value we implement in the config file
    theta = float("%.3f"%theta)
    phi = float("%.3f"%phi)
    source_dirx = - np.sin(theta) * np.cos(phi) #* np.sign(phi)
    source_diry = - np.cos(theta)
    source_dirz = - np.sin(theta) * np.sin(phi) #* np.sin(phi-np.pi/2)
    orientation_string = "/gun/direction %.5f %.5f %.5f\n"%(source_dirx, source_diry, source_dirz)
    #can be more precise once we save the source direction infp
    position_string = "/gun/position %.3f %.3f %.3f \n"%(source_xpos, source_ypos, source_zpos)

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


