#This is a python file with the function to write the config and tuning WCSim file and other useful functions 

#import uproot
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import getopt
#import pandas as pd
#import uproot

# #Get environment variables
WORKDIR = os.getenv('WCCALIB')

target_mPMT = 58
txtFile = WORKDIR+"/mPMTmapping/mPMT_positionAndDirection_new.txt"
#print(txtFile)
with open(txtFile) as f:
	#a = f.readlines()
	for line in f.readlines():
		print(line.split(" "))
		mPMT = line.split(" ")[0]
		if int(float(mPMT)) == target_mPMT:
			targetPMT_xpos = float(line.split(" ")[1]) 
			targetPMT_ypos = float(line.split(" ")[2]) 
			targetPMT_zpos = float(line.split(" ")[3])
			mPMTdirx = float(line.split(" ")[4]) 
			mPMTdiry = float(line.split(" ")[5]) 
			mPMTdirz = float(line.split(" ")[6]) 

#df = pd.DataFrame()
#table = rd.read_data3(txtFile)
#
#df["mPMT"] = np.array(table[0])
#df["x"] = np.array(table[1])
#df["y"] = np.array(table[2])
#df["z"] = np.array(table[3])
#df["dirx"] = np.array(table[4])
#df["diry"] = np.array(table[5])
#df["dirz"] = np.array(table[6])
#
#df = df[df['mPMT']==target_mPMT]
#
##Geometrical factors -  to position the source above and around the correct mPMT
centre_offset = 27.4 #cm - distance from the PMT posiiton to the centre of the mPMT sphere
#targetPMT_xpos = float(df["x"]) #-158.94 #-1.583940000000000055e+02 #1.868740000000000236e+02  #0.0
#targetPMT_ypos = float(df["y"]) #-1.554500000000000171e+02 #-2.673999999999999844e+01#2.673999999999999844e+01 #-128.05 - centre_offset
#targetPMT_zpos = float(df["z"])#-105.98#0.000000000000000+02 #3.731000000000000227e+01 #0.0
#mPMTdirx = float(df["dirx"]) #0  #0.829999999999
#mPMTdiry = float(df["diry"]) #-0.000000000000000000e+00
#mPMTdirz = float(df["dirz"]) #5.600000000000000533e-01
#
print("mPMT %i \n x %.2f \n y %.2f \n z %.2f \n dirx %.2f \n diry %.2f \n dirz %.2f"%(target_mPMT, targetPMT_xpos, targetPMT_ypos, targetPMT_zpos, mPMTdirx, mPMTdiry,mPMTdirz) )

#raise end

ref_dirx = 0
ref_diry = 1
ref_dirz = 0
mPMTdir_vector = np.array([mPMTdirx, mPMTdiry, mPMTdirz])
refdir_vector = np.array([ref_dirx, ref_diry, ref_dirz])


#need to calculate the rotation matrix from these 
rotation_axis = np.cross(refdir_vector, mPMTdir_vector)
#print(rotation_axis)
angle_of_rotation = np.arccos(np.dot(refdir_vector, mPMTdir_vector)) #* np.sign(np.dot(refdir_vector, mPMTdir_vector))
v = rotation_axis
skew_matrix = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
#print(skew_matrix)
rotation_matrix = np.identity(3) + np.sin(angle_of_rotation) * skew_matrix + (1- np.dot(refdir_vector, mPMTdir_vector)) * np.matmul(skew_matrix, skew_matrix)
if ref_diry == -mPMTdiry:
	rotation_matrix[1][1] = -rotation_matrix[1][1]

print("rotation matrix: ", rotation_matrix, "\nThis is the should be mPMT dir: ", np.matmul(rotation_matrix, refdir_vector), "\nThis is the true mPMT vector: ", mPMTdir_vector)


source_xpos = 0
source_ypos = 0
source_zpos = 0
mPMT_radius = 34.2 #2 * 34.2 #34.2 #need to put the source away from the centre of the mPMT sphere
#source_offset_from_dome = 0

#to calculate AR:
#targetPMT_xpos, targetPMT_ypos, targetPMT_zpos = targetPMT_xpos + mPMTdirx * mPMT_radius, targetPMT_ypos + mPMTdiry * mPMT_radius, targetPMT_zpos + mPMTdirz * mPMT_radius

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
    source_xpos = np.sin(theta) * np.cos(phi) #* np.sign(phi) #(R + mPMT_radius)
    source_ypos = np.cos(theta)
    source_zpos = np.sin(theta) * np.sin(phi) #* np.sin(phi-np.pi/2)
    #rotate w.r.t the mPMT direction
    array_source_pos_rotated = np.matmul(rotation_matrix, np.array([source_xpos, source_ypos, source_zpos])) 	
    source_xpos = array_source_pos_rotated[0]
    source_ypos = array_source_pos_rotated[1]
    source_zpos = array_source_pos_rotated[2]
    source_xpos = (source_xpos * (R + mPMT_radius) + targetPMT_xpos)
    source_ypos = (source_ypos * (R + mPMT_radius) + targetPMT_ypos)
    source_zpos = (source_zpos * (R + mPMT_radius) + targetPMT_zpos)
    print(source_xpos, source_ypos, source_zpos, 'source pos') 
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
    source_dirx = -np.sin(theta) * np.cos(phi) #* np.sign(phi)
    source_diry = -np.cos(theta)
    source_dirz = -np.sin(theta) * np.sin(phi) #* np.sin(phi-np.pi/2)
    #norm = np.sqrt(source_dirx*source_dirx + source_diry * source_diry + source_dirz * source_dirz)
    #source_dirx = source_dirx/norm
    #source_diry = source_diry/norm
    #source_dirz = source_dirz/norm
    array_source_dir_rotated = np.matmul(rotation_matrix, np.array([source_dirx, source_diry, source_dirz])) 	
    source_dirx = array_source_dir_rotated[0]
    source_diry = array_source_dir_rotated[1]
    source_dirz = array_source_dir_rotated[2]


    #print("The direction: ", array_source_dir_rotated)

    orientation_string = "/gun/direction %.5f %.5f %.5f\n"%(source_dirx, source_diry, source_dirz)
    #can be more precise once we save the source direction infp
    position_string = "/gun/position %.3f %.3f %.3f \n"%(source_xpos, source_ypos, source_zpos)

    with open("%s"%config_file_name, "w") as file:
        for line in template_txtFile:
            file.write(line)
        file.write(orientation_string)
        file.write(position_string + "\n")
        #file.write("/Tracking/fractionOpticalPhotonsToDraw 100.0 \n")
        file.write("/WCSimIO/RootFile %s \n"%(data_file_name))
        #now we are running in temp folder so no need to have the path to the file
        file.write("/WCSimIO/SaveRooTracker 0 \n")
        file.write("/run/beamOn %i"%run_beam_on)
    template_txtFile.close()
    #save the information about the source position and direction for debugging purposes
    writeOutSourceInfo(source_xpos, source_ypos, source_zpos, source_dirx, source_diry, source_dirz, FileID)

def makeConfigFile_LB(source_xpos, source_ypos, source_zpos, alpha_mode, theta, phi, R, FileID, nEvent = 10000):
    ''' 
    This function write a config file from the template we have and with the given inputs
    Craefull - now we have a 3 s.f. precision on the source position -> this will have to be updated everywhere
    Unlike before, now R is hte number of repetitions, i.e. the number of photons produced by the LB
    '''
    run_beam_on = nEvent
    template_txtFile = open(WORKDIR+"/config_files_prod/LaserBall_template.txt","r")
    data_saving_path = WORKDIR+"/mPMTmapping/data"
    config_saving_path = WORKDIR+"/WCSim_configFiles"
	#"/vols/t2k/users/ac4317/WCTE/WCSim/mPMTmapping/config_files"
    config_file_name = "%s/WCSim_config_mPMTmapping_401nm_FileID%i_%s_x%.3f_y%.3f_z%.3f_t%.3f_p%.3f_R%i.mac"%(config_saving_path, FileID,alpha_mode,source_xpos,source_ypos,source_zpos, theta, phi, R)
    data_file_name = "%s/wcsim_mPMTmapping_401nm_FileID%i_%s_x%.3f_y%.3f_z%.3f_t%.3f_p%.3f_R%i.root"%(data_saving_path,FileID, alpha_mode,source_xpos,source_ypos,source_zpos, theta, phi, R)
    #can be more precise once we save the source direction info
    position_string = "/gps/position %.3f %.3f %.3f \n"%(source_xpos, source_ypos, source_zpos)
    numberOfPhoton_string = "/gps/number %i \n"%(R)

    with open("%s"%config_file_name, "w") as file:
        for line in template_txtFile:
            file.write(line)
        file.write(numberOfPhoton_string)
        file.write(position_string + "\n")
        file.write("/Tracking/fractionOpticalPhotonsToDraw 100.0 \n")
        file.write("/WCSimIO/RootFile %s \n"%(data_file_name))
        #now we are running in temp folder so no need to have the path to the file
        file.write("/WCSimIO/SaveRooTracker 0 \n")
        file.write("/run/beamOn %i"%run_beam_on)
    template_txtFile.close()

def writeOutSourceInfo(source_xpos, source_ypos, source_zpos, source_dirx, source_diry, source_dirz, FileID):
	with open("sourceInfo_FileID%i.txt"%(FileID), 'a') as file:
		file.write("%.3f %.3f %.3f %.3f %.3f %.3f %i \n"%(source_xpos, source_ypos, source_zpos, source_dirx, source_diry, source_dirz, FileID))
	return 0
 
