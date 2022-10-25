#This is a python file to write multiple occurences of a WCSim file

#It runs with the following inputs:
# R - distance from the PMT to the source
# nTheta - the number of theta points we want
# nPhi - the number of phi points we want

#the total number of points is therefore nTheta*nPhi

import uproot
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import getopt

###first up - get the position of the PMT we are interested in - not needed anymore
#output_checks_file = str(sys.argv[1])
#tree_data = uproot.open("%s"%output_checks_file)["pmt_type0"]
#df_data = tree_data.arrays(library="pd")

#df_data = df_data[df_data["mPMT_id"]==58]
#df_data = df_data[df_data["mPMT_pmt_id"]==18]

#print("x, y,z position of the central PMT of the central bottom mPMT module: ", df_data["xpos"], "\n",df_data["ypos"], "\n",df_data["zpos"])

#Here are the default values for R, nTheta, nPhi
nTheta = 10
nPhi = 0
R = 10

theta_min = 0
theta_max = np.pi #up to 180 degrees for now

phi_min = 0
phi_max = np.pi #up to 180 degrees for now

FileID = 0

argv = sys.argv[1:]

opts, args = getopt.getopt(argv, "t:R:p:f:")
for opt, arg in opts:
        if opt in ['-R']:
            R = float(arg)
        elif opt in ['-p']:
            nPhi = int(arg)
        elif opt in ['-t']:
            nTheta = int(arg)
        elif opt in ['-f']:
            FileID = int(arg)




def makeConfigFile(source_xpos, source_ypos, source_zpos, alpha_mode,theta, phi, R):

    run_beam_on = 100


    template_txtFile = open("WCSim_template.txt","r")
    data_saving_path = "/vols/t2k/users/ac4317/WCTE/WCSim/mPMTmapping/data"
    config_saving_path = "/home/ac4317/Laptops/Year1/WCTE/wc_calibration/WCSim_configFiles"

    config_file_name = "%s/WCSim_config_mPMTmapping_401nm_FileID%i_%s_x%.2f_y%.2f_z%.2f_t%.2f_p%.2f_R%.2f.mac"%(config_saving_path, FileID,alpha_mode,source_xpos,source_ypos,source_zpos, theta, phi, R)
    data_file_name = "wcsim_mPMTmapping_401nm_%s_x%.2f_y%.2f_z%.2f_t%.2f_p%.2f_R%.2f.root"%(alpha_mode,source_xpos,source_ypos,source_zpos, theta, phi, R)

    position_string = "/gun/position %.2f %.2f %.2f \n"%(source_xpos, source_ypos, source_zpos)

    with open("%s"%config_file_name, "w") as file:
        for line in template_txtFile:
            file.write(line)
        file.write(position_string + "\n")
        file.write("/Tracking/fractionOpticalPhotonsToDraw 100.0 \n")
        file.write("/WCSimIO/RootFile %s/%s \n"%(data_saving_path, data_file_name))
        file.write("/WCSimIO/SaveRooTracker 0 \n")
        file.write("/run/beamOn %i"%run_beam_on)
    template_txtFile.close()


print(R)

targetPMT_xpos = 0.0
targetPMT_ypos = -128.05
targetPMT_zpos = 0.0

source_xpos = 0
source_ypos = 0
source_zpos = 0
alpha_mode = "noAlpha"

range_theta = np.linspace(theta_min, theta_max, nTheta)
range_phi = np.linspace(phi_min, phi_max, nPhi)

for theta in range_theta:
    for phi in range_phi:
        source_xpos = targetPMT_xpos + R * np.sin(theta) * np.cos(phi)
        source_ypos = targetPMT_ypos + R * np.cos(theta) #careful, y is the "typical" vertical z coordinate
        source_zpos = targetPMT_zpos + R * np.sin(theta) * np.sin(phi)
        makeConfigFile(source_xpos, source_ypos, source_zpos, alpha_mode,theta, phi, R)









