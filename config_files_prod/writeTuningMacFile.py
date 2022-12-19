#This is a python file to write multiple occurences of a WCSim file

#It runs with the following inputs:
# R - distance from the PMT to the source
# nTheta - the number of theta points we want
# nPhi - the number of phi points we want

#to use: python writeMacFile.py -R 10 -t 10 -p 10 -f 10  

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

argv = sys.argv[1:]

opts, args = getopt.getopt(argv, "t:R:p:f:e:a:r:")
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

alpha_mode = "Absff%.3e_Rayff%.3e"%(absff, rayff)

if absff>=10:
    alpha_mode = "Absff%.1e_Rayff%.3e"%(absff, rayff)
if rayff>=10:
    alpha_mode = "Absff%.3e_Rayff%.1e"%(absff, rayff)

def makeTuningConfigFile(FileID, abwff=10e10, rayff = 10e10):
    template_txtFile = open("WCSim_tuning_template.txt","r")
    config_saving_path = "/home/ac4317/Laptops/Year1/WCTE/wc_calibration/WCSim_tuningFiles"

    config_file_name = "%s/tuning_parameters_FileID%i.mac"%(config_saving_path, FileID)
    with open("%s"%config_file_name, "w") as file:
        file.write("/WCSim/tuning/abwff %s\n"%absff)
        file.write("/WCSim/tuning/rayff %s\n"%rayff)
        for line in template_txtFile:
            file.write(line)
    #print(config_file_name)
    template_txtFile.close()

makeTuningConfigFile(FileID, absff, rayff)







