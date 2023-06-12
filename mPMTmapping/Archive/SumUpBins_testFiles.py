import read_data_results3 as rd
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import getopt
import pandas as pd
import os
import glob

bin_number = sys.argv[1]

folder = "./Maps/ref_maps_oneBin/OneBin*bin%s*"%bin_number

all_test_folders =  glob.glob("%s"%folder)

for test_file in all_test_folders[:-1]:
    #print(test_file)
    PMT_position_file = '%s' %(test_file)
    #plt.style.use(["science", "notebook", "grid"])
    all_data  = np.array(rd.read_data3(PMT_position_file))
    #print(PMT_position_file)
    #print(all_data[6])
    rayff = all_data[-1]
    current_rayff = rayff[0]
    i = 0
    while (i<len(all_data[0])-1):
        count = 0
        photons = 0
        while (rayff[i] == current_rayff and i<len(all_data[0])-1):
            count = count + all_data[6][i]
            photons = photons + all_data[7][i]
            R = all_data[5][i]
            i = i+1
        #print("File ID: ", test_file, "hits per 1000 photons: %.2f"%(count/photons * 1000), "number of source positions simulated", photons/1000, "true rayff: ", current_rayff, "R = ", R)
        print(current_rayff, R, count/photons * 1000)
        current_rayff = rayff[i]





