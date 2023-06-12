import read_data_results3 as rd
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import getopt
import pandas as pd
import os
import glob

#This is reading and combining all of the reference files where there is both no absorption and no scattering to make the reference absorption map where each line is the max number of pe collected per photon in this bin 

bin_number = sys.argv[1]

folder = "./Maps/ref_noAttenuation/OneBin*bin%s_R120.00*"%bin_number

all_test_folders =  glob.glob("%s"%folder)

if len(all_test_folders)>0:
    test_file  = all_test_folders[0]
    PMT_position_file = '%s' %(test_file)
    all_data  = np.array(rd.read_data3(PMT_position_file))
    i = 0
    count = 0
    photons = 0
    while (i<len(all_data[0])-1):
        count = count + all_data[6][i]
        photons = photons + all_data[7][i]
        i = i+1
    if photons > 0:
        print(count/photons)
    else:
        print(0)

else:#still need to add a line for bookkeping even when there is no photon recorded in this bin
    print(0)



