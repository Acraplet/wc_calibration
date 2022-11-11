#This is a python file to read the maps in txt format and do the plotting faster
import read_data_results3 as rd
#read the .txt file:
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import getopt
import pandas as pd

filename = "No filename, please input one with the -f option"
comparision_file = 'none'
outputfile_name = 'none'

#plt.style.use(["science", "notebook", "grid"])

#filename = str(sys.argv[1])
phi_max = 90 #for plotting - get the max phi value to show only the right portion of the circle
argv = sys.argv[1:-1]
j = 0

#careful, the last input is always the output file name
outputfile_name = sys.argv[-1]


color = ['blue', 'green', 'black', 'lightgray', 'red', 'orange', 'cyan']

fig = plt.figure(figsize = (15, 8))
ax = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
for filename in argv:
    filename = str(filename)
    df = pd.DataFrame()
    table = rd.read_data3(filename)
    source_xyz = [table [0], table[1], table [2]]
    source_Rtp = [table [3], table [4], table[5]]

    x, y, z = np.array(table[0]), np.array(table[1]), np.array(table[2])
    R, theta, phi = np.array(table[5]), np.array(table[3]), np.array(table[4])

    phi_max = max(phi) * 180/np.pi

    Q_tot = np.array(table[6])
    nEvents = np.array(table[7])

    for i in range(len(x)):
        c = [x[i], y[i],  z[i],  theta[i], phi[i],  R[i],  Q_tot[i],  nEvents[i]]
        row =  pd.Series(data=c, index=['x', 'y', 'z', 'theta', 'phi', 'R', 'Q', 'events'], dtype=np.float64)
        df = df.append(row, ignore_index=True)

    for t in df["theta"].unique():
        df_buf = df[df["theta"]==t]
        ax.plot(t, (df_buf["Q"]/df_buf["events"]).mean(), 's', color = color[j])
        if len(df_buf)!=1:
            ax2.plot(t, (df_buf["Q"]/df_buf["events"]).std(), 's', color = color[j])
        else: #no std for the single point at the centre
            ax2.plot(t, 0, 's', color = color[j])

    #so the legend doesn't appear twenty times
    ax.plot(t, (df_buf["Q"]/df_buf["events"]).mean(), 's', color = color[j], label = '%s'%filename)
    ax2.plot(t, (df_buf["Q"]/df_buf["events"]).std(), 's', color = color[j], label = '%s'%filename)
    j += 1




ax.set_title('%s'%outputfile_name)
ax.set_ylabel('MEAN charge collected \n per photon over pi/2 in phi')
ax.set_xlabel('Theta')
ax2.set_ylabel('STD of charge collected \n per photon over pi/2 in phi')
ax2.set_xlabel('Theta')
ax.grid()
ax2.grid()
ax2.legend()
ax.legend()
plt.savefig('/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_pictures/Comparisions/Averaged-over-phi/%s.png'%outputfile_name)
plt.show()




