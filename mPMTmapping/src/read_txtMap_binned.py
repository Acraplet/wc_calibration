#This is a python file to read the maps and then show the data in pre-made bins
#This only works for a simgle file - no comparision possible at this stage
import read_data_results3 as rd
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import getopt
import pandas as pd
import os

filename = "No filename, please input one with the -f option"
comparision_file = 'none'
outputfile_name = 'none'
filename = []

# Get environment variables
WORKDIR = os.getenv('WCCALIB')

#plt.style.use(["science", "notebook", "grid"])

phi_max = 90 #for plotting - get the max phi value to show only the right portion of the circle
interpolation_mode = 'linear' # linear interpolation of the data
argv = sys.argv[1:]

#set up the geometry of the tank i.e. read the PMT positions
PMT_position_file = WORKDIR+'/mPMTmapping/PMT_positions.txt'
PMT = np.array(rd.read_data3(PMT_position_file)).T
PMT_mPMT = PMT[0]
PMT_mPMT_pmt = PMT[1]
PMT_x = PMT[2]
PMT_y = PMT[3]
PMT_z = PMT[4]
df_geom = pd.DataFrame()
for i in range(len(PMT[0])):
        c = [PMT[0][i],PMT[1][i],PMT[2][i],PMT[3][i],PMT[4][i]]
        row =  pd.Series(data=c, index=['mPMT', 'mPMT_pmt', 'x', 'y', 'z'], dtype=np.float64)
        df_geom = df_geom.append(row, ignore_index=True)

bins_position = WORKDIR+'/mPMTmapping/uniform_304_bins.txt'
bins = np.array(rd.read_data3(bins_position))

#Now look at the theta and phi - is a simple distance in theta and phi enough for the full distance?
bins_position = WORKDIR+'/mPMTmapping/uniform_top_bins_theta_phi.txt'
bins_theta = np.array(rd.read_data3(bins_position)).T

bins_position = WORKDIR+'/mPMTmapping/uniform_top_bins_withBinNumber.txt'
bins_all = np.array(rd.read_data3(bins_position))
#so this is a [[all x], [all y], [all z]] array with .T
#alternatively we can have an array
#[[x0,y0,z0], [x1, y1, z1], ...] which I think is better so that's bins without the .T at the end

#Read the user input
opts, args = getopt.getopt(argv, "m:o:p:f:c:a:")
for opt, arg in opts:
        if opt in ['-m']:
            phi_max = float(arg)
        elif opt in ['-f']:
            filename.append(str(arg))
        elif opt in ['-c']:
            filename.append(str(arg))
        elif opt in ['-o']:
            outputfile_name = str(arg)


if outputfile_name == 'none':
    print("No outputfile_name given - naming the output 'output.png'")
    outputfile_name = 'output'

#When we have a single file
df = pd.DataFrame()

for f in filename:
    table = rd.read_data3(f)
    print(f)
    source_xyz = [table [0], table[1], table [2]]
    source_Rtp = [table [3], table [4], table[5]]
    x, y, z = np.array(table[0]), np.array(table[1]), np.array(table[2])
    R, theta, phi = np.array(table[5]), np.array(table[3]), np.array(table[4])
    Q_tot =  np.array(table[6])
    nEvents = np.array(table[7])
    phi_max = 360
    for i in range(len(x)):
        c = [x[i], y[i],  z[i],  theta[i], phi[i],  R[i],  Q_tot[i],  nEvents[i]]
        row =  pd.Series(data=c, index=['x', 'y', 'z', 'theta', 'phi', 'R', 'Q', 'events'], dtype=np.float64)
        df = df.append(row, ignore_index=True)

#place each source position in a bin
count_uniform_bins = np.zeros(len(bins_all[0]))
count_uniform_bins_theta = np.zeros(len(bins_all[0]))
charge_uniform_bins = np.zeros(len(bins_all[0]))
charge_uniform_bins_theta = np.zeros(len(bins_all[0]))
for source_position_ID in range(len(df['x'])):
    source_position = np.array([df['x'][source_position_ID], df['y'][source_position_ID], df['z'][source_position_ID]])
    #Calculate the distance xyz bwteen the source position and each bin centre
    #distance = np.linalg.norm(bins - source_position, axis = 1)
    #id_of_closest_bin = np.argmin(distance)
    #now calculate the theta-phi distance 
    source_theta_phi = np.array([df['theta'][source_position_ID], df['phi'][source_position_ID]])
    distance_theta_phi = np.linalg.norm(np.array([bins_all[3], bins_all[4]]).T - source_theta_phi, axis = 1)
    id_of_closest_bin_theta = np.argmin(distance_theta_phi)
    #print(id_of_closest_bin_theta, id_of_closest_bin)
    #print(id_of_closest_bin_theta, bins_all[3][id_of_closest_bin_theta], bins_all[4][id_of_closest_bin_theta], source_theta_phi)

    #Store the data in the relevant bins
    #count_uniform_bins[id_of_closest_bin] += 1
    count_uniform_bins_theta[id_of_closest_bin_theta] += 1
    #charge_uniform_bins[id_of_closest_bin] += df['Q'][source_position_ID]
    charge_uniform_bins_theta[id_of_closest_bin_theta] += df['Q'][source_position_ID]
#print(count_uniform_bins_theta, charge_uniform_bins_theta)

ax = plt.axes(projection = 'polar')
ax.set_thetamin(0)
ax.set_ylabel("Q %s "%(filename))
ax.set_thetamax(phi_max)
for b in range(len(bins_all[0])):
    #if count_uniform_bins_theta[b]>0:
    #/( count_uniform_bins_theta[b]* nEvents[0]
    sc = ax.scatter(bins_all[4][b], bins_all[3][b], c = charge_uniform_bins_theta[b]/(count_uniform_bins_theta[b] * nEvents[0]), cmap = "Blues",vmin = 0, vmax = 0.33, s = 120)

plt.colorbar(sc)
ax.set_title('Charge collected per photon in this bin')
plt.savefig(WORKDIR+'/mPMTmapping/Maps/maps_pictures/Binned_visualisation/Non-interpolated/%s.png'%outputfile_name)
plt.show()
plt.close()
#print(count_uniform_bins_theta)

ax2 = plt.axes(projection = 'polar')
ax2.set_thetamin(0)
ax2.set_ylabel("Q %s "%(filename))
ax2.set_thetamax(phi_max)
for b in range(len(bins_all[0])):
    if count_uniform_bins_theta[b]>0:
    #/( count_uniform_bins_theta[b]* nEvents[0]
        print(bins_all[4][b], bins_all[3][b], count_uniform_bins_theta[b])
        sc2 = ax2.scatter(bins_all[4][b], bins_all[3][b], c = count_uniform_bins_theta[b], cmap = "Blues",vmin = 0, vmax = count_uniform_bins_theta.max(), s = 120)
plt.colorbar(sc2)
ax2.set_title('Number of source positions in this bin')
#plt.savefig('/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/Maps/maps_pictures/Binned_visualisation/Non-interpolated/%s.png'%outputfile_name)
plt.show()


ax = plt.axes(projection = '3d')
for b in range(len(bins_all[0])):
    #if count_uniform_bins_theta[b]>0:
        sc = ax.scatter(bins_all[0][b] ,  bins_all[2][b] ,bins_all[1][b], c = count_uniform_bins_theta[b], cmap = "Blues", vmin=0, vmax=count_uniform_bins_theta.max(), s = 120)
plt.colorbar(sc)
ax.set_title('Number of source positions in this bin')
plt.show()

ax = plt.axes(projection = '3d')
for b in range(len(bins_all[0])):
    #if count_uniform_bins_theta[b]>0:
        sc = ax.scatter(bins_all[0][b] ,  bins_all[2][b], bins_all[1][b], c = charge_uniform_bins_theta[b]/( count_uniform_bins_theta[b]* nEvents[0]), cmap = "Blues", vmin=0, vmax=0.4, s = 60)
plt.colorbar(sc)
ax.set_title('Charge collected per photon in this bin')
plt.show()




fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(projection='polar')
#only the points
ax.set_title("Q")
zr = df['Q']/df['events']
im = ax.scatter(df['phi'], df['theta'], c = zr, cmap = "nipy_spectral",s = 40)
col2 = plt.colorbar(im, label=f"Average number of *raw* P.E. in mPMT58 per photon", orientation="horizontal", format= "%.2f")
ax.set_thetamin(0)
ax.set_ylabel("Q %s "%(filename))
ax.set_thetamax(phi_max)
ax.set_xlabel(f'$\Theta$(rad)')
ax.xaxis.labelpad = 10
ax.yaxis.labelpad = 20
plt.savefig(WORKDIR+'/mPMTmapping/Maps/maps_pictures/Maps/Non-interpolated/%s.png'%outputfile_name)
plt.show()

#Now the interpolate
fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(projection='polar')
ax.set_title("Q")
ax.set_ylabel("Q %s "%(filename))
zr = df['Q']/10000
xi, yi = np.linspace(min(df['phi']), max(df['phi']), 100), np.linspace(min(df['theta']), max(df['theta']), 100)
xi, yi = np.meshgrid(xi, yi)
# Interpolate
rbf = scipy.interpolate.Rbf(df['phi'], df['theta'], zr, function=interpolation_mode, vmin = zr.min(), vmax = zr.max())
zi = rbf(xi, yi)
ax.contourf(xi, yi, zi, 500, cmap='nipy_spectral', vmin = zr.min(), vmax = zr.max())
#im = ax.scatter(phi, theta, c = zr, cmap = "nipy_spectral",s = 40)
col2 = plt.colorbar(im, label=f"Average number of *raw* P.E. in mPMT58 per photon", orientation="vertical", format= "%.2f")
ax.set_thetamin(0)
ax.set_thetamax(phi_max)
ax.xaxis.labelpad = 10
ax.yaxis.labelpad = 20
ax.set_xlabel(f'$\Theta$(rad)')
plt.savefig(WORKDIR+'/mPMTmapping/Maps/maps_pictures/Maps/Interpolated/%s.png'%outputfile_name)
plt.show()

#raise end
sys.exit()

#Now plot the source positions overlayed with the position of the maps
source_Rtp = np.array(source_Rtp)
source_xyz = np.array(source_xyz)
ax = plt.axes(projection ="3d")
#This is the mPMT (outer) shpere
r = 34.2
u, v = np.mgrid[0:2 * np.pi:30j, 0:np.pi:20j]
xs = r * np.cos(u) * np.sin(v)
ys = r * np.sin(u) * np.sin(v)
zs = r * np.cos(v) - 128.05 - 27.4
ax.plot_surface(xs, ys, zs, cmap=plt.cm.YlGnBu, alpha = 0.1)
#This is the mPMT inner sphere (10mm of thickness)
r = 33.2
u, v = np.mgrid[0:2 * np.pi:30j, 0:np.pi:20j]
xs = r * np.cos(u) * np.sin(v)
ys = r * np.sin(u) * np.sin(v)
zs = r * np.cos(v) - 128.05 - 27.4
ax.plot_surface(xs, ys, zs, cmap=plt.cm.YlGnBu, alpha = 0.1)
#this is the PMT sphere
r = 4.9
u, v = np.mgrid[0:2 * np.pi:30j, 0:np.pi:20j]
xs = r * np.cos(u) * np.sin(v)
ys = r * np.sin(u) * np.sin(v)
zs = r * np.cos(v) - 128.05
ax.plot_surface(xs, ys, zs, alpha = 0.1)
#This is the centre of the mPMT sphere
ax.scatter3D(0, 0, - 128.05 - 27.4, marker = 'x', color = 'k')
ax.scatter3D(x, z, y, marker = 'x', color ='r')
#these are the other PMTs - can be helpful for visualisation (uncomment)
#ax.scatter3D(PMT_x,PMT_z,PMT_y, marker = 'o')

#This is the position of the centre of the uniform bins
for b in range(len(bins_all)):
    ax.scatter3D(bins_all[0][b],bins_all[2][b], bins_all[1][b], 'o', color = 'k')

for i in range(len(x)):
    ax.plot([x[i], x[i] - R[i] *  np.sin(theta[i]) * np.cos(phi[i])], [z[i], z[i] - R[i] * np.sin(theta[i]) * np.sin(phi[i])], [y[i], y[i] - R[i] * np.cos(theta[i])], 'r-')
plt.show()





