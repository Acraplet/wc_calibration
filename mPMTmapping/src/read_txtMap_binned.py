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
PMT_position_file = WORKDIR+'/mPMTmapping/PMT_positions_with_direction.txt'
PMT = np.array(rd.read_data3(PMT_position_file)).T
PMT_mPMT = PMT[0]
PMT_mPMT_pmt = PMT[1]
PMT_x = PMT[2]
PMT_y = PMT[3]
PMT_z = PMT[4]
#test
PMT_dirx = PMT[5]
PMT_diry = PMT[6]
PMT_dirz = PMT[7]

print(PMT_diry)
df_geom = pd.DataFrame()
for i in range(len(PMT[0])):
        c = [PMT[0][i],PMT[1][i],PMT[2][i],PMT[3][i],PMT[4][i], PMT[5][i], PMT[6][i], PMT[7][i]]
        row =  pd.Series(data=c, index=['mPMT', 'mPMT_pmt', 'x', 'y', 'z', 'dirx', 'diry', 'dirz'], dtype=np.float64)
        df_geom = df_geom.append(row, ignore_index=True)

#Here check the cone of sight of the PMT
df_58 = df_geom[df_geom['mPMT']==58]
df_58 = df_58[df_58['mPMT_pmt']!=19]
df_59_0 = df_geom[df_geom['mPMT']==59]
df_59_0 = df_59_0[df_59_0['mPMT_pmt']==19]
df_58 = df_58.append(df_59_0, ignore_index = True)

print(df_58)

x_centre = 0
z_centre = 0
y_centre = -128.05-27.4
R_mPMT = - y_centre

df_58['OP_x'] = df_58['x'] - x_centre #teh centre of the bottom mPMT
df_58['OP_y'] = df_58['y'] - y_centre
df_58['OP_z'] = df_58['z'] - z_centre

norm = np.sqrt(df_58['OP_x']**2 + df_58['OP_y'] ** 2 + df_58['OP_z'] ** 2)
df_58['OP_x'] = df_58['OP_x']#/norm
df_58['OP_y'] = df_58['OP_y']#/norm
df_58['OP_z'] = df_58['OP_z']#/norm

norm_dirPMT = np.sqrt(df_58['x']**2 + df_58['y']**2 + df_58['z']**2)
df_58['x'] = df_58['x'] #/norm_dirPMT
df_58['y'] = df_58['y'] #/norm_dirPMT
df_58['z'] = df_58['z'] # /norm_dirPMT


#norm_dirPMT = np.sqrt(df_58['dirx']**2 + df_58['diry']**2 + df_58['dirz']**2)
df_58['dirx'] = df_58['dirx']
df_58['diry'] = df_58['diry']
df_58['dirz'] = df_58['dirz']

r = 34.2
u, v = np.mgrid[0: 2 * np.pi:30j, 0: 0.04 * np.pi:20j]


ax = plt.axes(projection = '3d')
for i in range(len(df_58)):
    # ax.plot([df_58['x'].iloc[i], df_58['x'].iloc[i]+df_58['dirx'].iloc[i]], [df_58['z'].iloc[i], df_58['z'].iloc[i]+df_58['dirz'].iloc[i]], [df_58['y'].iloc[i], df_58['y'].iloc[i]+df_58['diry'].iloc[i]], 'k-')

    # ax.plot([x_centre, x_centre+df_58['OP_x'].iloc[i] * R_mPMT], [z_centre, z_centre+df_58['OP_z'].iloc[i] * R_mPMT], [y_centre, y_centre+df_58['OP_y'].iloc[i] * R_mPMT], 'b--')

    theta_centre = np.arccos(-np.dot([df_58['OP_x'].iloc[i]/norm.iloc[i], df_58['OP_z'].iloc[i]/norm.iloc[i], df_58['OP_y'].iloc[i]/norm.iloc[i]], [df_58['x'].iloc[i]/norm_dirPMT.iloc[i], df_58['z'].iloc[i]/norm_dirPMT.iloc[i], df_58['y'].iloc[i]/norm_dirPMT.iloc[i]]))


    phi_centre = np.arccos(df_58['OP_x'].iloc[i]/(norm.iloc[i] * np.sin(theta_centre)) ) * np.sign(df_58['OP_z'].iloc[i])

    print(theta_centre, phi_centre)

    if np.sin(theta_centre) == 0:
        phi_centre = 0
    half_angle = 0.003

    # theta = np.linspace(0, , 100)
    # phi = np.linspace(phi_centre-half_angle, phi_centre+half_angle, 100)

    # y_shade = np.sqrt(half_angle**2 - (x_shade)**2) * np.sign(x_shade-theta_centre)
    #
    # theta = np.linspace(0, 2* np.pi, 100)
    # y = np.cos(theta) #* np.sign(x)
    # x = np.sin(theta)
    if np.sign(df_58['x'].iloc[i]) == 0:
        signx = 1
    else:
        signx = np.sign(df_58['x'].iloc[i])

    if np.sign(df_58['z'].iloc[i]) == 0:
        signz = 1
    else:
        signx = np.sign(df_58['z'].iloc[i])


    ax.scatter(np.cos(phi_centre) * np.sin(theta_centre) * 10 * signx, np.sin(phi_centre)*np.sin(theta_centre)* signz, np.cos(theta_centre))

    # ax.plot(x_shade, y_shade)
    # ax.plot(x, y)
    # ax.plot(np.cos(x_shade)*np.sin(y_shade), np.sin(x_shade)*np.sin(y_shade))

    # u, v = np.mgrid[phi_centre - half_angle: phi_centre + half_angle:30j, theta_centre - half_angle: theta_centre + half_angle:20j]
    # xs = r * np.cos(y) * np.sin(x)
    # ys = r * np.sin(y) * np.sin(x)
    # zs = r * np.cos(x) - 155.45
    # ax.scatter(xs, ys, zs, cmap=plt.cm.YlGnBu, alpha = 0.1)


    # u, v = np.mgrid[0: 2 * np.pi:30j, theta_centre - half_angle: theta_centre + half_angle:20j]
    # xs = r * np.cos(u) * np.sin(v)
    # ys = r * np.sin(u) * np.sin(v)
    # zs = r * np.cos(v) - 128.05 - 27.4
    # ax.plot_surface(xs, ys, zs, cmap=plt.cm.YlGnBu, alpha = 0.1)



plt.show()

#end for now
sys.exit()


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
list_safe = []

for f in filename:
    df = pd.DataFrame()
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
    list_safe.append(df)

#place each source position in a bin

def order_in_bins(df):
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
    return charge_uniform_bins_theta, count_uniform_bins_theta

charge_uniform_bins_theta, count_uniform_bins_theta = order_in_bins(df)

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

if len(list_safe) == 2:
    #plottong a binned comparision
    print('Now the comparision')
    charge_uniform_bins_theta1, count_uniform_bins_theta1 = order_in_bins(list_safe[0])
    charge_uniform_bins_theta2, count_uniform_bins_theta2 = order_in_bins(list_safe[1])
    charge_uniform_bins_theta_diff_tot = charge_uniform_bins_theta1/count_uniform_bins_theta1 - charge_uniform_bins_theta2/count_uniform_bins_theta2

    ax2 = plt.axes(projection = 'polar')
    ax2.set_thetamin(0)
    ax2.set_ylabel("Difference \n %s "%(filename))
    ax2.set_thetamax(phi_max)
    print(filename[0], filename[1], "difference")
    for b in range(len(bins_all[0])):
        if count_uniform_bins_theta1[b] * count_uniform_bins_theta2[b]!=0: #npt empty bins
        #/( count_uniform_bins_theta[b]* nEvents[0]
            charge_uniform_bins_theta_diff = charge_uniform_bins_theta1[b]/(count_uniform_bins_theta1[b]*1000)-charge_uniform_bins_theta2[b]/(count_uniform_bins_theta2[b]*1000)
            print(" File 800: ", charge_uniform_bins_theta1[b]/(count_uniform_bins_theta1[b]*1000), " file 805: ", charge_uniform_bins_theta2[b]/(count_uniform_bins_theta2[b]*1000), " difference : ", charge_uniform_bins_theta_diff)


            sc2 = ax2.scatter(bins_all[4][b], bins_all[3][b], c = charge_uniform_bins_theta_diff, cmap = "Blues",vmin = -0.1, vmax = 0.1, s = 120)
    plt.colorbar(sc2)
    ax2.set_title('Difference in number of pe per photon')
    #plt.savefig('/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/Maps/maps_pictures/Binned_visualisation/Non-interpolated/%s.png'%outputfile_name)
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
ax.scatter3D(PMT_x,PMT_z,PMT_y, marker = 'o')
for s in range(int(len(PMT_x))):
    ax.plot([PMT_x[s], PMT_dirx[s]+PMT_x[s]],[PMT_z[s], PMT_dirz[s]+PMT_z[s]],[PMT_y[s], PMT_diry[s]+PMT_y[s]], 'k-')

#This is the position of the centre of the uniform bins
for b in range(len(bins_all)):
    ax.scatter3D(bins_all[0][b],bins_all[2][b], bins_all[1][b], 'o', color = 'k')

for i in range(len(x)):
    ax.plot([x[i], x[i] - R[i] *  np.sin(theta[i]) * np.cos(phi[i])], [z[i], z[i] - R[i] * np.sin(theta[i]) * np.sin(phi[i])], [y[i], y[i] - R[i] * np.cos(theta[i])], 'r-')
plt.show()






