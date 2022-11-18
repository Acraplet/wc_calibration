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
argv = sys.argv[1:]

interpolation_mode = 'linear' # 'linear'

opts, args = getopt.getopt(argv, "m:o:p:f:c:a:")
for opt, arg in opts:
        if opt in ['-m']:
            phi_max = float(arg)
        elif opt in ['-f']:
            filename = str(arg)
        elif opt in ['-c']:
            comparision_file = str(arg)
        elif opt in ['-o']:
            outputfile_name = str(arg)


if outputfile_name == 'none':
    print("No outputfile_name given - naming the output 'output.png'")
    outputfile_name = 'output'


if comparision_file == 'none':
    table = rd.read_data3(filename)
    source_xyz = [table [0], table[1], table [2]]
    source_Rtp = [table [3], table [4], table[5]]
    x, y, z = np.array(table[0]), np.array(table[1]), np.array(table[2])
    R, theta, phi = np.array(table[5]), np.array(table[3]), np.array(table[4])
    Q_tot =  np.array(table[6])
    nEvents = np.array(table[7])
    phi_max = max(phi) * 180/np.pi

    #print(phi_max)

if comparision_file != 'none':
    df = pd.DataFrame()
    df_compa = pd.DataFrame()

    table = rd.read_data3(filename)
    table_compa = rd.read_data3(comparision_file)
    source_xyz = [table [0], table[1], table [2]]
    source_Rtp = [table [3], table [4], table[5]]

    source_xyz_compa = [table_compa [0], table_compa[1], table_compa[2]]
    source_Rtp_compa = [table_compa[3], table_compa[4], table_compa[5]]

    x, y, z = np.array(table[0]), np.array(table[1]), np.array(table[2])
    R, theta, phi = np.array(table[5]), np.array(table[3]), np.array(table[4])

    x_compa, y_compa, z_compa = np.array(table_compa[0]), np.array(table_compa[1]), np.array(table_compa[2])
    R_compa, theta_compa, phi_compa = np.array(table_compa[5]), np.array(table_compa[3]), np.array(table_compa[4])

    phi_max = max(phi) * 180/np.pi
    phi_max_compa = max(phi_compa) * 180/np.pi

    Q_tot_compa = np.array(table_compa[6])
    Q_tot = np.array(table[6])
    nEvents_compa = np.array(table_compa[7])
    nEvents = np.array(table[7])

    for i in range(len(x)):
        c = [x[i], y[i],  z[i],  theta[i], phi[i],  R[i],  Q_tot[i],  nEvents[i]]
        row =  pd.Series(data=c, index=['x', 'y', 'z', 'theta', 'phi', 'R', 'Q', 'events'], dtype=np.float64)
        df = df.append(row, ignore_index=True)

    for j in range(len(x_compa)):
        c = [x_compa[j],  y_compa[j], z_compa[j], theta_compa[j], phi_compa[j], R_compa[j], Q_tot_compa[j],  nEvents_compa[j]]
        row =  pd.Series(data=c, index=['x', 'y', 'z', 'theta', 'phi', 'R', 'Q', 'events'], dtype=np.float64, name = 'compa')
        df_compa = df_compa.append(row, ignore_index=True)



    if phi_max != phi_max_compa:
        print("\n\nWe are looking at different file format, need to cut to compatibility !! \n \n")

        print(df["phi"].max(), df_compa["phi"].max())

        df_compa=df_compa[df_compa["phi"]<=df["phi"].max()].reset_index()
        df=df[df["phi"]<=df_compa["phi"].max()].reset_index()

        print(df_compa, "\n", df)





if comparision_file == 'none':
    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(projection='polar')
    #only the points
    ax.set_title("Q")
    #plt.subplot(1,1,1)
    zr = Q_tot/nEvents
    im = ax.scatter(phi, theta, c = zr, cmap = "nipy_spectral",s = 40)
    col2 = plt.colorbar(im, label=f"Average number of *raw* P.E. in mPMT58 per photon", orientation="horizontal", format= "%.2f")

    ax.set_thetamin(0)
    ax.set_ylabel("Q %s "%(filename))
    ax.set_thetamax(phi_max)
    ##ax.set_phimin(0)
    ##ax.set_phimax(90)
    #im.figure.axes[1].tick_params(axis="x", labelsize=12)
    #ax.set_ylabel(f'$\phi()$')
    ax.set_xlabel(f'$\Theta$(rad)')
    ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 20
    plt.savefig('/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_pictures/Maps/Non-interpolated/%s.png'%outputfile_name)
    plt.show()

    #Now the interpolate
    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(projection='polar')
    ax.set_title("Q")
    ax.set_ylabel("Q %s "%(filename))
    zr = Q_tot/nEvents
    xi, yi = np.linspace(min(phi), max(phi), 100), np.linspace(min(theta), max(theta), 100)
    xi, yi = np.meshgrid(xi, yi)
    # Interpolate
    rbf = scipy.interpolate.Rbf(phi, theta, zr, function=interpolation_mode)
    zi = rbf(xi, yi)
    ax.contourf(xi, yi, zi, 500, cmap='nipy_spectral')

    im = ax.scatter(phi, theta, c = zr, cmap = "nipy_spectral",s = 40)
    col2 = plt.colorbar(im, label=f"Average number of *raw* P.E. in mPMT58 per photon", orientation="horizontal", format= "%.2f")

    ax.set_thetamin(0)
    ax.set_thetamax(phi_max)
    ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 20
    ax.set_xlabel(f'$\Theta$(rad)')
    plt.savefig('/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_pictures/Maps/Interpolated/%s.png'%outputfile_name)
    plt.show()


if comparision_file != 'none':



    fig = plt.figure(figsize=(20,10))

############## non-interpolated first plot ################
    ax = fig.add_subplot(221, projection='polar')
    ax.set_ylabel("Q %s "%(filename))
    zr = df["Q"]/df["events"]
    zr2 = df_compa["Q"]/df_compa["events"] #needed to have a common colorbar
    vmin = min(min(zr2), min(zr))
    vmax = max(max(zr2), max(zr))
    im = ax.scatter(df["phi"], df["theta"], c = zr, cmap = "nipy_spectral",s = 40,vmin = vmin, vmax = vmax)
    col2 = plt.colorbar(im, label=f"Average number of *raw* P.E. in mPMT58 per photon", orientation="vertical", format= "%.2f")
    ax.set_thetamax(phi_max)
    ax.xaxis.labelpad = 10
    ax.set_xlabel(f'$\Theta$(rad)')
############## non-interpolated second plot ################
    ax2 = fig.add_subplot(223, projection='polar')
    ax2.set_ylabel("Q %s"%(comparision_file))
    zr2 = df_compa["Q"]/df_compa["events"]
    im2 = ax2.scatter(df_compa["phi"], df_compa["theta"], c = zr2, cmap = "nipy_spectral",s = 40,vmin = vmin, vmax = vmax)
    col3 = plt.colorbar(im2, label=f"Average number of *raw* P.E. in mPMT58 per photon", orientation="vertical", format= "%.2f")
    #col3.set_clim(min(min(zr2), min(zr)), max(max(zr2), max(zr)))
    ax2.set_thetamin(0)
    ax2.set_thetamax(phi_max)
    ax2.xaxis.labelpad = 10
    ax2.set_xlabel(f'$\Theta$(rad)')
######need to re-order the dataframes so that we can make the comparision, this is for two files with the same id of the datapoints
    df = df.sort_values(["theta", "phi"], axis = 0, ascending = False).reset_index()
    df_compa = df_compa.sort_values(["theta", "phi"], axis = 0, ascending = False).reset_index()
    df_diff = df-df_compa
    if sum(df_diff['phi'])!=0:
        print(df_diff['phi'])
        print('It looks like the datasets do not have their points at the same position, please check what is up with this')
        raise MisleadingComparision

    print(df["theta"].unique())

############## non-interpolated third plot ################
    ax2 = fig.add_subplot(122, projection='polar')
    ax2.set_title("Q Comparision \n %s-%s"%(filename, comparision_file))

    zr2 = (df['Q']/df['events'] - df_compa['Q']/df_compa['events']) #/df['Q']/df['events']
    im2 = ax2.scatter(df['phi'], df['theta'], c = zr2, cmap = "nipy_spectral",s = 40)
    col3 = plt.colorbar(im2, label=f"Difference in number of *raw* P.E. in mPMT58 per photon", orientation="vertical", format= "%.2f")

    ax2.set_thetamin(0)
    ax2.set_thetamax(phi_max)
    ax2.xaxis.labelpad = 10

    ax2.set_xlabel(f'$\Theta$(rad)')
    plt.savefig('/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_pictures/Comparisions/Non-interpolated/%s.png'%outputfile_name)
    plt.show()

#### NOW with the interpolation

############## interpolated first plot ################
    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(221, projection='polar')
    ax.set_ylabel("Q %s "%(filename))
    zr =  df["Q"]/df["events"]

    xi, yi = np.linspace(min(phi), max(phi), 200), np.linspace(min(theta), max(theta), 200)
    xi, yi = np.meshgrid(xi, yi)
    # Interpolate
    rbf = scipy.interpolate.Rbf(df["phi"], df["theta"], zr, function=interpolation_mode, vmin = vmin, vmax = vmax)
    zi = rbf(xi, yi)
    ax.contourf(xi, yi, zi, 500, cmap='nipy_spectral',vmin = vmin, vmax = vmax)

    im = ax.scatter(df["phi"], df["theta"], c = zr, cmap = "nipy_spectral",s = 40, vmin = vmin, vmax = vmax)
    col2 = plt.colorbar(im, label=f"Average number of *raw* P.E. in mPMT58 per photon", orientation="vertical", format= "%.2f")
    ax.set_thetamin(0)
    ax.set_thetamax(phi_max)
    ax.xaxis.labelpad = 10
    ax.set_xlabel(f'$\Theta$(rad)')

############## interpolated second plot ################
    ax2 = fig.add_subplot(223, projection='polar')
    ax2.set_ylabel("Q %s"%(comparision_file))
    zr2 = df_compa["Q"]/df_compa["events"]

    xi, yi = np.linspace(min(phi), max(phi), 200), np.linspace(min(theta), max(theta), 200)
    xi, yi = np.meshgrid(xi, yi)
    # Interpolate
    rbf = scipy.interpolate.Rbf(df_compa["phi"], df_compa["theta"], zr2, function=interpolation_mode, vmin = vmin, vmax = vmax)
    zi = rbf(xi, yi)
    ax2.contourf(xi, yi, zi, 500, cmap='nipy_spectral',vmin = vmin, vmax = vmax)

    im2 = ax2.scatter(df_compa["phi"], df_compa["theta"], c = zr2, cmap = "nipy_spectral",s = 40, vmin = vmin, vmax = vmax)
    col3 = plt.colorbar(im2, label=f"Average number of *raw* P.E. in mPMT58 per photon", orientation="vertical", format= "%.2f")

    ax2.set_thetamin(0)
    ax2.set_thetamax(phi_max)
    ax2.xaxis.labelpad = 10
    ax2.set_xlabel(f'$\Theta$(rad)')

############## interpolated third plot ################
    ax3 = fig.add_subplot(122, projection='polar')
    ax3.set_title("Q Comparision \n %s-%s"%(filename, comparision_file))

    zr2 = df['Q']/df['events'] - df_compa['Q']/df_compa['events'] ##the difference in number of p.e. per event
    #zr2 = (df['Q']/df['events'] - df_compa['Q']/df_compa['events'])/df['Q']/df['events']
    xi, yi = np.linspace(min(df['phi']), max(df['phi']), 200), np.linspace(min(df['theta']), max(df['theta']), 200)
    xi, yi = np.meshgrid(xi, yi)
    # Interpolate
    rbf = scipy.interpolate.Rbf(df['phi'], df['theta'], zr2, function=interpolation_mode)
    zi = rbf(xi, yi)
    ax3.contourf(xi, yi, zi, 500, cmap='nipy_spectral')
    im2 = ax3.scatter(df['phi'], df['theta'], c = zr2, cmap = "nipy_spectral",s = 40)
    col3 = plt.colorbar(im2, label=f"Difference in number of *raw* P.E. in mPMT58 per photon", orientation="vertical", format= "%.2f")

    ax3.set_thetamin(0)
    ax3.set_thetamax(phi_max)
    ax3.xaxis.labelpad = 10

    ax3.set_xlabel(f'$\Theta$(rad)')
    plt.savefig('/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_pictures/Comparisions/Interpolated/%s.png'%outputfile_name)
    plt.show()

source_Rtp = np.array(source_Rtp)
source_xyz = np.array(source_xyz)

#the figure
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
#x.plot()

ax.scatter3D(0, 0, - 128.05 - 27.4, marker = 'x', color = 'k')
ax.scatter3D(x, z, y, marker = 'x', color ='r')


for i in range(len(x)):
    ax.plot([x[i], x[i] - R[i] *  np.sin(theta[i]) * np.cos(phi[i])], [z[i], z[i] - R[i] * np.sin(theta[i]) * np.sin(phi[i])], [y[i], y[i] - R[i] * np.cos(theta[i])], 'r-')
plt.show()





