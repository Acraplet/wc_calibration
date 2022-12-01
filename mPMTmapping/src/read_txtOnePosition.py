#This is a python file to read the maps in txt format and do the plotting faster
import read_data_results3 as rd
#read the .txt file:
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import scipy.optimize
import getopt
import pandas as pd
import iminuit as im
import scipy.interpolate as spi


filename = "No filename, please input one with the -f option"
comparision_file = 'none'
outputfile_name = 'none'

plt.style.use(["science", "notebook", "grid"])

#filename = str(sys.argv[1])
phi_max = 90 #for plotting - get the max phi value to show only the right portion of the circle
argv = sys.argv[1:]

PMT_position_file = '/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/PMT_positions.txt'
PMT = np.array(rd.read_data3(PMT_position_file)).T
PMT_mPMT = PMT[0]
PMT_mPMT_pmt = PMT[1]
PMT_x = PMT[2]
PMT_y = PMT[3]
PMT_z = PMT[4]

filename_ref = '/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_txtFiles/mPMT_map_ID116.txt'

def exponential(alpha_abs, A, R):
    #print(u)
    #print(alpha_abs)
    return A  * np.exp(-R / alpha_abs)

def least_squares_np(params): #a= amplitude, o=offset, p=phase, l=wavelength xx = portion to fit in x and yy in y
  sigma=0.01
  return sum((yy-exponential(xx,*params))**2/sigma)



def dist(df_geom_buf, df_buf):
    #df_geom_buf = df_geom[df_geom['mPMT_pmt'] == PMT]
    #print(df_geom_buf['x'])
    R = df_buf['R']
    theta = df_buf['theta']
    phi = df_buf['phi']
    x = df_buf['x'] - R *  np.sin(theta)* np.cos(phi)
    z = df_buf['z'] - R * np.sin(theta) * np.sin(phi)
    y = df_buf['y'] - R * np.cos(theta)
    return (x - df_geom_buf['x'])**2 + (y - df_geom_buf['y'])**2 + (z - df_geom_buf['z'])**2

df_geom = pd.DataFrame()
for i in range(len(PMT[0])):
        c = [PMT[0][i],PMT[1][i],PMT[2][i],PMT[3][i],PMT[4][i]]
        row =  pd.Series(data=c, index=['mPMT', 'mPMT_pmt', 'x', 'y', 'z'], dtype=np.float64)
        df_geom = df_geom.append(row, ignore_index=True)


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


theta_target = float(filename.split("_theta")[1].split("_")[0])
phi_target = float(filename.split("_phi")[1].split("_")[0])
R_target = float(filename.split("_R")[1].split(".")[0])

color = ['black', 'red', 'blue', 'green', 'darkorange', 'cyan', 'magenta']

a = 0
plt.figure(0, figsize =(25,20) )
base = filename.split("_R")[0]
outname = 'OnePosition_theta%.2f_phi%.2f_allR'%(theta_target, phi_target)
R_list = [5.00, 10.00, 20.00, 40, 80, 160, 320]
alpha_list = [10, 20, 40, 60, 100,150,220]

table_ref = rd.read_data3(filename_ref )
source_xyz_ref  = [table_ref[0], table_ref[1], table_ref[2]]
source_Rtp_ref  = [table_ref[3], table_ref[4], table_ref[5]]
x_ref , y_ref , z_ref  = np.array(table_ref [0]), np.array(table_ref [1]), np.array(table_ref [2])
R_ref , theta_ref , phi_ref  = np.array(table_ref [5]), np.array(table_ref [3]), np.array(table_ref [4])
Q_tot_ref  =  np.array(table_ref[6])
nEvents_ref  = np.array(table_ref[7])

df_ref  = pd.DataFrame()
for i in range(len(x_ref)):
    c = [x_ref [i], y_ref [i],  z_ref [i],  theta_ref [i], phi_ref [i],  R_ref [i],  Q_tot_ref [i],  nEvents_ref [i]]
    row =  pd.Series(data=c, index=['x', 'y', 'z', 'theta', 'phi', 'R', 'Q', 'events'], dtype=np.float64)
    df_ref  = df_ref.append(row, ignore_index=True)

df_ref = df_ref[df_ref['theta']==theta_target]
df_ref = df_ref[df_ref['phi'] == phi_target]

print( df_ref[df_ref['phi'] == phi_target])


for j in R_list:
    print(j)
    filename = "%s_R%.2f.txt"%(base, float(j))
    table = rd.read_data3(filename)
    source_xyz = [table [0], table[1], table [2]]
    source_Rtp = [table [3], table [4], table[5]]
    x, y, z = np.array(table[0]), np.array(table[1]), np.array(table[2])
    R, theta, phi = np.array(table[5]), np.array(table[3]), np.array(table[4])
    Q_tot =  np.array(table[6])
    nEvents = np.array(table[7])
    abwff =  np.array(table[8])
    rayff = np.array(table[9])

    phi_max = max(phi) * 180/np.pi
    df = pd.DataFrame()
    for i in range(len(x)):
        c = [x[i], y[i],  z[i],  theta[i], phi[i],  R[i],  Q_tot[i],  nEvents[i], abwff[i], rayff[i]]
        row =  pd.Series(data=c, index=['x', 'y', 'z', 'theta', 'phi', 'R', 'Q', 'events', 'abwff', 'rayff'], dtype=np.float64)
        df = df.append(row, ignore_index=True)

    df_abwff = df[df['rayff'] >= 10000]
    df_rayff = df[df['abwff'] >= 10000]
    df_abwff = df_abwff.sort_values(["abwff"], axis = 0, ascending = True).reset_index()
    df_rayff = df_rayff.sort_values(["rayff"], axis = 0, ascending = True).reset_index()

    xs = np.linspace(10, 220, 1000)

    if float(df_ref['Q']) == 0:
        df_ref['Q'] = 0.00001
    #plt.figure(1)
    #plt.subplot(2,1,1)
    #plt.plot(df_abwff['abwff'], alpha_list, 'x')
    #res = scipy.optimize.curve_fit(exponential, df_abwff['abwff'],  alpha_list)
    #plt.plot(xs,exponential(xs, *res[0]), label="R = %icm"%j,color = color[a])
    #plt.xlabel('abwff')
    #plt.ylabel('alpha (attenuation length)')
    #plt.subplot(2,1,2)
    #xs = np.linspace(0, 0.013, 1000)
    #plt.plot(df_rayff['rayff'], alpha_list, 'x')
    #res = scipy.optimize.curve_fit(exponential, df_rayff['rayff'],  alpha_list)
    #plt.plot(xs,exponential(xs, *res[0]), label="R = %icm"%j,color = color[a])
    #plt.show()


    plt.figure(0)
    #df_abwff['abwff']
    plt.subplot(2,1,1)
    plt.title('Charge as a function of scattering and abosption - spline interpolation \n parameters at position: theta= %.2frad phi= %.2frad R= %icm'%(theta_target, phi_target, R_target))
    plt.errorbar(np.array(alpha_list), df_abwff['Q']/float(df_ref['Q']), yerr = np.sqrt(df_abwff['Q'] * (1-df_abwff['Q']/df_abwff['events']))/float(df_ref['Q']), fmt = 'o', color = color[a], markersize = 10)


    y = df_abwff['Q']/float(df_ref['Q'])
    y_err = np.sqrt(df_abwff['Q'] * (1-df_abwff['Q']/df_abwff['events']))/float(df_ref['Q'])
    print(y)
    #cs = spi.CubicSpline(alpha_list, y )

    #print(df_abwff['abwff'], df_abwff['Q']/df_ref['Q'])
    #res = scipy.optimize.curve_fit(exponential, df_abwff['abwff'], df_abwff['Q']/df_ref['Q'])
    xx = np.array(alpha_list) #df_abwff['abwff']
    yy =  df_abwff['Q']/float(df_ref['Q'])
    po = [0.3, j]
    m=im.Minuit(least_squares_np, (po[0], po[1]))
    m.migrad()
    res=m.values
    print(res)

    #plt.plot(xs,cs(xs), label="R = %icm"%j,color = color[a], alpha = 0.5)
    #plt.yscale('log')
    plt.plot(xs,exponential(xs, *res),label="R = %icm"%j,color = color[a], alpha = 0.5)



    plt.xlabel(f'1/Alpha Absorption')
    plt.ylabel('Charge collected\nin fraction of the max charge')
    plt.subplot(2,1,2)
    plt.errorbar(df_rayff['rayff'], df_rayff['Q']/float(df_ref['Q']), yerr = np.sqrt(df_rayff['Q'] * (1-df_rayff['Q']/df_rayff['events']))/float(df_ref['Q']), fmt = 'o', color = color[a], markersize = 10)
    plt.xlabel(f'Rayleigh scattering parameter')
    plt.ylabel('Charge collected\nin fraction of the max charge')



    xs = np.linspace(0, 0.013, 1000)
    cs = spi.CubicSpline(df_rayff['rayff'], df_rayff['Q']/float(df_ref['Q']))
    plt.plot(xs,cs(xs), label="R = %icm"%j,color = color[a])



    a = a+1
    #plt.yscale('log')






plt.subplot(2,1,2)
plt.legend()
plt.subplot(2,1,1)
plt.legend()
plt.savefig('/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_pictures/Maps/OnePosition/%s.png'%outname)
#plt.show()
#plt.show()

raise End



























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
    plt.style.use(["science", "notebook", "grid"])
    distance_to_closest_PMT = []
    Q_array = []
    #check the closet distance to PMT
    for index, row in df.iterrows():
        dist_min = 100000
        for index_geom, row_geom in df_geom.iterrows():
        #for PMT in df_geom['mPMT_pmt'].unique():
           #print()
           distance = dist(row_geom, row)
           if distance <= dist_min:
               dist_min = distance
        distance_to_closest_PMT.append(dist_min)
        Q_array.append(row['Q']/row['events'])
    plt.title('Charge collected vs distance to closest PMT\nFileID%s'%filename)
    plt.xlabel('Distance to the closest PMT (cm)')
    plt.ylabel('Charge collected per photon')
    plt.plot(distance_to_closest_PMT, Q_array, 'x')
    plt.xscale('log')
    plt.show()

    raise End
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
ax.scatter3D(PMT_x,PMT_z,PMT_y, marker = 'o')

for i in range(len(x)):
    ax.plot([x[i], x[i] - R[i] *  np.sin(theta[i]) * np.cos(phi[i])], [z[i], z[i] - R[i] * np.sin(theta[i]) * np.sin(phi[i])], [y[i], y[i] - R[i] * np.cos(theta[i])], 'r-')
plt.show()





