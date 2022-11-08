#This is a python file to read the maps in txt format and do the plotting faster
import read_data_results3 as rd
#read the .txt file:
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

filename = str(sys.argv[1])
name = filename
table = rd.read_data3(filename)
source_xyz = [table [0], table[1], table [2]]
source_Rtp = [table [3], table [4], table[5]]


x, y, z = np.array(table[0]), np.array(table[1]), np.array(table[2])
R, theta, phi = np.array(table[5]), np.array(table[3]), np.array(table[4])

Q_tot =  np.array(table[6])
nEvents = np.array(table[7])

#print(table[0])

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
R = R
print(R)
for i in range(len(x)):
    ax.plot([x[i], x[i] - R[i] *  np.sin(theta[i]) * np.cos(phi[i])], [z[i], z[i] - R[i] * np.sin(theta[i]) * np.sin(phi[i])], [y[i], y[i] - R[i] * np.cos(theta[i])], 'r-')
plt.show()


fig = plt.figure()
ax = fig.add_subplot(projection='polar')
ax.set_title("Q")
##plt.subplot(1,1,1)
zr = Q_tot/nEvents
xi, yi = np.linspace(min(phi), max(phi), 100), np.linspace(min(theta), max(theta), 100)
xi, yi = np.meshgrid(xi, yi)
# Interpolate
rbf = scipy.interpolate.Rbf(phi, theta, zr, function='cubic')
zi = rbf(xi, yi)
ax.contourf(xi, yi, zi, 500, cmap='nipy_spectral')
#ax.clim(0,1)

im = ax.scatter(phi, theta, c = zr, cmap = "nipy_spectral",s = 40)
col2 = plt.colorbar(im, label=f"Average number of *raw* P.E. in mPMT58 per photon", orientation="horizontal", format= "%.1f")
#col2.set_ticks(levels)
#ax.clim(0,1)

ax.set_thetamin(0)
ax.set_thetamax(90)
##ax.set_phimin(0)
##ax.set_phimax(90)
#im.figure.axes[1].tick_params(axis="x", labelsize=12)
#ax.set_ylabel(f'$\phi()$')
ax.set_xlabel(f'$\Theta$(rad)')
#plt.savefig("../maps_pictures/mPMT58_Efficiency_Map_FileID%i.png"%run_id)
plt.show()

