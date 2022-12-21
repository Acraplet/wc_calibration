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

interpolate = True

source_xyz, source_Rtp, source_faction_evDigiHits, source_faction_evRawHits, source_fraction_rawQ = [],[],[],[],[]
source_NHits, source_rawNHits, source_fraction_Q, source_mean_evRawQ, source_mean_evDigiQ = [],[],[],[],[]
print(len(table))

for ev in range(len(table)):
    t = table[ev]
    source_xyz.append([t[0], t[1], t[2]])
    source_Rtp.append([t[3], t[4], t[5]])
    source_faction_evDigiHits.append(t[6])

    source_faction_evRawHits.append(t[7])
    source_fraction_rawQ.append(t[8])
    source_fraction_Q.append(t[9])
    source_NHits.append(t[10])
    source_rawNHits.append(t[11])
    source_mean_evRawQ.append(t[12])
    source_mean_evDigiQ.append(t[13])

source_xyz = np.array(source_xyz)
source_Rtp = np.array(source_Rtp)

    #this is what was saved
    #np.array(source_xyz).T[:][0], np.array(source_xyz).T[:][1], np.array(source_xyz).T[:][2], np.array(source_Rtp).T[:][0], np.array(source_Rtp).T[:][1], np.array(source_Rtp).T[:][2], source_faction_evDigiHits, source_faction_evRawHits, source_fraction_rawQ, source_fraction_Q, source_NHits, source_rawNHits, source_mean_evRawQ, source_mean_evDigiQ
#Here are the plotting functions from before:

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
#print(np.array(source_xyz).T[0])
R = source_Rtp.T[0][1]
print(R)
theta = source_Rtp.T[1]
phi = source_Rtp.T[2]

x, y, z = source_xyz.T[0], source_xyz.T[1], source_xyz.T[2]

ax.scatter3D(0, 0, - 128.05 - 27.4, marker = 'x', color = 'k')
ax.scatter3D(x, z, y, marker = 'x', color ='r')
for i in range(len(x)):
    ax.plot([x[i], x[i] - R *  np.sin(theta[i]) * np.cos(phi[i]) ], [z[i], z[i] - R * np.sin(theta[i]) * np.sin(phi[i])], [y[i], y[i] - R * np.cos(theta[i])], 'r-')
plt.show()


#Plot the fraction of all the detector hits that fall onto the 58th mPMT
plt.subplot(2,2,1)
plt.suptitle("%s"%name)
im = plt.scatter(np.array(source_Rtp).T[:][2],np.array(source_Rtp).T[:][1], c = source_fraction_rawQ)
col = plt.colorbar(im, label="Fraction of digi Q in mPMT 58", orientation="vertical", format= "%.1e")
im.figure.axes[1].tick_params(axis="y", labelsize=12)

plt.xlabel(f'$\phi(rad)$')
plt.ylabel(f'$\Theta$(rad)')
plt.subplot(2,2,2)
im2 = plt.scatter(np.array(source_Rtp).T[:][2],np.array(source_Rtp).T[:][1], c = source_fraction_Q)
col2 = plt.colorbar(im2, label="Fraction of raw Q in mPMT 58", orientation="vertical", format= "%.1e")
im2.figure.axes[3].tick_params(axis="y", labelsize=12)

plt.xlabel(f'$\phi(rad)$')
plt.ylabel(f'$\Theta$(rad)')
#plt.subplot(2,2,3)
#im3 = plt.scatter(np.array(source_Rtp).T[:][2],np.array(source_Rtp).T[:][1], c = source_fraction_NHits)
#col3 = plt.colorbar(im3, label="Fraction of digi hits in mPMT 58", orientation="vertical", format= "%.1e")
#im3.figure.axes[5].tick_params(axis="y", labelsize=12)

#plt.xlabel(f'$\phi(rad)$')
#plt.ylabel(f'$\Theta$(rad)')
#plt.subplot(2,2,4)
#im4 = plt.scatter(np.array(source_Rtp).T[:][2],np.array(source_Rtp).T[:][1], c = source_fraction_rawNHits)
#col4 = plt.colorbar(im4, label="Fraction of raw hits in mPMT 58", orientation="vertical", format= "%.1e")
#im4.figure.axes[7].tick_params(axis="y", labelsize=12)

#plt.xlabel(f'$\phi(rad)$')
#plt.ylabel(f'$\Theta$(rad)')
plt.show()


#Plot the fraction of events out of which the 58th mPMT sees at least one hit
levels = np.linspace(0,1, 10)

plt.subplot(1,2,1)
plt.suptitle("%s"%name)
plt.title("DigiNHits")

if (interpolate):
    x = np.array(source_Rtp).T[:][2]
    y = np.array(source_Rtp).T[:][1]
    z = source_faction_evDigiHits
    xi, yi = np.linspace(min(x), max(x), 100), np.linspace(min(y), max(y), 100)
    xi, yi = np.meshgrid(xi, yi)
    # Interpolate
    rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
    zi = rbf(xi, yi)
    plt.contourf(xi, yi, zi, 500, cmap='jet')




im = plt.scatter(np.array(source_Rtp).T[:][2],np.array(source_Rtp).T[:][1], c = source_faction_evDigiHits, cmap = "nipy_spectral", s = 40)
col = plt.colorbar(im, label="Fraction of events for which there is at least 1 *digi* hit in mPMT58", orientation="horizontal", format= "%.1f")
col.set_ticks(levels)
im.figure.axes[1].tick_params(axis="x", labelsize=12)
plt.clim(0,1)
plt.xlabel(f'$\phi(rad)$')
plt.ylabel(f'$\Theta$(rad)')



plt.subplot(1,2,2)
plt.title("NHits")
im = plt.scatter(np.array(source_Rtp).T[:][2],np.array(source_Rtp).T[:][1], c = source_faction_evRawHits, cmap = "nipy_spectral",s = 40)
col2 = plt.colorbar(im, label=f"Fraction of events for which there is at least 1 *raw* hit in mPMT58", orientation="horizontal", format= "%.1f")
col2.set_ticks(levels)
plt.clim(0,1)
im.figure.axes[3].tick_params(axis="x", labelsize=12)

plt.xlabel(f'$\phi(rad)$')
plt.ylabel(f'$\Theta$(rad)')

plt.show()




#Now interpolate Hits
levels = np.linspace(0,1, 10)

plt.subplot(1,2,1)
plt.suptitle("%s"%name)
plt.title("DigiNHits")

x = np.array(source_Rtp).T[:][2]
y = np.array(source_Rtp).T[:][1]
z = source_faction_evDigiHits
xi, yi = np.linspace(min(x), max(x), 100), np.linspace(min(y), max(y), 100)
xi, yi = np.meshgrid(xi, yi)
# Interpolate
rbf = scipy.interpolate.Rbf(x, y, z, function='cubic')
zi = rbf(xi, yi)
plt.contourf(xi, yi, zi, 500, cmap='nipy_spectral')
plt.clim(0,1)


im = plt.scatter(np.array(source_Rtp).T[:][2],np.array(source_Rtp).T[:][1], c = source_faction_evDigiHits, cmap = "nipy_spectral", s = 40)
col = plt.colorbar(im, label="Fraction of events for which there is at least 1 *digi* hit in mPMT58", orientation="horizontal", format= "%.1f")
col.set_ticks(levels)
im.figure.axes[1].tick_params(axis="x", labelsize=12)
plt.clim(0,1)
plt.xlabel(f'$\phi(rad)$')
plt.ylabel(f'$\Theta$(rad)')



plt.subplot(1,2,2)

x = np.array(source_Rtp).T[:][2]
y = np.array(source_Rtp).T[:][1]
z = source_faction_evRawHits
xi, yi = np.linspace(min(x), max(x), 100), np.linspace(min(y), max(y), 100)
xi, yi = np.meshgrid(xi, yi)
# Interpolate
rbf = scipy.interpolate.Rbf(x, y, z, function='cubic')
zi = rbf(xi, yi)
plt.contourf(xi, yi, zi, 500, cmap='nipy_spectral')
plt.clim(0,1)

plt.title("NHits")
im = plt.scatter(np.array(source_Rtp).T[:][2],np.array(source_Rtp).T[:][1], c = source_faction_evRawHits, cmap = "nipy_spectral",s = 40)
col2 = plt.colorbar(im, label=f"Fraction of events for which there is at least 1 *raw* hit in mPMT58", orientation="horizontal", format= "%.1f")
col2.set_ticks(levels)
plt.clim(0,1)
im.figure.axes[3].tick_params(axis="x", labelsize=12)

plt.xlabel(f'$\phi(rad)$')
plt.ylabel(f'$\Theta$(rad)')

plt.show()

#Now interpolate Q
levels = np.linspace(0,1, 10)

plt.subplot(1,2,1)
plt.suptitle("%s"%name)
plt.title("DigiQ")

x = np.array(source_Rtp).T[:][2]
y = np.array(source_Rtp).T[:][1]
z = source_mean_evDigiQ
xi, yi = np.linspace(min(x), max(x), 100), np.linspace(min(y), max(y), 100)
xi, yi = np.meshgrid(xi, yi)
# Interpolate
if (interpolate):
    rbf = scipy.interpolate.Rbf(x, y, z, function='cubic')
    zi = rbf(xi, yi)
    plt.contourf(xi, yi, zi, 500, cmap='nipy_spectral')
#plt.clim(0,1)


im = plt.scatter(np.array(source_Rtp).T[:][2],np.array(source_Rtp).T[:][1], c = source_mean_evDigiQ, cmap = "nipy_spectral", s = 40)
col = plt.colorbar(im, label="Average number of *digitised* charge in mPMT58", orientation="horizontal", format= "%.1f")
#col.set_ticks(levels)
im.figure.axes[1].tick_params(axis="x", labelsize=12)
#plt.clim(0,1)
plt.xlabel(f'$\phi(rad)$')
plt.ylabel(f'$\Theta$(rad)')



plt.subplot(1,2,2)

x = np.array(source_Rtp).T[:][2]
y = np.array(source_Rtp).T[:][1]
z = source_mean_evRawQ
xi, yi = np.linspace(min(x), max(x), 100), np.linspace(min(y), max(y), 100)
xi, yi = np.meshgrid(xi, yi)
# Interpolate
if (interpolate):
    rbf = scipy.interpolate.Rbf(x, y, z, function='cubic')
    zi = rbf(xi, yi)
    plt.contourf(xi, yi, zi, 500, cmap='nipy_spectral')
#plt.clim(0,1)

plt.title("Q")
im = plt.scatter(np.array(source_Rtp).T[:][2],np.array(source_Rtp).T[:][1], c = source_mean_evRawQ, cmap = "nipy_spectral",s = 40)
col2 = plt.colorbar(im, label=f"Average number of *raw* charge in mPMT58", orientation="horizontal", format= "%.1f")
#col2.set_ticks(levels)
#plt.clim(0,1)
im.figure.axes[3].tick_params(axis="x", labelsize=12)

plt.xlabel(f'$\phi(rad)$')
plt.ylabel(f'$\Theta$(rad)')

plt.show()
