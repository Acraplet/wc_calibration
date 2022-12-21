#Same analysis of the WCSim files but in python - much easier!
#to use from the data folder python ../src/make_raw_mPMTmap.py $(ls *ID6*)
#here, only to look at the simple case of the ras data

import uproot
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
#import pandas.dtat

printEventByEvent = True #True
interpolate = False #True

plt.style.use(["science", "notebook", "grid"])
#Storage for the source position:
source_xyz =[]
source_Rtp = []
source_fraction_rawQ = []
source_fraction_Q = []
source_rawNHits = []
source_NHits = []
source_faction_evDigiHits = []
source_mean_evDigiQ = []
source_faction_evRawHits = []
source_mean_evRawQ = []
source_fraction_NHits = []
source_fraction_rawNHits = []


#the figure
ax = plt.axes(projection ="3d")

#All of the inputs
for name in sys.argv[1:]:
    name = str(name)

    #We can read a lot of information from the file name
    run_id = int(name.split("FileID")[1].split("_")[0])
    x = float(name.split("_x")[1].split("_")[0])
    y = float(name.split("_y")[1].split("_")[0])
    z = float(name.split("_z")[1].split("_")[0])
    theta = float(name.split("_t")[1].split("_")[0])
    phi = float(name.split("_p")[1].split("_")[0])
    R = 220#float(name.split("_R")[1].split("_")[0])

    print("source x, y, z: ", x, y, z, "R, theta, phi : ", theta, phi, R)

    #Open all the relevant branches of the files in the given config
    #tree_data = uproot.open("%s"%name)["CherenkovDigiHits"]
    #df_data = tree_data.arrays(library="pd")

    #To check the parentID of the particles I am simulating
    #tree_data = uproot.open("%s"%name)["Tracks"]
    #df_data = tree_data.arrays(library="pd")
    #df_data = df_data[abs(df_data["ParentID"])<=1000]

    #print(df_data["ParentID"].unique())
    #raise END

    tree_raw = uproot.open("%s"%name)["CherenkovHits"]
    df_raw = tree_raw.arrays(library="pd")
    #This can be usefull to check that the source is correctly positioned
    #tree_source = uproot.open("%s"%name)["EventInfo"]
    #df_source = tree_source.arrays(library="pd")
    tree_geometry = uproot.open("%s"%name)["Geometry"]
    df_geometry = tree_geometry.arrays(library="pd")
    table = [df_geometry["mPMT"],df_geometry["mPMT_pmt"], df_geometry["x"],df_geometry["y"],df_geometry["z"]]
    np.savetxt('/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/PMT_positions.txt',table, fmt = '%.2f')
    #now make a new dataset with only the information of the 58mPMT - the one of interest here
    #We have an issue that the 19th PMT is matched with the wrong mPMT - need to correct for that
    #pretty slow process -  could optimise with a separate function
    #df_data_58 = df_data[df_data["mPMT"]==58]
    #df_data_58 = df_data_58[df_data_58["mPMT_pmt"]!=19]
    #df_data_59 = df_data[df_data["mPMT"]==59]
    #df_data_58 = df_data_58.append(df_data_59[df_data_59["mPMT_pmt"]==19])

    df_raw_58 = df_raw[df_raw["mPMT"]==58]
    df_raw_58 = df_raw_58[df_raw_58["mPMT_pmt"]!=19]
    df_raw_59 = df_raw[df_raw["mPMT"]==59]
    df_raw_58 = df_raw_58.append(df_raw_59[df_raw_59["mPMT_pmt"]==19])

    df_geometry_58 = df_geometry[df_geometry["mPMT"]==58]
    df_geometry_58 = df_geometry_58[df_geometry_58["mPMT_pmt"]!=19]
    df_geometry_59 = df_geometry[df_geometry["mPMT"]==59]
    df_geometry_58 = df_geometry_58.append(df_geometry_59[df_geometry_59["mPMT_pmt"]==19])

    nEvents = df_raw["Event"].max()
    print(type(nEvents))
    if type(nEvents) != float:
        nEvents = int(nEvents)
    else:
        nEvents = 0

    nRecordedEvents = nEvents
    nEvents = 1000

    if (printEventByEvent):
        #here if we want an output of the charge and number of hits collected by our mPMT for each event to compare variations
        #for ev in range(1, 11):
            ##df_data_buf = df_data[df_data["Event"] == ev]

            ##print("EVENT NUMBER: %i Source R theta(deg) phi(deg) [%.2f, %.2f, %.2f]"%(ev, R, theta * 180/np.pi, phi * 180/np.pi))
            ##print("Theta = ", theta)
            ##df_data_buf_58 = df_data_58[df_data_58["Event"] == ev]
            ##print("NHits DIGI:", sum(df_data_buf["NDigiHits"]))
            ##print("Sum of NHits DIGI: ", sum(np.array(df_data_buf_58["NDigiHits"])))
            ##print("QTotDigi :", df_data_buf["QTotDigi"].unique())
            ##print("Sum of Q: ", sum(np.array(df_data_buf_58["Q"])))
            ##print("Percentage of the total charge in the 58th mPMT: ", sum(np.array(df_data_buf_58["Q"]))*100/df_data_buf["QTotDigi"].unique(), "\n")


            #df_raw_buf = df_raw[df_raw["Event"] == ev]
            #df_raw_buf_58 = df_raw_58[df_raw_58["Event"] == ev]
        print("nRecorded events", nRecordedEvents)
        df_raw_buf = df_raw
        df_raw_buf_58 = df_raw_58

        print("Total NHits RAW:", sum(df_raw_buf["NHits_noDN"]))

        print("Sum of NHits RAW in mPMT 58: ", sum(np.array(df_raw_buf_58["NHits_noDN"])))

        print("Total QTot RAW:", sum(df_raw_buf["PMT_QTot"]))

        print("Sum of Q RAW in mPMT 58: ", sum(np.array(df_raw_buf_58["PMT_QTot"])))
        print("Percentage of the RAW total charge in the 58th mPMT: ", sum(np.array(df_raw_buf_58["PMT_QTot"]))*100/max(sum(df_raw_buf["PMT_QTot"]), 0.001), "\n\n")


    #nEvent_DigiNHits_nonZero = 0
    #nEvent_DigiNHits58_nonZero = 0



    if nRecordedEvents != 0 :
        nEvent_RawNHits_nonZero = 0
        nEvent_RawNHits58_nonZero = 0

        for ev in range(1, nEvents+1):
            #loop over events - would be faster in cpp...
                #df_data_buf = df_data[df_data["Event"] == ev]
                #df_data_buf_58 = df_data_58[df_data_58["Event"] == ev]
                #if sum(df_data_buf["NDigiHits"]) != 0:
                    #nEvent_DigiNHits_nonZero += 1
                #if sum(df_data_buf_58["NDigiHits"]) != 0:
                    #nEvent_DigiNHits58_nonZero += 1

                df_raw_buf = df_raw[df_raw["Event"] == ev]
                df_raw_buf_58 = df_raw_58[df_raw_58["Event"] == ev]

                if sum(df_raw_buf["NHits"]) != 0:
                    nEvent_RawNHits_nonZero += 1
                if sum(df_raw_buf_58["NHits"]) != 0:
                    nEvent_RawNHits58_nonZero += 1

        source_xyz.append([x, y, z])
        source_Rtp.append([R,theta, phi])


        #The fraction of events where we see (digitised) hits
        #source_faction_evDigiHits.append(nEvent_DigiNHits58_nonZero/nEvent_DigiNHits_nonZero)
        source_faction_evRawHits.append(nEvent_RawNHits58_nonZero/nEvent_RawNHits_nonZero)

        #fraction of (digitised) charge that was detected by mPMT 58
        source_fraction_rawQ.append(sum(df_raw_58["PMT_QTot"])/sum(df_raw["PMT_QTot"]))
        #source_fraction_Q.append(sum(np.array(df_data_58["Q"]))/sum(df_data["QTotDigi"]))

        #fraction of (digitised) hits that was detected by mPMT 58
        #source_fraction_NHits.append(sum(df_data_58["NDigiHits"])/sum(df_data["NDigiHits"]))
        source_fraction_rawNHits.append(sum(df_raw_58["NHits"])/sum(df_raw["NHits"]))

        #number of (digitised) hits that was detected by mPMT 58
        #source_NHits.append(sum(df_data_58["NDigiHits"]))
        source_rawNHits.append(sum(df_raw_58["NHits"]))

        #mean number of (digitised) hits that was detected by mPMT 58
        if len(df_raw_58["PMT_QTot"]) != 0:
            source_mean_evRawQ.append(sum(df_raw_58["PMT_QTot"])/(nEvents))
        else:
            source_mean_evRawQ.append(0)

        #if len(df_data_58["Q"]) != 0:
            #source_mean_evDigiQ.append(sum(df_data_58["Q"])/(nEvents))
        #else:
            #source_mean_evDigiQ.append(0)


        ax.scatter3D(x, z, y, marker = "x", color = 'r')
        r = R
        #need to put the light rays - following the direction the way it is calculated for the WCSim .mac file:
        ax.plot([x, x - r *  np.sin(theta) * np.cos(phi) ], [z, z - r * np.sin(theta) * np.sin(phi)], [y, y - r * np.cos(theta)], 'r-')

    if nRecordedEvents == 0:
        source_xyz.append([x, y, z])
        source_Rtp.append([R,theta, phi])
        source_faction_evRawHits.append(0)
        source_fraction_rawQ.append(0)
        source_fraction_rawNHits.append(0)
        source_rawNHits.append(0)
        source_mean_evRawQ.append(0)
        ax.scatter3D(x, z, y, marker = "x", color = 'r')
        r = R
        ax.plot([x, x - r *  np.sin(theta) * np.cos(phi) ], [z, z - r * np.sin(theta) * np.sin(phi)], [y, y - r * np.cos(theta)], 'r-')


#For now, save the data in a txt file for ease of access
table = [np.array(source_xyz).T[:][0], np.array(source_xyz).T[:][1], np.array(source_xyz).T[:][2], np.array(source_Rtp).T[:][0], np.array(source_Rtp).T[:][1], np.array(source_Rtp).T[:][2], source_faction_evRawHits, source_fraction_rawQ, source_rawNHits,source_mean_evRawQ]

np.savetxt("../maps_txtFiles/map_raw_FileID%i.txt"%run_id, table, fmt = "%.2f")


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
#these are the mPMTs
ax.scatter3D(df_geometry_58["x"],df_geometry_58["z"],df_geometry_58["y"], marker = 'o')
plt.savefig("../maps_pictures/Source_positions_FileID%i.png"%run_id)

plt.show()


fig = plt.figure()
ax = fig.add_subplot(projection='polar')
ax.set_title("Q")
#plt.subplot(1,1,1)
x = np.array(source_Rtp).T[:][2]
y = np.array(source_Rtp).T[:][1]
z = source_mean_evRawQ
xi, yi = np.linspace(min(x), max(x), 100), np.linspace(min(y), max(y), 100)
xi, yi = np.meshgrid(xi, yi)
# Interpolate
rbf = scipy.interpolate.Rbf(x, y, z, function='cubic')
zi = rbf(xi, yi)
ax.contourf(xi, yi, zi, 500, cmap='nipy_spectral')
#ax.clim(0,1)

im = ax.scatter(np.array(source_Rtp).T[:][2],np.array(source_Rtp).T[:][1], c = source_mean_evRawQ, cmap = "nipy_spectral",s = 40)
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
plt.savefig("../maps_pictures/mPMT58_Efficiency_Map_FileID%i.png"%run_id)
plt.show()




fig = plt.figure()

ax = fig.add_subplot(projection='polar')
ax.set_title("NHits")
#plt.subplot(1,1,1)
x = np.array(source_Rtp).T[:][2]
y = np.array(source_Rtp).T[:][1]
z = source_faction_evRawHits
xi, yi = np.linspace(min(x), max(x), 100), np.linspace(min(y), max(y), 100)
xi, yi = np.meshgrid(xi, yi)
# Interpolate
rbf = scipy.interpolate.Rbf(x, y, z, function='cubic')
zi = rbf(xi, yi)
ax.contourf(xi, yi, zi, 500, cmap='nipy_spectral')
#ax.clim(0,1)
im = ax.scatter(np.array(source_Rtp).T[:][2],np.array(source_Rtp).T[:][1], c = source_faction_evRawHits, cmap = "nipy_spectral",s = 40)
col2 = plt.colorbar(im, label=f"Fraction of events for which there is at least 1 *raw* hit (i.e. charge) in mPMT58", orientation="horizontal", format= "%.1f")
#col2.set_ticks(levels)
#ax.clim(0,1)

ax.set_thetamin(0)
ax.set_thetamax(90)
##ax.set_phimin(0)
##ax.set_phimax(90)
#im.figure.axes[1].tick_params(axis="x", labelsize=12)
#ax.set_ylabel(f'$\phi()$')
ax.set_xlabel(f'$\Theta$(rad)')
plt.show()

#Plot the fraction of all the detector hits that fall onto the 58th mPMT

fig = plt.figure()
ax = fig.add_subplot(projection='polar')
ax.set_title("rawQ")
#plt.subplot(1,1,1)
x = np.array(source_Rtp).T[:][2]
y = np.array(source_Rtp).T[:][1]
z = source_fraction_rawQ
xi, yi = np.linspace(min(x), max(x), 100), np.linspace(min(y), max(y), 100)
xi, yi = np.meshgrid(xi, yi)
# Interpolate
rbf = scipy.interpolate.Rbf(x, y, z, function='cubic')
zi = rbf(xi, yi)
ax.contourf(xi, yi, zi, 500, cmap='nipy_spectral')
#ax.clim(0,1)

im = ax.scatter(np.array(source_Rtp).T[:][2],np.array(source_Rtp).T[:][1], c = source_fraction_rawQ, cmap = "nipy_spectral",s = 40)
col2 = plt.colorbar(im, label=f"Fraction of raw Q in mPMT 58", orientation="horizontal", format= "%.1f")
#col2.set_ticks(levels)
#ax.clim(0,1)

ax.set_thetamin(0)
ax.set_thetamax(90)
##ax.set_phimin(0)
##ax.set_phimax(90)
#im.figure.axes[1].tick_params(axis="x", labelsize=12)
#ax.set_ylabel(f'$\phi()$')
ax.set_xlabel(f'$\Theta$(rad)')
plt.show()


