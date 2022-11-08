# wc_calibration
WCTE calibration code - PMT timing response, angular response, water attenuation length etc... 

First step to this analysis is to obtain the mPMT maps. This is done through a combinaison of WCSim and the python analyis code. So far we are using the PMT raw data as a reference, calculating for each source position the ratio of the total charge collected by the 58th mPMT (at the centre of the bottom end cap of WCTE) over the total number of photons sent (usually 1,000).

The config files for WCSim are stored using a FileID which indicates which group of runs it belongs to. 

It also has a label for the attenuation length we simulated in WCSim:
- noAlpha: Absff = 10e10, Rayff = 10e10, Mieff = 0.0
- 10cmScat: 
- 10cmAbs:
- 10cmScat10cmAbs: 

Then the file name saves the source position x, y, z coordinates in space as well as its relative angles theta and phi with respect to the mPMT centre and distance R to the mPMT surface

The config files are produced in wc\_calibration/WCSim\_configFiles/ by calling python /config\_files\_prod/writeMacFile.py with the options :

- f the (integer) ID of the run
- t the number of points that we want to simualte equally spaced in sin(theta) between the limits (0-1) -> this can be changed in the file
- p the number of points that we want to simulate equally spaced in phi between the limits (0, pi/2) -> this can be changed in the file
- e the number of events that we want to simulate

The photon number shot (1) per event, its wavelength (401.9nm = 3.08945eV), the detector type, QE, triggering process etc.. is all pre-set in the /config\_files\_prod/WCSim\_template.txt which is the base for the .mac files production - it can be modified.


The config files should then be copied to my /vols/t2k/user/ac4317/WCTE/mPMTmapping/config\_files folder. 
Then, from /vols/t2k/user/ac4317/WCTE/, modify the submit\_mPMTmapping.sh file to get the right mPMT ID. Be careful to re-compile the WCSim code by running source /vols/t2k/user/ac4317/WCTE/source\_at\_start\_test.sh and then make in the main WCTE folder. Be especially careful of this if you change things in the tunning\_paramters.mac file as these will only be applied then. Wait for the batch jobs to be in the reading mode r before running make WCSim again.

*PYTHON*
After the WCSim files have been produced, copy the \_flat files in your home /wc\_calibration/mPMTmapping/data folder where you can then run python /wc\_calibration/mPMTmapping/src/make\_raw\_mPMTmap.py $(ls \*IDxx\*) replacing xx by your desired file ID. This will plot and save in the /wc\_calibration/mPMTmapping/maps\_txtFiles folder the positions of the source in a 3D plot, and a map of the average recorded number of p.e. per event [be careful to modify the code if you do not run with the default 1000 events with a single photon sent each time] the mPMT then save a .txt file with the following entries:
1. Source x position
2. Source y position
3. Source z position
4. Source R distance from the mPMT surface
5. Theta angle from the mPMT axis of the source
6. Phi angle from the mPMT axis of the source
7. Fraction of events that recieved at least one raw hit
8. Fraction of the total detected charge that was detected in the 58th mPMT
9. Mean charge collected by the 58th mPMT per event

Once these maps have been made, they can be plotted at will with the following command python /wc\_calibration/mPMTmapping/src/read\_raw\_mPMTmap.py /wc\_calibration/mPMTmapping/maps\_txtFiles/map\_raw\_FileIDxx.txt where again the xx need to be swapped for the run ID you are interested in. Some examples of maps are already present in that folder for testing/understanding purposes

*C++ (recommended)*
After the WCSim files have been produced, copy the \_flat files in your home /wc\_calibration/mPMTmapping/data folder. In the /wc_calibration/mPMTmapping folder you can then call bash make_mPMTmap.sh _ID replacing _ID with the run ID you want. This saved a .txt file named /wc_calibration/mPMTmapping/maps\_txtFiles/map\_FileIDxx.txt with the following entries:

1. Source x position
2. Source y position
3. Source z position
4. Theta angle from the mPMT axis of the source
5. Phi angle from the mPMT axis of the source
6. Source R distance from the mPMT surface
7. Total charge collected in mPMT 58
8. Total number of events




Once these maps have been made, they can be plotted at will with the following command python /wc\_calibration/mPMTmapping/src/read\_txtMap.py /wc\_calibration/mPMTmapping/maps\_txtFiles/map\_FileIDxx.txt where again the xx need to be swapped for the run ID you are interested in. Some examples of maps are already present in that folder for testing/understanding purposes

