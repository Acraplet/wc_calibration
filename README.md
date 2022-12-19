# wc_calibration
_"Les mPMTs sont des fleurs dont les PMTs seraient les pétales"_

WCTE calibration code - PMT timing response, angular response, water attenuation length etc... 

## Make config files
First step to this analysis is to obtain the mPMT maps. This is done through a combinaison of WCSim and the python analyis code. So far we are using the PMT raw data as a reference, calculating for each source position the ratio of the total charge collected by the 58th mPMT (at the centre of the bottom end cap of WCTE) over the total number of photons sent (usually 1,000).

The config files for WCSim are stored using a FileID which indicates which group of runs it belongs to. 
In the file name is stored the absoption (abwff) and Rayleigh scattering coefficients that are input into WCSim.
Then the file name saves the source position x, y, z coordinates in space as well as its relative angles theta and phi with respect to the mPMT centre and distance R to the mPMT surface

The config files are produced in wc_calibration/WCSim_configFiles/ and the corresponding tuning files produced in wc_calibration/WCSim_tuningFiles/ by calling 
```
python /config_files_prod/writeMacFile.py
```
with the options :

- f the (integer) ID of the run
- t the number of points that we want to simualte equally spaced in sin(theta) between the limits (0 rad < theta < 1.1 rad) -> this can be changed in the file
- p the number of points that we want to simulate equally spaced in phi between the limits (0, pi/2) -> this can be changed in the file
- e the number of events that we want to simulate
- a the absorption coefficient (abwff in WCSim)
- r the Rayleigh scattering coefficient (rayff in WCSim)
- R the distance between the source position and the mPMT dome surface
- d the number of points to simulate if we want a random draw of points in sin(theta), phi instead of having them equally spaced
- u the number of points to simulate uniformly accross the sphere 

Careful, -d and -u shouldn't be used together! 

The photon number shot per event, its wavelength (401.9nm = 3.08945eV), the detector type, QE, triggering process etc..  are all pre-set in the /config_files_prod/WCSim_template.txt which is the base for the .mac files production - it can be modified.

Alternatively, if you want to make many similar config files, use

```
bash make_all_config.sh abwff rayff ID0
```
which will make 7 config files with the input absorption and scattering params and IDs going from ID0 to ID0 + 7 with 20 points in theta, and phi and with R = 5, 10, 20, 40, 80, 160, 250 (used to be 320)cm. 

## Run WCSim on the batch system (efficiently) 

The config files should then be copied to my /vols/t2k/user/ac4317/WCTE/mPMTmapping/config_files folder, and the tuning files to /vols/t2k/user/ac4317/WCTE/mPMTmapping/tuning_files. Then follow the steps presented in the image below. 

<img src="https://github.com/Acraplet/wc_calibration/blob/main/personnal-notes/Batch_submit_sketch.jpg" width="400" />

## Extract the maps as .txt file from WCSim 

*C++ (recommended)*
After the WCSim files have been produced, copy the \_flat files in your home /wc\_calibration/mPMTmapping/data folder. In the /wc_calibration/mPMTmapping folder you can then call 
```
bash make_mPMTmap.sh ID
```
replacing ID with the run ID you want. This saved a .txt file named /wc_calibration/mPMTmapping/maps\_txtFiles/map\_FileIDxx.txt with the following entries:

1. Source x position
2. Source y position
3. Source z position
4. Theta angle from the mPMT axis of the source
5. Phi angle from the mPMT axis of the source
6. Source R distance from the mPMT surface
7. Total charge collected in mPMT 58
8. Total number of events
9. Abwff - absorption coefficient
10. Rayff - Rayleigh scattering coeffcient

## Plot and compare the maps


Once these maps have been made, they can be plotted individually from the mPMTmapping folder with 
```
bash plotting_mPMTmap.sh ID
```
Or compared with
```
bash plotting_mPMTmap.sh ID1 ID2
```
It can be usefull to have a look at the profile of the map averaged over phi to compare more than two files. This is done with
```
bash compare_mPMTmap.sh ID1 ... IDn 
```

Note: the source coordinates need to be identical for the comparision to be made. TODO: Add a way to check the mPMT symmetry. 
Some examples of maps are already present in the maps_txtFiles folder for testing/understanding purposes.

Some left-over python codes are available but not maintained to do similar things as the c++ code. 


