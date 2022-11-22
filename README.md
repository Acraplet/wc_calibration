# wc_calibration
_"Les mPMTs sont des fleurs dont les PMTs seraient les pétales"_

WCTE calibration code - PMT timing response, angular response, water attenuation length etc... 

## Make config files
First step to this analysis is to obtain the mPMT maps. This is done through a combinaison of WCSim and the python analyis code. So far we are using the PMT raw data as a reference, calculating for each source position the ratio of the total charge collected by the 58th mPMT (at the centre of the bottom end cap of WCTE) over the total number of photons sent (usually 1,000).

The config files for WCSim are stored using a FileID which indicates which group of runs it belongs to. 

It also has a label for the attenuation length we simulated in WCSim:
- noAlpha: Absff = 10e10, Rayff = 10e10, Mieff = 0.0
- 10cmScat: Absff = 10e10, Rayff = 0.000555, Mieff = 0.0
- 10cmAbs: Absff = 0.000243, Rayff = 10e10, Mieff = 0.0
- 10cmScat10cmAbs: Absff = 0.000243, Rayff = 0.000555, Mieff = 0.0
- 20cmScat20cmAbs: Absff = 0.000486, Rayff = 0.00111, Mieff = 0.0

Then the file name saves the source position x, y, z coordinates in space as well as its relative angles theta and phi with respect to the mPMT centre and distance R to the mPMT surface

The config files are produced in wc_calibration/WCSim_configFiles/ by calling 
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

The photon number shot (1) per event, its wavelength (401.9nm = 3.08945eV), the detector type, QE, triggering process etc.. is all pre-set in the /config_files_prod/WCSim_template.txt which is the base for the .mac files production - it can be modified.

Alternatively, if you want to make many similar config files, use

```
bash make_all_config.sh abwff rayff ID0
```
which will make 7 config files with the input absorption and scattering params and IDs going from ID0 to ID0 + 7 with 20 points in theta, and phi and with R = 5, 10, 20, 40, 80, 160, 320cm. 

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


## Python alternative

*PYTHON (not recommended)*
After the WCSim files have been produced, copy the \_flat files in your home /wc_calibration/mPMTmapping/data folder where you can then run 
```
python /wc_calibration/mPMTmapping/src/make_raw_mPMTmap.py $(ls \*IDxx\*)
```
replacing xx by your desired file ID. This will plot and save in the /wc_calibration/mPMTmapping/maps_txtFiles folder the positions of the source in a 3D plot, and a map of the average recorded number of p.e. per event [be careful to modify the code if you do not run with the default 1000 events with a single photon sent each time] the mPMT then save a .txt file with the following entries:
1. Source x position
2. Source y position
3. Source z position
4. Source R distance from the mPMT surface
5. Theta angle from the mPMT axis of the source
6. Phi angle from the mPMT axis of the source
7. Fraction of events that recieved at least one raw hit
8. Fraction of the total detected charge that was detected in the 58th mPMT
9. Mean charge collected by the 58th mPMT per event

Once these maps have been made, they can be plotted at will with the following command 
```
python /wc_calibration/mPMTmapping/src/read_raw_mPMTmap.py /wc_calibration/mPMTmapping/maps_txtFiles/map_raw_FileIDxx.txt
```
where again the xx need to be swapped for the run ID you are interested in. Some examples of maps are already present in that folder for testing/understanding purposes


