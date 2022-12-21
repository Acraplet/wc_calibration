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

## Extract and look at the maps as .txt file from WCSim 

*C++ (recommended)*
After the WCSim files have been produced, copy the \_flat files in your home /wc\_calibration/mPMTmapping/data folder. In the /wc_calibration/mPMTmapping folder you can then call 
```
bash make_mPMTmap.sh ID
```
replacing ID with the run ID you want. This saved a .txt file named /wc_calibration/mPMTmapping/Maps/maps\_txtFiles/map\_FileIDxx.txt with the following entries:

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
Note: the source coordinates need to be identical for the comparision to be made.

In later version of the code we will be looking at binned response to get away from edge effects. To have a look at the binned charge distribution for a given configuration do
```
bash plot_binned.sh ID
``` 
where the bins are stored in uniform_304_bins.txt and made via 
```
python src/get_uniform_bins.py
```
So far we are only looking at a quarter of the sphere but is is easy to modify this code to extrand the range of uniform bins.
Some examples of maps are already present in the maps_txtFiles folder for testing/understanding purposes.

Some left-over python codes are available but not maintained to do similar things as the c++ code. 

## Extract the absorption length

### Get the reference max amplitude
For absorption, the max charge that gets collected at a given position when there is no scattering and no absorption is saved in Maps/maps_Reference by the code
```
bash make_reference ID_ref
```
### Fit the response
The typical response function for absorption is an exponential function Q_pred = A_max * exp(-R/abs_length). A_max is the reference value that we have saved just before and the R value (dome-source distance) is known and lifted from the dataset. abs_length is obtained from minimising the chi2 between Q_true and Q_pred. The code that does this is
```
./bin/AbsorptionFitter ID1 ID2 ... IDn
```
The accuracy and precision of the fit is greatly improved if we use multiple maps (at different Rs but same abwff) together. Intuitively, the no-absorption, zero-scattering maps should be independant of R however it is not quite the case (mainly because of rounding errors). This issue should be resolved by using a binned approach but until then we need to have a reference max charge for each source position at each R that we want to use. 


## Extract the scattering length

### Obtain the data of the reference scattering behaviour

To understand how the charge collected by the mPMT is dependant on the scattering we need to collect reference data, at each (R, theta, phi) positions we have simulated data with a bunch of scattering lengths. These serve as reference points for the behaviour and are stacked together in Maps/maps_onePosition with 

```
bash make_onePositionMap.sh ID_ref
```
as well as files with a give rayff and an infinite abwff, we can also use this code on files with an infinite rayff and a given abwff this is useful later on to verify the hypothesis we made of there being an exponential response for the absorption coefficient 

### Obtain the response function of the mPMT to the scattering length

The code
```
./bin/Fitter ID
```
looks up the reference values in the OnePosition files and then fits separately the attenuation and scattering response to this data. For the absorption fit, an exponential with two free parameters is fitted to these reference points with Q_pred = A * exp(-R/abs_length_true). the output of this fit is saved in the mPMTmapping/reference_root file corresponding to this (R, theta, phi) position. A quick comparision between the mean estimate of R and A over all of the source positions and the corresponding true values, respectively input to WCSim and the A_max obtained from datasets without any attenuation are a check of the hypothesis we made previously. 

This code looks in a second time at the scattering length and fits to the reference datapoints a piece-wise continuous polynomial with fixed boundaries. To do this the code fits a fake data TSpline3 of which the nodes are fixed at a given x and the nodes' y are the fitted parameters for this fit. The density of nodes follows the magnitude of the gradient of the response which has a steeper evolution at low scattering lengths. the overall number of nodes is always smaller than the number of datapoints. The start- and end-points gradient to the spline are fixed respectively as the gradient between the two first data points and 0 as can be easily justified physically. 

### Use the response function to extract the scattering length
After the source position-dependant mPMT response to scattering has been fitted, we use it to extract the scattering length of any test mPMT map with

```
./bin/ScatteringFitter ID1 ID2 ... IDn
```

just like with absorption we try to minimise the difference between Q_pred and Q_true as a function of the scattering length but this time Q_pred = TSpline3_fakedata(*y_of_nodes_from above) at scattering_length_pred. Again we can use multiple files to refine the approach. unlike for absorption the results are not very convincing so far...


