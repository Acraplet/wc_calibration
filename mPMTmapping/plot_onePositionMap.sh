#this is a called code to plot and save the plot of the effeciency map and comparision between different files

for name in $(ls /home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_onePosition/*R40.00*)
do
    echo "Plotting the map of the file "$name""
    python src/read_txtOnePosition.py -f $name 
    #/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_onePosition/$name
done
