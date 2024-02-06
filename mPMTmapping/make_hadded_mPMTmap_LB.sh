#This is a batch file to run a complete map for a given run ID
#called with bash make_mPMTmap.sh ID where ID is the FileID that we want and the flat files are stored in the
#/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/data/ folder
#now we are adding to the same map all of the sub files that compose a bigger one (stats!) 
ID=$1
#VAR="$(ls ./data/data/*FileID""$ID""_*_flat*)"
#make

rm Maps/maps_txtFiles/mPMT_map_ID"$ID".txt
for subID in "${@:2}";
do
	#to make sure we don't add to an existing file
	VAR="$(ls /vols/t2k/users/ac4317/WCTE/wc_calibration/mPMTmapping/data/*FileID""$subID""_*_flat*)"
	for name in $VAR
	do
    		echo $name
    		./bin/make_mPMTmap_binned_LB -f $name -o "Maps/maps_txtFiles/mPMT_map_ID"$ID".txt"
    		echo
	done
done
