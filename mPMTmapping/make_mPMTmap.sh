#This is a batch file to run a complete map for a given run ID
#called with bash make_mPMTmap.sh ID where ID is the FileID that we want and the flat files are stored in the
#/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/data/ folder
ID=$1
VAR="$(ls ./data/*FileID""$ID""_*_flat*)"
#VAR="$(/vols/t2k/users/ac4317/WCTE/WCSim/mPMTmapping/data/*FileID""$ID""_*_flat*)"
make

#rm Maps/maps_txtFiles/mPMT_map_ID"$ID".txt
#to make sure we don't add to an existing file
for name in $VAR
do
    echo $name
    ./bin/make_mPMTmap -f $name -o "Maps/maps_txtFiles/mPMT_map_ID"$ID".txt"
    echo
done
