#This is a batch file to run a complete map for a given run ID
#called with bash make_mPMTmap.sh ID where ID is the FileID that we want and the flat files are stored in the
#/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/data/ folder
ID=$1
cone=$2
radius=$3
multiplicatorSurface=$4
effectiveRadius=$5
multiplicatorReflector=$6
echo "Please use this code by calling bash make_solidAngleMap_LB.sh runID coneSize PMTtrueRadius multiplicatorSurface PMTeffectiveRadius multiplicatorReflector"
#VAR="$(ls ./data/data/*FileID""$ID""_*_flat*)"
VAR="$(ls /vols/t2k/users/ac4317/WCTE/wc_calibration/mPMTmapping/data/*FileID""$ID""_*_flat*)"
#make

rm Maps/Expected_Number_Photons_FileID"$ID"_"$radius".txt
rm Maps/Expected_Number_Photons_FileID"$ID"_r"$radius"_e"$effectiveRadius"_m"$multiplicatorSurface"_n"$multiplicatorReflector".txt 
#to make sure we don't add to an existing file
for name in ${VAR[0]}
do
    echo $name
    ./bin/make_predNbPhotons_LB_withAR -f $name -c $cone -e $effectiveRadius -m $multiplicatorSurface  -n $multiplicatorReflector -r $radius -o "Maps/Expected_Number_Photons_FileID"$ID"_r"$radius"_e"$effectiveRadius"_m"$multiplicatorSurface"_n"$multiplicatorReflector".txt"
    echo
done
