#this is a called code to plot and save the plot of the effeciency map and comparision between different files

ID_1=$1
mode=binned

if [[ $# -eq 2 ]]
then
    ID_2=$2
    echo "Plotting comparision of the two files "$ID_1" and "$ID_2""
    python src/read_txtMap_binned.py -f Maps/maps_txtFiles/mPMT_map_ID$ID_1.txt -f Maps/maps_txtFiles/mPMT_map_ID$ID_2.txt -o "comparision_FileID"$ID_1"_FileID"$ID_2""$mode""
fi

if [[ $# -eq 3 ]]
then
    ID_2=$2
    ID_3=$3
    echo "Plotting comparision of the two files "$ID_1" and "$ID_2""
    python src/read_txtMap_binned.py -f Maps/maps_txtFiles/mPMT_map_ID$ID_1.txt -f Maps/maps_txtFiles/mPMT_map_ID$ID_2.txt -f Maps/maps_txtFiles/mPMT_map_ID$ID_3.txt -o "comparision_FileID"$ID_1"_FileID"$ID_2"_FileID"$ID_3""$mode""
fi

if [[ $# -eq 4 ]]
then
    ID_2=$2
    ID_3=$3
    ID_4=$4
    echo "Plotting comparision of the two files "$ID_1" and "$ID_2""
    python src/read_txtMap_binned.py -f Maps/maps_txtFiles/mPMT_map_ID$ID_1.txt -f Maps/maps_txtFiles/mPMT_map_ID$ID_2.txt -f Maps/maps_txtFiles/mPMT_map_ID$ID_3.txt -f Maps/maps_txtFiles/mPMT_map_ID$ID_4.txt -o "comparision_FileID"$ID_1"_FileID"$ID_2"_FileID"$ID_3"_FileID"$ID_4""$mode""
fi

if [[ $# -eq 1 ]]
then
    echo "Plotting the map of the file "$ID_1""
    python src/read_txtMap_binned.py -f Maps/maps_txtFiles/mPMT_map_ID$ID_1.txt -o "plotting_FileID"$ID_1""$mode""
fi
