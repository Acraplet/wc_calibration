#this is a called code to plot and save the plot of the effeciency map and comparision between different files

ID_1=$1
mode=cubic

if [[ $# -eq 2 ]]
then
    ID_2=$2
    echo "Plotting comparision of the two files "$ID_1" and "$ID_2""
    python src/read_txtMap.py -f maps_txtFiles/mPMT_map_ID$ID_1.txt -c maps_txtFiles/mPMT_map_ID$ID_2.txt -o "comparision_FileID"$ID_1"_FileID"$ID_2""$mode""
fi

if [[ $# -eq 1 ]]
then
    echo "Plotting the map of the file "$ID_1""
    python src/read_txtMap.py -f maps_txtFiles/mPMT_map_ID$ID_1.txt -o "plotting_FileID"$ID_1""$mode""
fi
