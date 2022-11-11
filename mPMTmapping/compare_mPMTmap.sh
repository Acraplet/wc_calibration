#this is a called code to plot and save the plot of the effeciency map and comparision between different files
output_name='comparision_1d'

for var in "$@" #cycle through the imputs that we want
do
    name="$name"" maps_txtFiles/mPMT_map_ID""$var"".txt"
    echo $name
    output_name="$output_name""_""$var"
done

#output_name="$output_name"

python src/compare_txtMap.py $name $output_name
