make_onePositionMap: make_onePositionMap.cc
	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_onePositionMap make_onePositionMap.cc
#this is a file that distributes the read-out of the map in a source position (Rtp)-dependant file
#the file just appends the information and this is then used for checking how this specific position
#behaves as the scattering and absorption lengths are varied
#Please be careful not to add configurations that are not "reference" 
#(i.e. good stats, only scattering or absorption toggled at once, sensible R distance)
#later on we will move away from this to got to a binned approach
ID=$1
VAR="$(ls /home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/data/*FileID""$ID""_*_flat*)"
make
echo "Please be sure that "$ID" is a reference file as it might impact the quality of your fit later on to have non-reference file in your look up table"
sleep 2s
for name in $VAR
do
    echo $name
    ./bin/make_oneBin_test -f $name -o "Maps/maps_oneBin/mPMT_map_ID"$ID".txt"
    echo
done
