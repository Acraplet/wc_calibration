#This is a usefull caller code for making a lot of mPMT maps at once - from i to f

export WCCALIB=$PWD # the wc_calibration directory
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WCCALIB/mPMTmapping/lib
source ../source_at_start_test.sh
cd mPMTmapping/
unset WCSIMDIR

make

for a in 1210 1211 1212 1213
#1201 1202 1203 1204 1205 1206 1207 1208 1209
do
    echo $a
    bash make_mPMTmap.sh $a
done
