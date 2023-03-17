#This is a batch file that stores as a one liner the max value fo the charge collected if there is no
#absorption or scattering  we need to be very careful to use the correct reference otherwise it will
#mess up our absorption length estimate
ID=$1
VAR="$(ls ./data/*FileID""$ID""_*_flat*)"
#make
echo "Please be sure that "$ID" is a reference file as it might impact the quality of your fit later on to have non-reference file in your look up table - This file shouldn't have either scattering or attenuation"
sleep 2s
for name in $VAR
do
    echo $name
    ./bin/makeRef_maxNbHitAtPos -f $name
    echo
done
