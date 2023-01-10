#This is a usefull caller code for making a lot of mPMT maps at once
for a in {300..399}
do
    echo $a
    bash make_mPMTmap.sh $a
done
