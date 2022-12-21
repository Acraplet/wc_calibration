#This is a usefull caller code for making a lot of mPMT maps at once
for a in {175..202}
do
    echo $a
    bash make_mPMTmap.sh $a
done
