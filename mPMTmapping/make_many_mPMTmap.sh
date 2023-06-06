#This is a usefull caller code for making a lot of mPMT maps at once - from i to f
for a in {1176..1199}
do
    echo $a
    bash make_mPMTmap.sh $a
done
