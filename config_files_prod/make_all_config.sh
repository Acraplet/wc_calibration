abs=$1
ray=$2
w=$3
for R in 5 10 20 40 80 160 320
do
    python writeMacFile.py -R $R -a $abs -r $ray -e 1000 -p 20 -t 20 -f $w
    ((w+=1))
done
