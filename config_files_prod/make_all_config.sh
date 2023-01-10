abs= 100000000000
w=$1
#now we simulate 100 source positions (i.e. about 60 source positions in a given bin)
for ray in 0.00222 0.00555 0.01111 0.02222 0.03888 0.05555 0.1111 0.1666 0.2777 0.3333
do
	echo $ray
	for R in 10 20 40 80 120 140 160 180 210 250
	do
    		python writeMacFile.py -R $R -a 100000000000 -r $ray -e 1000 -f $w -u 100 
    		((w+=1))
	done
done
