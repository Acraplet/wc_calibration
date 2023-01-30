abs= 100000000000
w=$1
#now we simulate 50000 source positions (i.e. about 60 source positions in a given bin)
for ray in 0.1388 0.1666 0.19435 0.2222 0.24987 0.2777 0.30535
	#0.00555 0.01111 0.02222 0.03888 0.05555 0.1111 0.1666 0.2222 0.2777 0.3333 100000000000 
	#0.1388 0.1666 0.19435 0.2222 0.24987 0.2777 0.30535 #test
	#100000000000 0.00555 0.01111 0.02222 0.03888 0.05555 0.1111 0.1666 0.2222 0.2777 0.3333 #ref
	#0.08322 0.13888 0.19444 0.24999 #these are the tests (other tests)
#0.00222 0.00555 0.01111 0.02222 0.03888 0.05555 0.1111 0.1666 0.2777 0.3333 #these are the references
do
	echo $ray
	for R in 15 30 60 100 130 150 170 190 200 230 #- these are the test distances
		#10 20 40 80 120 140 160 180 210 250
		#15 30 60 100 130 150 170 190 200 230 - these are the test distances
		#10 20 40 80 120 140 160 180 210 250 - these are the reference R distances
	do
    		python writeMacFile.py -R $R -a 100000000000 -r $ray -e 1000 -f $w -u 100 
    		((w+=1))
	done
done
