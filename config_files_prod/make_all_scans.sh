R=$1
w=$2
abs=10e10
#these are the rayleigh coefficients that I use for reference 
for ray in 0.00055 0.00111 0.002222 0.003333 0.00555 0.008325 0.01221 0.00445 0.0067 0.00945 0.0106 
do
    python writeScanMacFile.py -R $R -a $abs -r $ray -e 10000 -f $w
    ((w+=1))
done
