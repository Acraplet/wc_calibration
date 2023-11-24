w=$1
abs=10e10
#these are the rayleigh coefficients that I use for reference 
for R in 10 40 60 80 100 120 140 160 180 210    
#0.00055 0.00111 0.002222 0.003333 0.00555 0.008325 0.01221 0.00445 0.0067 0.00945 0.0106 
do
    python writeMacFile.py -R $R -a 10e10 -r 10e10 -e 1000 -f $w -u 4000
    ((w+=1))
done
