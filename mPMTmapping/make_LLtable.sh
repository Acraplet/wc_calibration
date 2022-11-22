#This is a caller file to store all the files comparisions
#it is based on the src/calc_LL.py code
rm table.txt
for ID in {14..37}
do
    for ID2 in {14..37}
    do
        if [ $ID -ne 20 ]
        then
            if [ $ID2 -ne 20 ]
            then
		echo "making: ref: "$ID"" ""$ID2""
		echo " "
                echo "$(python src/calc_LL.py -f $ID -c $ID2)" >> table.txt
            fi
        fi
    done
done

