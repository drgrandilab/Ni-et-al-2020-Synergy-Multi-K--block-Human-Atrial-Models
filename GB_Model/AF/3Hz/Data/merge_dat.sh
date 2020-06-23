

for M in `seq 0 1 11`; do 
for j in `seq 0 1 7`; do

for i in `seq 1 1 28`; do
	# echo $i

	cat ../AF_model_ID_$M/CPU_ID_$i\_test_drug_$j.dat >> merged_data.dat.model.$M.drug_$j

	# echo $i\_test.dat 
done

done


done 
