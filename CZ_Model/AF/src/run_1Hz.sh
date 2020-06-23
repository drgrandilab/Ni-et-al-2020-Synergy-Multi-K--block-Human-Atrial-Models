
run() {
echo $1
python run_AF_1Hz.py $1
sleep 5
}


for (( i = 0; i < 12; i++ )); do
	# if  (($i % 6 == 0 )); then
	# 	wait
	# fi
	# statements
	mkdir AF_Model_$i
	# mkdir AF_Model_$i/Block_Model_$j
	cp NCZ_Model run_AF_1Hz.py AF_Model_$i/
	cd AF_Model_$i/
	run $i &
	cd -

done
wait


