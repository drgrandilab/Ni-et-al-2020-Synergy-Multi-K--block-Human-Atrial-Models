
for i in `seq 0 1 11`  # test AF model = 0 only
do
	matlab -nodisplay -nodesktop -nosplash -r "run($i);exit"
done 
