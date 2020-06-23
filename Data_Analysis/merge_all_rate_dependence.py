
# merge existing .pkl files into a single pkl file.


import pandas as pd
import numpy as np


dat = []

for i in range(12):
	dat.append(pd.read_pickle('APD_Data/Rate_dependence_1_%d.pkl'%i))

res= pd.concat(dat, axis=0)



res.to_pickle('APD_Data/Merged_Rate_dependence_1.pkl')


res.to_csv('APD_Data/Merged_Rate_dependence_1.csv', sep=',',na_rep='NAN')


dat = []

for i in range(12):
	dat.append(pd.read_pickle('APD_Data/Rate_dependence_3_%d.pkl'%i))

res= pd.concat(dat, axis=0)



res.to_pickle('APD_Data/Merged_Rate_dependence_3.pkl')


res.to_csv('APD_Data/Merged_Rate_dependence_3.csv', sep=',',na_rep='NAN')


# for i in range(12):
# 	print 

# print np.mean(res.loc['Baseline']), np.mean(res.loc['Baseline'])



print res.describe()
