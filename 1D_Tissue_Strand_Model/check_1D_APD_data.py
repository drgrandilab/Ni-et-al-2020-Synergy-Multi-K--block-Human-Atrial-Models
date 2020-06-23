

from subprocess import call
import math
# S1 = 500
# S2 = 250
import sys


import numpy as np

import pandas as pd


def get_data(model):
	data = np.loadtxt('OneD.log.dat.AF_Model.1.block.%d'%model)
	print data[1,1], len(data[:,1])
	index = []
	data_out = []
	for i in range(len(data[:,1])):
		tmp = data[i, 0];
		if tmp not in index:
			index.append(tmp)
			data_out.append(data[i, :])
	print index 
	print len(index)
	print data_out[2:5][:]

	return data_out

data_apd = []

for i in range(8):
	data_pd = pd.DataFrame(get_data(i))
	data_pd.columns = [ 'SIM_ID',
	'cv1',
	'cv2',
	'AP_Number1',
	'AP_Number2',
	'BCL_1',
	'APD_out_end',
	'APD75_out',
	'APD50_out',
	'APD30_out',
	'APD20_out',
	'Vmax_out',
	'Vmin_out',
	'dVdtmax_out',
	'plateau_potential_out',
	'INa_max_out',
	'dVdt_Time_RP',
	'V_m70_Time_RP',
	'Cai_max_out',
	'Cai_min_out',
	'BCL_2',
	'APD_out_end_2',
	'APD75_out_2',
	'APD50_out_2',
	'APD30_out_2',
	'APD20_out_2',
	'Vmax_out_2',
	'Vmin_out_2',
	'dVdtmax_out_2',
	'plateau_potential_out_2',
	'INa_max_out_2',
	'dVdt_Time_RP_2',
	'V_m70_Time_RP_2',
	'Cai_max_out_2',
	'Cai_min_out_2',
	'BCL_3',
	'APD_out_end_3',
	'APD75_out_3',
	'APD50_out_3',
	'APD30_out_3',
	'APD20_out_3',
	'Vmax_out_3',
	'Vmin_out_3',
	'dVdtmax_out_3',
	'plateau_potential_out_3',
	'INa_max_out_3',
	'dVdt_Time_RP_3',
	'V_m70_Time_RP_3',
	'Cai_max_out_3',
	'Cai_min_out_3',
	'BCL_4',
	'APD_out_end_4',
	'APD75_out_4',
	'APD50_out_4',
	'APD30_out_4',
	'APD20_out_4',
	'Vmax_out_4',
	'Vmin_out_4',
	'dVdtmax_out_4',
	'plateau_potential_out_4',
	'INa_max_out_4',
	'dVdt_Time_RP_4',
	'V_m70_Time_RP_4',
	'Cai_max_out_4',
	'Cai_min_out_4']
	data_pd.index = (data_pd['SIM_ID'])
	print data_pd.index.values
	print data_pd.columns
	data_pd.to_pickle('OneD.log.dat.AF_Model.1.block.%d.pkl'%i)
	data_pd.to_csv('OneD.log.dat.AF_Model.1.block.%d.csv'%i)

	data_apd.append(data_pd['APD_out_end_3'])


res = pd.concat (data_apd,axis=1, keys=['Baseline', 'IK2p', 'ISK', 'IK2p_ISK', 'IKur', 'IK2p_IKur', 'ISK_IKur','IK2P_ISK_IKur'])
res.to_pickle('Summarise_data_Model_1.pkl')
res.to_csv('Summarise_data_Model_1.csv')