
# merge existing .pkl files into a single pkl file.


import pandas as pd
import numpy as np
import numpy as np
import collections


import matplotlib

# matplotlib.use('Qt4Agg') 

import matplotlib.gridspec as gridspec

import matplotlib.pyplot as plt
from matplotlib import rc
# import pandas as pd

from random import *

# import csv
import sys



def analyse_data(Model):
	data_1Hz = pd.read_pickle('../1Hz/Summarise_data_Model_%d.pkl'%Model)
	data_3Hz = pd.read_pickle('../3Hz/Summarise_data_Model_%d.pkl'%Model)
	data_1Hz = data_1Hz.transpose()
	data_3Hz = data_3Hz.transpose()


	# print data_1Hz.describe()
	# print data_1Hz['Baseline']/data_1Hz['Baseline']
	# print data_3Hz.describe()

	APD_R1 = data_1Hz['Baseline']

	APD_R1_Baseline      = data_1Hz['Baseline']/data_1Hz['Baseline']
	APD_R1_IK2p          = data_1Hz['IK2p']/data_1Hz['Baseline']
	APD_R1_ISK           = data_1Hz['ISK']/data_1Hz['Baseline']
	APD_R1_IKur          = data_1Hz['IKur']/data_1Hz['Baseline']
	APD_R1_IK2p_ISK      = data_1Hz['IK2p_ISK']/data_1Hz['Baseline']
	APD_R1_IK2p_IKur     = data_1Hz['IK2p_IKur']/data_1Hz['Baseline']
	APD_R1_ISK_IKur      = data_1Hz['ISK_IKur']/data_1Hz['Baseline']
	APD_R1_ISK_IKur_Over_IKur      = data_1Hz['ISK_IKur']/data_1Hz['IKur']
	APD_R1_ISK_IKur_Over_ISK      = data_1Hz['ISK_IKur']/data_1Hz['ISK']
	APD_R1_IK2P_ISK_IKur = data_1Hz['IK2P_ISK_IKur']/data_1Hz['Baseline']

	# print APD_R1_Baseline.describe()


	f1 = [APD_R1_Baseline ,APD_R1_IK2p, APD_R1_ISK, APD_R1_IK2p_ISK, APD_R1_IKur,APD_R1_IK2p_IKur,APD_R1_ISK_IKur ,APD_R1_IK2P_ISK_IKur,APD_R1,APD_R1_ISK_IKur_Over_IKur,APD_R1_ISK_IKur_Over_ISK]

	for i in range(len(f1)):
		f1[i] [f1[i]<0] = np.nan
	res1= pd.concat( f1, axis=1,
	keys=['Baseline', 'IK2p', 'ISK', 'IK2p_ISK', 'IKur', 'IK2p_IKur', 'ISK_IKur','IK2P_ISK_IKur', 'APD_baseline', 'ISK_IKur_Over_IKur','ISK_IKur_Over_ISK']
	)

	# print res1.describe()

	APD_R3 = data_3Hz['Baseline']


	APD_R3_Baseline      = data_3Hz['Baseline']/data_3Hz['Baseline']
	APD_R3_IK2p          = data_3Hz['IK2p']/data_3Hz['Baseline']
	APD_R3_ISK           = data_3Hz['ISK']/data_3Hz['Baseline']
	APD_R3_IKur          = data_3Hz['IKur']/data_3Hz['Baseline']
	APD_R3_IK2p_ISK      = data_3Hz['IK2p_ISK']/data_3Hz['Baseline']
	APD_R3_IK2p_IKur     = data_3Hz['IK2p_IKur']/data_3Hz['Baseline']
	APD_R3_ISK_IKur      = data_3Hz['ISK_IKur']/data_3Hz['Baseline']
	APD_R3_IK2P_ISK_IKur = data_3Hz['IK2P_ISK_IKur']/data_3Hz['Baseline']
	APD_R3_ISK_IKur_Over_IKur      = data_3Hz['ISK_IKur']/data_3Hz['IKur']
	APD_R3_ISK_IKur_Over_ISK      = data_3Hz['ISK_IKur']/data_3Hz['ISK']
	# print APD_R1_Baseline.describe()
	f2 = [APD_R3_Baseline ,APD_R3_IK2p, APD_R3_ISK, APD_R3_IK2p_ISK, APD_R3_IKur,APD_R3_IK2p_IKur,APD_R3_ISK_IKur ,APD_R3_IK2P_ISK_IKur,APD_R3, APD_R3_ISK_IKur_Over_IKur,APD_R3_ISK_IKur_Over_ISK]
	for i in range(len(f1)):
		f2[i] [f2[i]<0] = np.nan
	res3= pd.concat(f2 , axis=1,
	keys=['Baseline', 'IK2p', 'ISK', 'IK2p_ISK', 'IKur', 'IK2p_IKur', 'ISK_IKur','IK2P_ISK_IKur', 'APD_baseline', 'ISK_IKur_Over_IKur','ISK_IKur_Over_ISK']
	)


	res = res3/res1;

	print res.describe()
	res.to_csv('APD_Data/Rate_dependence_%d.csv'%Model, sep=',',na_rep='NAN')
	res1.to_csv('APD_Data/Rate_dependence_1_%d.csv'%Model, sep=',',na_rep='NAN')
	res3.to_csv('APD_Data/Rate_dependence_3_%d.csv'%Model, sep=',',na_rep='NAN')
	res.to_pickle('APD_Data/Rate_dependence_%d.pkl'%Model)
	res1.to_pickle('APD_Data/Rate_dependence_1_%d.pkl'%Model)
	res3.to_pickle('APD_Data/Rate_dependence_3_%d.pkl'%Model)


for i in range(12):
	analyse_data(i)
