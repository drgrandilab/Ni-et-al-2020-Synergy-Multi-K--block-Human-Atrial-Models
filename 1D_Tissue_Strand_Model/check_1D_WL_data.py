

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
data_cv = []

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
	data_cv.append(data_pd['cv2'])


res = pd.concat (data_apd,axis=1, keys=['Baseline', 'IK2p', 'ISK', 'IK2p_ISK', 'IKur', 'IK2p_IKur', 'ISK_IKur','IK2P_ISK_IKur'])
res_cv = pd.concat (data_cv,axis=1, keys=['Baseline', 'IK2p', 'ISK', 'IK2p_ISK', 'IKur', 'IK2p_IKur', 'ISK_IKur','IK2P_ISK_IKur'])
res.to_pickle('Summarise_APD_Model_1.pkl')
res.to_csv('Summarise_APD_Model_1.csv')
res_cv.to_pickle('Summarise_CV_Model_1.pkl')
res_cv.to_csv('Summarise_CV_Model_1.csv')

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

import matplotlib as mpl
rc('mathtext',default='regular')
fonts = 20-0;
fonts_title = 20-0
leng_fontsize = 17-0
mpl.rcParams['font.size'] = fonts
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['axes.labelsize'] = fonts
mpl.rcParams['xtick.labelsize'] = fonts
mpl.rcParams['ytick.labelsize'] = fonts
mpl.rcParams['axes.linewidth'] = 0
mpl.rcParams['ytick.major.pad']='8'
mpl.rcParams['xtick.major.pad']='10'
lnwidth = 0.6*np.array( [3.5, 		3.5 ,     3.5,       3.5,      3.5,      3.5,          3.5]);
colors = ['k',     'r',    'b',     'b',   'm', 'cyan',  'magenta'];
markers = 7*'v'

fig = plt.figure(figsize=(1.1*0.5*19*1.5/2.54, 19*1.2/2/2.54))
num_row = 1;
num_col = 1;
gs1 = gridspec.GridSpec(num_row,num_col
	# width_ratios=[1.15,1],
 #    height_ratios=[1,1]
	);
panel = {};

for i in np.arange(num_row):
	for j in np.arange(num_col):
		if j>0:
			panel[i,j] = plt.subplot(gs1[i,j], sharey = panel[i,0]);
		else :
			panel[i,j] = plt.subplot(gs1[i,j]);

		ax = panel[i,j];

# label = ['I$_{K2P}$','I$_{K,Ca}$','I$_{Kur}$', 'I$_{K2P}$ + I$_{K,Ca}$', 'I$_{K2P}$ + I$_{Kur}$', 'I$_{K,Ca}$ + I$_{Kur}$', 'I$_{K2P}$ + I$_{K,Ca}$ + I$_{Kur}$']
label = ['Ctl', 'I$_{K2P}$', 'I$_{KCa}$', 'I$_{Kur}$', 'I$_{K2P}$+I$_{KCa}$', 'I$_{K2P}$+I$_{Kur}$', 'I$_{KCa}$+I$_{Kur}$', 'I$_{K2P}$+I$_{KCa}$+I$_{Kur}$']

# label = ['I$_{K2P}$ + I$_{SK}$', 'I$_{K2P}$ + I$_{Kur}$', 'I$_{SK}$ + I$_{Kur}$', 'I$_{K2P}$ + I$_{SK}$ + I$_{Kur}$']
col = 0;
row = 0;

# for i in range(3):
	# col=i
ax = panel[row,col];

keys=['Baseline', 'IK2p', 'ISK',  'IKur', 'IK2p_ISK', 'IK2p_IKur', 'ISK_IKur','IK2P_ISK_IKur']


bp = ax.boxplot([res[col] for col in keys],showfliers=False)
numBoxes=len(keys)
for i in range(numBoxes):

    y = res[keys[i]]
    x = np.random.normal(1+i, 0.04, size=len(y))
    ax.plot(x, y, 'b.', alpha=0.3)
for box in bp['boxes']:
	 box.set(  linewidth=2)
for whisker in bp['whiskers']:
	whisker.set(  linewidth=2)


for cap in bp['caps']:
    cap.set( linewidth=2)
for median in bp['medians']:
    median.set(linewidth=2)





for i in np.arange(num_col*num_row):
	ax = fig.axes[i]
	# ax.autoscale(enable=True, axis='x')  
	# ax.autoscale(enable=True, axis='r')  
	ax.spines['right'].set_visible(False)
	#ax.spines['left'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(True)
	ax.xaxis.set_ticks_position('bottom')
	ax.spines['bottom'].set_linewidth(1.5)
	ax.spines['left'].set_linewidth(1.5)
	ax.spines['right'].set_linewidth(1.5)
	ax.yaxis.tick_left()
	# ax.set_xlim(0.5+0.2,7.8-0.1)
	# # plt.xticks(np.linspace(0,600, 5))
	# ax.set_yticks(np.linspace(-5, 35, 9))
	ax.yaxis.set_tick_params(width=1.5)
	ax.xaxis.set_tick_params(width=1.5)
	plt.subplots_adjust(left=0.13, bottom=0.15, right=0.95, top=0.95,
		wspace=0.32, hspace=0.16)
	# ax.yaxis.set_label_coords(-0.1, 0.5)

ax.set_xticks(np.linspace(1,8,8))
ax.get_xaxis().set_ticklabels(label, rotation = 15)



plt.show()
