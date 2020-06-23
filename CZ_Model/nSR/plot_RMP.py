
import numpy as np
import collections


import matplotlib

# matplotlib.use('Qt4Agg') 

import matplotlib.gridspec as gridspec

import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd

from random import *

# import csv



# para_set = np.loadtxt('para_log.dat')





Data = np.loadtxt('1Hz/AP.log.dat.1Hz', unpack=True)
# Data2 = np.loadtxt('AF_Model_0/AP.log.dat.AF_Model.0.block.3', unpack=True)


N = len(Data)
print N




def get_difference(ori, changed):
	avg_by_cell = 0.5*(changed[16]+ changed[1]) - 0.5*(ori[16]+ ori[1])
	avg_by_cell_percentage = avg_by_cell / (0.5*(ori[16]+ ori[1]))

	return avg_by_cell, avg_by_cell_percentage


def check_alternans(Data):
	alternans= alternans = np.abs(Data[1] - Data[16])
	j=0
	for i in alternans:
		if i > 5:
			j=j+1
	print j
	return alternans




# plt.hist(Data[8], 25)
# plt.show()

sel =[]
Data2 = np.loadtxt('3Hz/AP.log.dat.3Hz', unpack=True)

for i in range(len(Data[0])):
	if Data[1, i] > 317.41 - 43.19*3 and Data[1, i] < 317.41 + 43.19*3:   # APD90
		if Data[7, i] > -73.98 - 3.99*3 and Data[7, i] < -73.98 + 3.99*3:   # RMP
			if Data[3, i] > 138.09 - 45.14*3 and Data[3, i] < 138.09 + 45.14*3:   # APD50
				if  Data2[7, i] > -73.98 - 3.99*3 and Data2[7, i] < -73.98 + 3.99*3:   # RMP

					# print i
					sel.append(i)




rejected =  []; #range(300)


# rejected = list(set(sel).difference(set(range(300))))

# rejected = list(set(sel) - set(range(300)))
for i in range(1200):
	if i not in sel:
		rejected.append(i)


print len(sel)
print len(rejected)
print rejected

import numpy as np
import collections


import matplotlib

# matplotlib.use('Qt4Agg') 

import matplotlib.gridspec as gridspec

import matplotlib.pyplot as plt
import csv
from matplotlib import rc
import pandas as pd

## rcParams are the default parameters for matplotlib
import matplotlib as mpl
rc('mathtext',default='regular')
fonts = 20;
fonts_title = 20
leng_fontsize = 17
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

# fig = plt.figure(figsize=(19*0.8/2.54,0.5*19*1.5/2.54))
fig = plt.figure(figsize=(19*0.7*0.8/2.54,0.5*19*1.5/2.54))

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



for i in np.arange(num_col*num_row):
	ax = fig.axes[i]
	# ax.autoscale(enable=True, axis='x')  
	# ax.autoscale(enable=True, axis='r')  
	ax.spines['right'].set_visible(False)
	#ax.spines['left'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(True)
	ax.xaxis.set_ticks_position('bottom')
	ax.spines['bottom'].set_linewidth(2)
	ax.spines['left'].set_linewidth(2)
	ax.spines['right'].set_linewidth(2)
	ax.yaxis.tick_left()
	# ax.set_xlim(308950+30,309350-20)
	ax.set_ylim(-5,400)
	ax.set_xlim(-87,-40)
	# # plt.xticks(np.linspace(0,600, 5))
	# ax.set_xticks([0, 250, 500])
	ax.yaxis.set_tick_params(width=2)
	ax.xaxis.set_tick_params(width=2)
	plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95,
		wspace=0.32, hspace=0.16)
	ax.yaxis.set_label_coords(-0.13, 0.5)
# ax.hist(Data[1,sel],20, rwidth=0.5, color='k')
ax.hist([Data[7, sel], Data[7, rejected]],20, rwidth=0.8, color=['r', 'k'],histtype='bar', stacked=True)
# ax.hist([ Data2[1, sel],Data2[1, rejected]],20, rwidth=0.8, color=['r', 'k'],histtype='bar', stacked=True)
ax.plot([-74-3*4, -74,-74+3*4], [30, 30, 30], '*')

plt.xlabel('RMP (mV)')
plt.ylabel('Count')
plt.savefig('RMP1_New.pdf', dpi=300)
plt.show()

# fig = plt.figure(figsize=(19*0.8/2.54,0.5*19*1.5/2.54))

fig = plt.figure(figsize=(19*0.7*0.8/2.54,0.5*19*1.5/2.54))
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



	
# plt.hist(Data2[7, sel],25)
# ax.hist([Data2[7, rejected], Data2[7, sel]],20, rwidth=0.9, color=['b', 'k'])
ax.hist([Data2[7, sel], Data2[7, rejected]],20, rwidth=0.8, color=['r', 'k'],histtype='bar', stacked=True)

for i in np.arange(num_col*num_row):
	ax = fig.axes[i]
	# ax.autoscale(enable=True, axis='x')  
	# ax.autoscale(enable=True, axis='r')  
	ax.spines['right'].set_visible(False)
	#ax.spines['left'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(True)
	ax.xaxis.set_ticks_position('bottom')
	ax.spines['bottom'].set_linewidth(2)
	ax.spines['left'].set_linewidth(2)
	ax.spines['right'].set_linewidth(2)
	ax.yaxis.tick_left()
	# ax.set_xlim(308950+30,309350-20)
	ax.set_ylim(-5, 355)
	ax.set_xlim(-87,-40)
	# ax.set_ylim(0, 355)
	# # plt.xticks(np.linspace(0,600, 5))
	# ax.set_xticks([0, 250, 500])
	ax.yaxis.set_tick_params(width=2)
	ax.xaxis.set_tick_params(width=2)
	plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95,
		wspace=0.32, hspace=0.16)
	ax.yaxis.set_label_coords(-0.13, 0.5)
	# ax.yaxis.set_label_coords(-0.1, 0.5)
plt.xlabel('RMP (mV)')
plt.ylabel('Count')
plt.savefig('RMP3.pdf', dpi=300)

plt.show()


print np.mean(Data[7, sel]),np.std(Data[7, sel]), 'DATA'
print np.mean(Data2[7, sel]),np.std(Data2[7, sel]), 'Data2'

# print sel, len(sel)