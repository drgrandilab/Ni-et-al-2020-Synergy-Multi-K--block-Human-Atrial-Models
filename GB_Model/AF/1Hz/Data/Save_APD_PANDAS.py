
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


# para_set = np.loadtxt('para_log.dat')


def get_aligned_data(filename,sep=' '):    # to read data if there is any non-accepted data
	Data =  np.genfromtxt(filename, delimiter=sep)

	width = len(Data[0,:])


	struct_data = []  #*600;

	struct_data = -1 * np.ones((600, 33))  # size of results at most
	# for i in range(600):
	# 	struct_data.append([])
	# print len(Data[:,0]), 'length'
	print len(struct_data)
	for i in range(len(Data[:, 0])):
		index = int(Data[i,0] - 1);  # index of model ID starts from 1
		# print index
		struct_data[index] = np.copy(Data[i,:])



	nondata = np.ones((1, width))*np.nan
	# print nondata
	for i in range(600):

															#these are excluded models.
		if struct_data[i][0] == -1 or struct_data[i][0] in [266,335,112,410,522]:
			struct_data[i] = nondata   # assign nondata if entry is empty, indicating that the data is empty!!!
			# print 'nnnn'
			# print i + 1

	return struct_data.transpose()
	# return np.copy(np.array(struct_data))






# def get_difference(ori, changed):
# 	avg_by_cell = 0.5*(changed[APDIndex2]+ changed[1]) - 0.5*(ori[APDIndex2]+ ori[1])
# 	avg_by_cell_percentage = avg_by_cell / (0.5*(ori[APDIndex2]+ ori[1]))
# 	return avg_by_cell, avg_by_cell_percentage


# def check_alternans(Data):
# 	alternans= np.abs(Data[1] - Data[APDIndex2])
# 	j=0
# 	for i in alternans:
# 		if i > 5:
# 			j=j+1
# 	print j
# 	return alternans


from scipy import stats



def get_max(A, B):
	value = []
	for i in range(len(A)):
		value.append(np.max(A[i], B[i]))
	return value


def get_min(A, B):
	value = []
	for i in range(len(A)):
		value.append(np.min(A[i], B[i]))
	return value

def data_ana(Model):
	print Model

	block_data = []


	Data0_1Hz = get_aligned_data('../../1Hz/Data/merged_data.dat.model.'+Model+'.drug_'+str(0))
	print(Model)
	Data0_3Hz = get_aligned_data('../../3Hz/Data/merged_data.dat.model.'+Model+'.drug_'+str(0))  # this is toy play, need to be updated with real life analysis
	# Data0_3Hz = get_aligned_data('merged_data.dat.model.'+Model+'.drug_'+str(0))  # this is toy play, need to be updated with real life analysis



	print Data0_1Hz[2,1]
	# load all data under each conditions 
	for i in range(8):
		# Data = np.loadtxt('Data/AF_model_ID_'+ Model + '.Drug.'+str(i), unpack=True)#, delimiter='\t', dtype=np.float)

		Data = get_aligned_data('merged_data.dat.model.'+Model+'.drug_'+str(i))
		block_data.append(Data)



	Data0 = get_aligned_data('merged_data.dat.model.'+Model+'.drug_'+str(0)) #np.loadtxt('Data/AF_model_ID_'+ Model + '.Drug.0', unpack=True)#, delimiter='\t', dtype=np.float)
	Data1 = get_aligned_data('merged_data.dat.model.'+Model+'.drug_'+str(1)) #np.loadtxt('Data/AF_model_ID_'+ Model + '.Drug.1', unpack=True)
	Data2 = get_aligned_data('merged_data.dat.model.'+Model+'.drug_'+str(2)) #np.loadtxt('Data/AF_model_ID_'+ Model + '.Drug.2', unpack=True)
	Data3 = get_aligned_data('merged_data.dat.model.'+Model+'.drug_'+str(3)) #np.loadtxt('Data/AF_model_ID_'+ Model + '.Drug.3', unpack=True)
	Data4 = get_aligned_data('merged_data.dat.model.'+Model+'.drug_'+str(4)) #np.loadtxt('Data/AF_model_ID_'+ Model + '.Drug.4', unpack=True)
	Data5 = get_aligned_data('merged_data.dat.model.'+Model+'.drug_'+str(5)) #np.loadtxt('Data/AF_model_ID_'+ Model + '.Drug.5', unpack=True)
	Data6 = get_aligned_data('merged_data.dat.model.'+Model+'.drug_'+str(6)) #np.loadtxt('Data/AF_model_ID_'+ Model + '.Drug.6', unpack=True)
	Data7 = get_aligned_data('merged_data.dat.model.'+Model+'.drug_'+str(7)) #np.loadtxt('Data/AF_model_ID_'+ Model + '.Drug.7', unpack=True)



	APDIndex = 2;
	APDIndex2 = 2 + 14;

	#check AF model, if within range
	# remove models outside range
	alternans_cutoff = 5;
	print len(Data0[0])

	print Data0_1Hz[APDIndex,0], Data0_1Hz[APDIndex,1]
	for i in range(len(Data0[0])):

		if not (Data0_1Hz[APDIndex,i] > 112 and Data0_1Hz[APDIndex,i] < 217+3*35.74 and Data0_1Hz[9,i] < -76.85+3*3.61 and Data0_1Hz[9,i] > -76.85-3*3.61 \
		          and Data0_1Hz[6,i] > 100.41-3*27.68 and Data0_1Hz[6, i] < 100.41+3*27.68 \
		          and Data0_3Hz[9,i] < -76.85+3*3.61 and Data0_3Hz[9,i] > -76.85-3*3.61 \
		          and np.abs(Data0_1Hz[APDIndex, i] -Data0_1Hz[APDIndex2, i] )<alternans_cutoff \
				  and np.abs(Data0_3Hz[APDIndex, i] -Data0_3Hz[APDIndex2, i] )< alternans_cutoff) :
			Data0[APDIndex2, i] = np.nan
			Data0[APDIndex, i] = np.nan


	A = (Data0[APDIndex2] + Data0[APDIndex]) /2.0
	# A = train_data[:, 16]
	C = (Data7[APDIndex2]-Data0[APDIndex2] + Data7[APDIndex]-Data0[APDIndex]) /2.0
	D = (Data1[APDIndex2] +Data2[APDIndex2]+Data4[APDIndex2]-3*Data0[APDIndex2] + Data1[APDIndex] +Data2[APDIndex]+Data4[APDIndex]-3*Data0[APDIndex])/2.0

	B = C - D



	# create local copy, because python only assign by reference here ... 

	A = np.copy(A)
	B = np.copy(B)
	C = np.copy(C)
	D = np.copy(D)



	print stats.ttest_ind(C, D,equal_var=False)
	print stats.ttest_rel(C, D), len(C)


	# relative change in APD
	Model_1_prolong = 100.0 * 0.5*(Data1[APDIndex2] - Data0[APDIndex2] + Data1[APDIndex] - Data0[APDIndex] )/ (Data0[APDIndex2] +Data0[APDIndex]) *2.0;
	Model_2_prolong = 100.0 * 0.5* (Data2[APDIndex2] - Data0[APDIndex2] + Data2[APDIndex] - Data0[APDIndex])/ (Data0[APDIndex2] +Data0[APDIndex]) *2.0;
	Model_3_prolong = 100.0 * 0.5* (Data3[APDIndex2] - Data0[APDIndex2] + Data3[APDIndex] - Data0[APDIndex])/ (Data0[APDIndex2] +Data0[APDIndex]) *2.0;
	Model_4_prolong = 100.0 * 0.5* (Data4[APDIndex2] - Data0[APDIndex2] + Data4[APDIndex] - Data0[APDIndex])/ (Data0[APDIndex2] +Data0[APDIndex]) *2.0;
	Model_5_prolong = 100.0 * 0.5* (Data5[APDIndex2] - Data0[APDIndex2] + Data5[APDIndex] - Data0[APDIndex])/ (Data0[APDIndex2] +Data0[APDIndex]) *2.0;
	Model_6_prolong = 100.0 * 0.5* (Data6[APDIndex2] - Data0[APDIndex2] + Data6[APDIndex] - Data0[APDIndex])/ (Data0[APDIndex2] +Data0[APDIndex]) *2.0;
	Model_7_prolong = 100.0 * 0.5* (Data7[APDIndex2] - Data0[APDIndex2] + Data7[APDIndex] - Data0[APDIndex])/ (Data0[APDIndex2] +Data0[APDIndex]) *2.0;


	APD_Copy = [];
	APD_Copy.append(np.copy((Data0[APDIndex2] +Data0[APDIndex]) /2.0))
	APD_Copy.append(np.copy((Data1[APDIndex2] +Data1[APDIndex]) /2.0))
	APD_Copy.append(np.copy((Data2[APDIndex2] +Data2[APDIndex]) /2.0))
	APD_Copy.append(np.copy((Data3[APDIndex2] +Data3[APDIndex]) /2.0))
	APD_Copy.append(np.copy((Data4[APDIndex2] +Data4[APDIndex]) /2.0))
	APD_Copy.append(np.copy((Data5[APDIndex2] +Data5[APDIndex]) /2.0))
	APD_Copy.append(np.copy((Data6[APDIndex2] +Data6[APDIndex]) /2.0))
	APD_Copy.append(np.copy((Data7[APDIndex2] +Data7[APDIndex]) /2.0))

	# cutoff_APD_AF = 112

	cutoff_APD_AF = 112

	for i in range(1,8):
		APD_Copy[i][np.isnan(APD_Copy[0])] = np.nan;

	alternans_cutoff = 5;

	APD_Copy[1][np.abs(Data1[APDIndex] - Data1[APDIndex2] )> 5] = np.nan
	APD_Copy[2][np.abs(Data2[APDIndex] - Data2[APDIndex2] )> 5] = np.nan
	APD_Copy[3][np.abs(Data3[APDIndex] - Data3[APDIndex2] )> 5] = np.nan
	APD_Copy[4][np.abs(Data4[APDIndex] - Data4[APDIndex2] )> 5] = np.nan
	APD_Copy[5][np.abs(Data5[APDIndex] - Data5[APDIndex2] )> 5] = np.nan
	APD_Copy[6][np.abs(Data6[APDIndex] - Data6[APDIndex2] )> 5] = np.nan
	APD_Copy[7][np.abs(Data7[APDIndex] - Data7[APDIndex2] )> 5] = np.nan
	import pandas as pd

	# APD_Copy[0][(Data0[16] +Data0[1]) /2.0<cutoff_APD_AF] = np.nan;

	# get PANDAS data frame here
	# selected_para=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 57, 58, 59, 60, 61, 63, 64, 65, 67, 69, 70, 72, 73, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 86, 87, 88, 90, 91, 92, 93, 94, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 114, 115, 116, 117, 120, 121, 122, 123, 124, 125, 126, 127, 128, 130, 131, 132, 134, 135, 137, 138, 139, 140, 141, 143, 144, 145, 146, 147, 148, 149, 150, 152, 153, 154, 155, 157, 158, 159, 161, 162, 165, 166, 167, 168, 169, 170, 171, 172, 175, 176, 177, 178, 179, 180, 181, 182, 183, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 200, 201, 203, 206, 208, 210, 212, 213, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 242, 243, 244, 245, 247, 248, 249, 250, 251, 252, 253, 254, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 329, 330, 331, 332, 333, 334, 335, 336, 337, 339, 340, 341, 342, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 360, 361, 362, 364, 366, 367, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 384, 385, 386, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 410, 411, 413, 414, 415, 416, 420, 421, 422, 423, 424, 425, 426, 427, 428, 430, 431, 432, 434, 435, 436, 437, 438, 440, 441, 442, 444, 446, 447, 448, 450, 452, 453, 454, 456, 457, 458, 459, 460, 461, 462, 463, 464, 466, 467, 469, 470, 471, 472, 474, 475, 476, 477, 478, 479, 480, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 506, 507, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 521, 522, 523, 524, 525, 526, 527, 529, 531, 533, 534, 535, 537, 538, 539, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 564, 565, 566, 568, 569, 570, 572, 573, 575, 576, 577, 578, 579, 580, 582, 583, 584, 585, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600]
	selected_para = range(600)
	new_data_frame = pd.DataFrame(APD_Copy, columns=selected_para, index = ['Baseline',  'IK2p', 'ISK', 'IK2p_ISK',  'IKur', 'IK2p_IKur', 'ISK_IKur', 'IK2P_ISK_IKur'])

	new_data_frame.to_csv('Summarise_data_Model_'+Model+'.csv', sep=',',na_rep='NAN')

	new_data_frame.to_pickle('Summarise_data_Model_'+Model+'.pkl')


if __name__ == '__main__':
	# file=open('Model_stats.dat.1', 'a+')
	# print >> file, 'Model', 
	# for i in range(8):
	# 	print >> file, 'length', 'mean_%i'% i, 'std_%i'% i,
	# print >> file, '\n'
	# file.close()

	for i in range(12):
		data_ana(str(i))
