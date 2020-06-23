

from subprocess import call
import math
# S1 = 500
# S2 = 250
import sys


import numpy as np

BCL=333.333;



def run_Normal_For_IC(Mode, Para, AF_model, block_para, Drug_ID, f):
	# global mut
	#call(["./main","BCL", str(S1), "S2", str(S2), "Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./SAN', str(P1), str(P2)]) #"Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./NCZ_Model_New_Ito', '250', '14', '250', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	# call(['./NCZ_Model_New_Ito', '500', '14', '500', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	P=['./NCZ_Model_for_IC', str(BCL), '14', str(BCL), 'WT', Mode, '0', '1']
	for i in range(len(Para)):
		P.append(str(Para[i]))
	for i in range(len(AF_model)):
		P.append(str(AF_model[i]))
	for i in range(len(block_para)):
		P.append(str(block_para[i]))

	P.append(str(Drug_ID));
	call(P, stdout=f)



def run_Normal(Para, AF_model, block_para,f,Ini_cond_File, ID=0):
	# global mut
	#call(["./main","BCL", str(S1), "S2", str(S2), "Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./SAN', str(P1), str(P2)]) #"Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./NCZ_Model_New_Ito', '250', '14', '250', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	# call(['./NCZ_Model_New_Ito', '500', '14', '500', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	P=['./ONE_D_NCZ_VEC_Ini_Cond', 'BCL', str(BCL), 'S2', str(BCL-100), 'S1_number', str(100+int(1010000/BCL)-1), 'AF', '4', 'Popl_SF_Number', '15', 'Initial_conditions', Ini_cond_File, 'Time_Start', str(1010000), 'Total_time', str(1010000+100*BCL) ]
	P.append('Popul_scaling_factors')
	for i in range(len(Para)):
		P.append(str(Para[i]))

	P.append('Remodelling_F_Number')
	P.append('3')
	P.append('Remodelling_Factors')

	for i in range(len(AF_model)):
		P.append(str(AF_model[i]))

	P.append('Modulation_F_Number')
	P.append('3')
	P.append('Modulation_Factors')
	for i in range(len(block_para)):
		P.append(str(block_para[i]))


	P.append('OneD_OutFile')
	P.append(f)
	P.append('Sim_ID')
	P.append(str(ID))
	print P
	call(P)





def run_VW(Para, AF_model, block_para,f,Ini_cond_File, S2, ID=0):
	# global mut
	#call(["./main","BCL", str(S1), "S2", str(S2), "Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./SAN', str(P1), str(P2)]) #"Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./NCZ_Model_New_Ito', '250', '14', '250', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	# call(['./NCZ_Model_New_Ito', '500', '14', '500', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	P=['./ONE_D_NCZ_VEC', 'BCL', str(BCL), 'S2', str(S2), 'S1_number', str(100+int(1010000/BCL)-1), 'AF', '4', 'Popl_SF_Number', '15', 'Initial_conditions', Ini_cond_File, 'Time_Start', str( 1010000+100*BCL -2*BCL - 10 ), 'Total_time', str(1010000+100*BCL) ]
	P.append('Popul_scaling_factors')
	for i in range(len(Para)):
		P.append(str(Para[i]))

	P.append('Remodelling_F_Number')
	P.append('3')
	P.append('Remodelling_Factors')

	for i in range(len(AF_model)):
		P.append(str(AF_model[i]))

	P.append('Modulation_F_Number')
	P.append('3')
	P.append('Modulation_Factors')
	for i in range(len(block_para)):
		P.append(str(block_para[i]))


	P.append('OneD_OutFile')
	P.append(f)
	P.append('Sim_ID')
	P.append(str(ID))
	print P
	call(P)
	with open('conduction.log',"r") as f:
		line1 = f.readline().strip();
		line2 = f.readline().strip();
	if int(line1) not in [0,1,2] or int(line2) not in [0,1,2]:
		print sys.stderr, '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> error in reading ''conduction.log'', ID = ', ID, 'line1 = ', line1, 'line2 = ', line2
		sys.exit(0);


	return int(line1), int(line2);

# def run(S2):
# 	call(["./ONE_D_CNZ_C_VEC","BCL", str(S1), "S2", str(S2), "Mutation", "D322H", "S1_number", "20"])  #
# 	with open('conduction.log',"r") as f:
# 		line = f.readline().strip();
# 		if line != "2":
# 			return False;
# 		line = f.readline().strip();
# 		if line != "2":
# 			return False;
# 	return True


def check_if_bilateral_cond((a,b)):
	if (a==2 and b == 2):
		go=True
	else:
		go=False
	return go

def check_if_bilateral_block((a,b)):
	if (a==1 and b == 1):
		go=True
	else:
		go=False
	return go



def run_SimDrug(Para, AF_model, block_para,f):
	# global mut
	#call(["./main","BCL", str(S1), "S2", str(S2), "Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./SAN', str(P1), str(P2)]) #"Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./NCZ_Model_New_Ito', '250', '14', '250', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	# call(['./NCZ_Model_New_Ito', '500', '14', '500', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	P=['./NCZ_Model', '1000', '14', '1000', 'WT', 'SimDrug', '0', '1']
	for i in range(len(Para)):
		P.append(str(Para[i]))
	for i in range(len(AF_model)):
		P.append(str(AF_model[i]))
	for i in range(len(block_para)):
		P.append(str(block_para[i]))
	call(P, stdout=f)

def get_upper_bound(para_set, AF_set, block_set,OneD_files,Cell_ID, S2):

	f = open('workfile_upper', 'a')
	
	go = check_if_bilateral_cond(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))
	lastgo = go
	# go = reduce_S2(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))

	# run_VW(Para, AF_model, block_para,f,Ini_cond_File, S2, ID=0)
	# call('mv OneD_output.dat S2.BCL.3Hz.Popul.Index.'+str(Cell_ID)+'.AF_Model.' + str(AF_ID) + '.block.'+str(Drug_ID), shell=True);
	compute = True
	d_S2 = 30;
	while compute:
			if(lastgo == (not go) ):
				d_S2=d_S2/2.0

			if go:
				S2 =S2-d_S2 +0.001
			else:
				S2 =S2+d_S2 -0.001
			lastgo=go
			go=check_if_bilateral_cond(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))
			print  >>  f, "S1= ", BCL, "S2= ", S2, "d_S2 = ", d_S2, "go = ", go, ' Cell_ID = ', Cell_ID
			if d_S2<0.05 and go:
				compute = False
			if S2> BCL or S2 < 50:
				compute= False
				S2 = np.nan
			if np.isnan(S2):
				compute=False
	return S2,d_S2

def get_lowwer_bound(para_set, AF_set, block_set,OneD_files,Cell_ID, S2):

	f = open('workfile_lowwer', 'a')
	
	go = check_if_bilateral_block(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))
	lastgo = go
	# go = reduce_S2(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))

	# run_VW(Para, AF_model, block_para,f,Ini_cond_File, S2, ID=0)
	# call('mv OneD_output.dat S2.BCL.3Hz.Popul.Index.'+str(Cell_ID)+'.AF_Model.' + str(AF_ID) + '.block.'+str(Drug_ID), shell=True);
	compute = True
	d_S2 = 5;
	while compute:
			if(lastgo == (not go) ):
				d_S2=d_S2/2.0
			if go:
				S2 =S2+d_S2 -0.001
			else:
				S2 =S2-d_S2 +0.001
			lastgo=go
			go=check_if_bilateral_block(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))
			print  >>  f, "S1= ", BCL, "S2= ", S2, "d_S2 = ", d_S2, "go = ", go, ' Cell_ID = ', Cell_ID
			if d_S2<0.05 and go:
				compute = False
			if S2> BCL or S2 < 50:
				compute= False
				S2 = np.nan
			if np.isnan(S2):
				compute=False
	return S2, d_S2



#AF_Model 1, 3Hz
#758

#AF_Model 5, 3Hz
#769


# if __name__ == '__main__':
# 	main()

selected_para = [0, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 63, 64, 65, 66, 67, 68, 69, 71, 72, 73, 74, 75, 76, 77, 79, 80, 81, 82, 83, 84, 85, 86, 87, 89, 90, 91, 92, 94, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 145, 147, 148, 149, 150, 151, 152, 153, 154, 156, 157, 159, 160, 161, 162, 164, 165, 166, 167, 169, 171, 172, 173, 175, 176, 177, 179, 180, 181, 182, 184, 185, 186, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 202, 203, 204, 205, 206, 209, 210, 211, 212, 213, 214, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 235, 236, 238, 239, 240, 241, 242, 243, 244, 245, 247, 248, 249, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 262, 263, 264, 265, 266, 267, 269, 270, 271, 272, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 291, 292, 293, 294, 295, 296, 297, 298, 300, 301, 302, 303, 304, 305, 306, 308, 309, 310, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 338, 339, 344, 345, 346, 348, 349, 350, 351, 352, 353, 355, 356, 357, 358, 359, 361, 363, 365, 366, 367, 368, 369, 370, 371, 372, 374, 376, 377, 378, 379, 380, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 397, 398, 399, 400, 401, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 423, 424, 425, 426, 427, 428, 429, 431, 432, 434, 435, 436, 437, 438, 439, 440, 441, 443, 444, 445, 446, 447, 448, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 461, 463, 464, 466, 467, 468, 469, 471, 474, 475, 476, 477, 479, 480, 481, 483, 484, 485, 487, 488, 489, 490, 491, 493, 494, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 510, 511, 512, 513, 514, 516, 517, 519, 520, 521, 522, 523, 524, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 537, 538, 539, 540, 541, 542, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 560, 561, 562, 563, 564, 565, 566, 568, 569, 570, 571, 572, 573, 574, 575, 576, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 595, 596, 597, 600, 601, 602, 604, 605, 606, 607, 608, 609, 610, 612, 613, 614, 616, 618, 619, 620, 621, 622, 623, 624, 626, 627, 629, 631, 632, 633, 634, 635, 636, 637, 638, 641, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 659, 660, 661, 663, 664, 665, 667, 669, 670, 672, 673, 674, 675, 676, 677, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 690, 691, 692, 694, 695, 696, 697, 698, 700, 701, 702, 704, 705, 706, 708, 709, 710, 711, 712, 713, 715, 716, 718, 719, 722, 723, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 736, 737, 738, 742, 743, 744, 745, 746, 747, 749, 751, 752, 753, 754, 757, 759, 760, 761, 762, 763, 764, 765, 766, 767, 768, 769, 770, 771, 773, 775, 776, 777, 778, 779, 780, 782, 783, 784, 786, 787, 788, 789, 790, 791, 792, 793, 794, 795, 796, 797, 798, 800, 801, 802, 803, 804, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 818, 819, 820, 821, 822, 823, 824, 825, 826, 827, 828, 829, 830, 831, 833, 834, 837, 838, 840, 841, 842, 843, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874, 875, 877, 878, 879, 880, 881, 882, 884, 885, 887, 889, 890, 891, 892, 893, 894, 895, 896, 898, 900, 901, 902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921, 922, 924, 927, 928, 929, 931, 932, 934, 935, 936, 937, 938, 939, 940, 941, 942, 944, 945, 946, 947, 949, 950, 951, 952, 953, 954, 955, 956, 957, 959, 960, 961, 962, 963, 965, 966, 967, 968, 969, 970, 971, 972, 973, 974, 975, 976, 977, 978, 979, 980, 981, 982, 985, 986, 987, 988, 989, 990, 991, 992, 993, 994, 995, 996, 997, 998, 999, 1000, 1001, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1026, 1027, 1029, 1030, 1031, 1032, 1034, 1035, 1036, 1037, 1039, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 1074, 1075, 1076, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1089, 1091, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1099, 1100, 1103, 1104, 1105, 1106, 1108, 1109, 1110, 1111, 1112, 1115, 1116, 1118, 1119, 1120, 1122, 1123, 1124, 1126, 1127, 1129, 1130, 1131, 1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1141, 1142, 1144, 1145, 1146, 1147, 1148, 1149, 1150, 1152, 1153, 1154, 1155, 1158, 1159, 1160, 1161, 1162, 1163, 1164, 1165, 1167, 1168, 1169, 1170, 1171, 1172, 1173, 1174, 1175, 1176, 1177, 1178, 1180, 1181, 1182, 1183, 1184, 1185, 1186, 1187, 1188, 1189, 1190, 1191, 1192, 1193, 1195, 1196, 1197, 1198];

AF_ID = int(sys.argv[1])
if (AF_ID == 1):
	selected_para = [0, 3, 6, 7, 8, 10, 13, 14, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 33, 34, 35, 36, 37, 38, 42, 44, 45, 46, 47, 48, 50, 53, 54, 57, 58, 59, 60, 63, 64, 65, 69, 71, 72, 73, 75, 76, 79, 80, 81, 82, 83, 85, 86, 87, 89, 90, 91, 94, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 117, 120, 122, 125, 126, 127, 130, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 145, 147, 150, 151, 153, 154, 156, 159, 160, 162, 164, 165, 166, 169, 171, 172, 173, 175, 176, 177, 179, 180, 184, 186, 189, 191, 192, 193, 196, 197, 198, 199, 203, 205, 206, 209, 212, 214, 216, 219, 220, 221, 222, 223, 224, 225, 226, 227, 229, 230, 231, 233, 235, 239, 240, 241, 242, 243, 244, 245, 248, 249, 251, 252, 253, 254, 255, 256, 259, 260, 263, 264, 265, 267, 269, 270, 271, 276, 277, 278, 279, 280, 281, 282, 283, 284, 286, 287, 288, 289, 292, 293, 294, 295, 296, 297, 301, 303, 305, 309, 310, 312, 313, 314, 315, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 331, 332, 333, 335, 336, 339, 345, 346, 348, 350, 351, 353, 355, 356, 357, 358, 361, 363, 366, 367, 368, 369, 370, 371, 372, 374, 377, 378, 379, 380, 383, 384, 385, 386, 387, 389, 390, 392, 394, 397, 398, 399, 400, 401, 403, 405, 407, 408, 409, 410, 411, 412, 413, 415, 416, 418, 419, 420, 423, 424, 426, 427, 428, 429, 431, 432, 434, 435, 437, 438, 439, 441, 443, 444, 445, 447, 448, 451, 452, 454, 455, 457, 458, 459, 463, 464, 466, 467, 468, 469, 474, 475, 477, 479, 480, 481, 483, 484, 485, 487, 488, 489, 490, 491, 493, 498, 499, 500, 502, 504, 506, 507, 508, 511, 513, 516, 519, 521, 522, 523, 524, 526, 527, 528, 529, 530, 531, 533, 535, 537, 538, 540, 542, 545, 546, 547, 548, 549, 551, 554, 557, 560, 561, 563, 564, 565, 568, 570, 571, 572, 573, 574, 575, 581, 582, 584, 585, 587, 588, 590, 592, 593, 596, 597, 600, 601, 602, 604, 606, 608, 609, 610, 612, 613, 618, 619, 620, 621, 622, 623, 624, 626, 627, 629, 632, 633, 634, 635, 636, 637, 638, 645, 646, 649, 651, 652, 653, 654, 655, 656, 657, 659, 663, 664, 665, 667, 669, 670, 672, 673, 674, 675, 676, 677, 679, 680, 681, 682, 684, 685, 686, 687, 688, 691, 692, 694, 695, 696, 697, 698, 702, 704, 706, 708, 710, 711, 713, 715, 718, 719, 722, 723, 725, 726, 727, 728, 734, 736, 737, 738, 742, 743, 744, 745, 746, 747, 752, 753, 754, 759, 761, 765, 766, 767, 768, 769, 770, 771, 773, 775, 776, 777, 778, 779, 782, 783, 786, 791, 792, 795, 798, 800, 801, 802, 803, 804, 806, 807, 808, 810, 811, 812, 813, 815, 816, 818, 819, 820, 822, 823, 824, 825, 826, 827, 828, 829, 830, 833, 834, 837, 842, 843, 849, 850, 852, 853, 857, 858, 860, 861, 862, 865, 866, 867, 868, 869, 870, 871, 873, 874, 875, 877, 879, 880, 881, 882, 884, 885, 887, 889, 893, 894, 898, 900, 901, 904, 907, 908, 909, 910, 912, 913, 914, 915, 916, 918, 919, 920, 921, 922, 924, 928, 929, 931, 932, 934, 935, 936, 937, 938, 939, 940, 941, 944, 946, 947, 949, 951, 952, 953, 954, 955, 956, 959, 960, 961, 963, 965, 967, 968, 969, 970, 971, 972, 973, 974, 975, 976, 977, 978, 979, 981, 982, 986, 987, 988, 989, 990, 992, 993, 994, 995, 996, 997, 998, 999, 1001, 1003, 1004, 1005, 1006, 1008, 1009, 1010, 1011, 1014, 1015, 1016, 1017, 1018, 1019, 1021, 1022, 1027, 1029, 1032, 1039, 1041, 1042, 1043, 1045, 1046, 1047, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1063, 1064, 1065, 1066, 1067, 1069, 1070, 1071, 1072, 1073, 1076, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1087, 1089, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1103, 1105, 1106, 1108, 1110, 1111, 1112, 1116, 1118, 1120, 1122, 1124, 1126, 1127, 1129, 1130, 1133, 1134, 1135, 1136, 1137, 1139, 1141, 1144, 1145, 1146, 1147, 1149, 1150, 1153, 1154, 1155, 1158, 1160, 1162, 1163, 1164, 1165, 1168, 1171, 1172, 1174, 1176, 1177, 1178, 1181, 1183, 1184, 1186, 1187, 1188, 1189, 1190, 1191, 1192, 1193, 1195, 1196, 1198] 
#AF_Model 5, 1Hz
elif (AF_ID == 5):
	selected_para = [0, 3, 6, 7, 8, 10, 13, 14, 15, 19, 20, 21, 22, 23, 24, 25, 26, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 42, 45, 47, 48, 50, 52, 53, 57, 58, 59, 60, 61, 63, 64, 65, 66, 69, 71, 72, 73, 75, 79, 80, 81, 82, 83, 85, 86, 87, 89, 90, 91, 92, 94, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 117, 118, 120, 121, 122, 123, 125, 126, 127, 129, 130, 131, 132, 133, 134, 135, 137, 138, 139, 140, 141, 142, 143, 145, 147, 150, 151, 154, 156, 159, 160, 162, 164, 165, 169, 171, 172, 173, 175, 176, 177, 180, 182, 184, 186, 189, 191, 192, 193, 194, 196, 197, 198, 199, 203, 204, 205, 206, 209, 212, 214, 216, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 233, 235, 239, 240, 241, 242, 243, 244, 248, 249, 251, 252, 253, 254, 255, 256, 259, 260, 262, 264, 269, 270, 271, 272, 276, 277, 278, 279, 280, 281, 283, 284, 287, 288, 291, 292, 293, 294, 295, 296, 297, 298, 301, 302, 303, 305, 309, 310, 312, 313, 314, 315, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 328, 331, 332, 335, 336, 339, 345, 346, 348, 350, 351, 353, 355, 356, 357, 358, 361, 363, 366, 367, 368, 369, 370, 371, 372, 374, 376, 377, 378, 379, 380, 383, 384, 385, 386, 387, 389, 390, 392, 394, 397, 398, 399, 400, 401, 403, 404, 405, 406, 407, 408, 409, 410, 412, 413, 414, 415, 416, 418, 419, 420, 423, 424, 426, 428, 429, 431, 432, 434, 435, 437, 439, 441, 443, 445, 447, 448, 451, 452, 454, 455, 457, 458, 459, 463, 466, 467, 468, 469, 474, 475, 476, 477, 479, 480, 481, 483, 484, 485, 487, 488, 489, 490, 491, 493, 498, 499, 500, 501, 504, 507, 508, 511, 512, 513, 516, 519, 521, 522, 523, 524, 526, 527, 528, 529, 531, 532, 533, 535, 537, 539, 540, 541, 542, 544, 546, 547, 548, 550, 551, 557, 558, 560, 561, 562, 563, 564, 565, 566, 568, 569, 570, 571, 572, 573, 574, 575, 581, 584, 585, 587, 588, 592, 593, 596, 597, 600, 601, 602, 604, 606, 607, 608, 609, 610, 612, 613, 618, 619, 620, 621, 622, 623, 624, 626, 627, 629, 632, 633, 634, 635, 636, 638, 645, 646, 648, 649, 650, 651, 652, 653, 654, 655, 657, 659, 661, 663, 664, 665, 667, 669, 670, 672, 673, 674, 675, 676, 677, 679, 680, 681, 682, 684, 685, 686, 687, 688, 690, 691, 692, 694, 695, 696, 697, 698, 702, 704, 706, 708, 711, 713, 715, 718, 719, 722, 723, 725, 726, 727, 728, 729, 734, 736, 737, 738, 742, 743, 744, 745, 746, 747, 749, 753, 759, 761, 764, 765, 766, 767, 768, 769, 770, 773, 775, 776, 777, 778, 779, 782, 783, 786, 789, 790, 791, 792, 795, 796, 798, 800, 801, 802, 803, 804, 806, 808, 811, 812, 813, 815, 816, 818, 819, 820, 821, 822, 823, 824, 825, 826, 827, 828, 829, 830, 831, 833, 834, 837, 838, 840, 841, 842, 843, 846, 847, 849, 850, 852, 853, 857, 858, 860, 862, 865, 866, 867, 868, 869, 870, 871, 873, 874, 875, 877, 879, 880, 881, 882, 884, 885, 887, 889, 893, 894, 898, 901, 902, 905, 907, 908, 909, 910, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921, 922, 924, 928, 929, 931, 932, 934, 935, 936, 938, 939, 942, 944, 946, 947, 949, 951, 952, 954, 955, 956, 957, 959, 960, 961, 963, 965, 966, 967, 968, 969, 970, 971, 972, 973, 974, 975, 976, 977, 978, 979, 981, 982, 986, 988, 990, 991, 992, 993, 994, 995, 996, 997, 998, 999, 1000, 1001, 1003, 1004, 1005, 1006, 1008, 1009, 1010, 1011, 1014, 1015, 1016, 1017, 1018, 1019, 1022, 1023, 1027, 1029, 1030, 1032, 1034, 1039, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 1074, 1076, 1081, 1082, 1083, 1084, 1085, 1087, 1089, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1103, 1106, 1108, 1110, 1111, 1112, 1115, 1116, 1118, 1119, 1120, 1122, 1124, 1127, 1129, 1130, 1133, 1134, 1135, 1136, 1137, 1139, 1141, 1142, 1144, 1145, 1146, 1147, 1149, 1150, 1153, 1154, 1155, 1158, 1160, 1162, 1163, 1164, 1169, 1171, 1172, 1176, 1177, 1181, 1183, 1186, 1187, 1188, 1189, 1190, 1191, 1192, 1193, 1195, 1196, 1198] 
else:
	print sys.stderr, '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> input AF_ID error, please check'
	sys.exit(0);

call(['make'])
# para_set = np.loadtxt('../pre_optim.dat')
# para_set = np.loadtxt('../../para_log.dat')
# AF_set = np.loadtxt('../../AF_models.dat')
# block_set = np.loadtxt('../../block_para.dat')
# call(['rm', 'log.dat'])
para_set = np.loadtxt('para_log.dat')
AF_set = np.loadtxt('AF_models.dat')
block_set = np.loadtxt('block_para.dat')

files=[]
singlecell_file=[]
OneD_files=[]
for i in range(8):
	singlecell_file.append(open("AP.log.dat.AF_Model." + sys.argv[1] + '.block.'+str(i), "w+"))
	files.append("AP.log.dat.AF_Model." + sys.argv[1] + '.block.'+str(i))
	OneD_files.append("OneD.log.dat.AF_Model." + sys.argv[1] + '.block.'+str(i))
	call('mkdir Block_Model_'+str(i), shell=True)
# f2 = open("para_log.dat", "w+")
singlecell_file_pre = open("Singlecell.AP.log.dat.AF_Model." + sys.argv[1] + '.block.'+str(-1), "w+")
vw_results=[]
for i in range(8):
	vw_results.append(open("vw.results..AF_Model." + sys.argv[1] + '.block.'+str(i), "w+"))
# sys.stdout = open('file', 'w')sys
# AF_ID = int(sys.argv[1])
# num2 = int(sys.argv[2])




# Cell_ID =3;
# S2 = 179.9999999+0.001
# Drug_ID = 0;
# run_VW(para_set[Cell_ID], AF_set[AF_ID], block_set[Drug_ID],OneD_files[Drug_ID], 'Specific', S2, Cell_ID )
# # S2_upper, d_S2_up = get_upper_bound(para_set[Cell_ID], AF_set[AF_ID], block_set[Drug_ID],OneD_files[Drug_ID],Cell_ID, S2)


for ii  in range(len(selected_para)):
	Cell_ID = selected_para[ii];  # ID 
	'''
	run_Normal_For_IC('Normal', para_set[Cell_ID], AF_set[AF_ID], block_set[0], 0, singlecell_file_pre)  # run initial conditions 
	call('cp Restart_ICs/ICs.bin Restart_ICs/AF_ID_%d_CellID_%d_ICs.bin' % (AF_ID, Cell_ID), shell=True);

	for Drug_ID in range(8):
		run_Normal_For_IC('SimDrug', para_set[Cell_ID], AF_set[AF_ID], block_set[Drug_ID], Drug_ID, singlecell_file[Drug_ID])  # run initial conditions 
		call('cp SingleCell_Restart_ICs/ICs_Drug_%d.bin SingleCell_Restart_ICs/AF_ID_%d_CellID_%d_ICs_Drug_%d.bin' % (Drug_ID, AF_ID, Cell_ID, Drug_ID), shell=True);

	'''

	for Drug_ID in range(8):

		Ini_cond_File = 'SingleCell_Restart_ICs/AF_ID_%d_CellID_%d_ICs_Drug_%d.bin' % (AF_ID, Cell_ID, Drug_ID)
		run_Normal(para_set[Cell_ID], AF_set[AF_ID], block_set[Drug_ID],OneD_files[Drug_ID], Ini_cond_File, Cell_ID)
		# call('mv OneD_output.dat S1.BCL.3Hz.Popul.Index.'+str(Cell_ID)+'.AF_Model.' + str(AF_ID) + '.block.'+str(Drug_ID), shell=True);
		S2 = 200
		# # go = True
		
		S2_upper, d_S2_up = get_upper_bound(para_set[Cell_ID], AF_set[AF_ID], block_set[Drug_ID],OneD_files[Drug_ID],Cell_ID, S2)

		# S2 = 200

		S2_lowwer, d_S2_low = get_lowwer_bound(para_set[Cell_ID], AF_set[AF_ID], block_set[Drug_ID],OneD_files[Drug_ID],Cell_ID, S2_upper-10)

		print >> vw_results[Drug_ID], Cell_ID, S2_upper,d_S2_up, S2_lowwer, d_S2_low
		# S2=131.98359375 +0.2

		# run_VW(para_set[Cell_ID], AF_set[AF_ID], block_set[Drug_ID],OneD_files[Drug_ID], 'Specific', S2, Cell_ID )

		# f = open('workfile', 'w')
		
		# go = reduce_S2(run_VW(para_set[Cell_ID], AF_set[AF_ID], block_set[Drug_ID],OneD_files[Drug_ID], 'Specific', S2, Cell_ID ))
		# lastgo = go
		# go = reduce_S2(run_VW(para_set[Cell_ID], AF_set[AF_ID], block_set[Drug_ID],OneD_files[Drug_ID], 'Specific', S2, Cell_ID ))

		# run_VW(Para, AF_model, block_para,f,Ini_cond_File, S2, ID=0)
		# call('mv OneD_output.dat S2.BCL.3Hz.Popul.Index.'+str(Cell_ID)+'.AF_Model.' + str(AF_ID) + '.block.'+str(Drug_ID), shell=True);
		# compute = False
		# d_S2 = 50;
		# while compute:
		# 		if(lastgo == (not go) ):
		# 			d_S2=d_S2/2.0

		# 		if go:
		# 			S2 =S2-d_S2-0.5
		# 		else:
		# 			S2 =S2+d_S2+0.5
		# 		lastgo=go
		# 		go=check_if_bilateral_cond(run_VW(para_set[Cell_ID], AF_set[AF_ID], block_set[Drug_ID],OneD_files[Drug_ID], 'Specific', S2, Cell_ID ))
		# 		print  >>  f, "S1= ", BCL, "S2= ", S2, "d_S2 = ", d_S2, "go = ", go, ' Cell_ID = ', Cell_ID
		# 		if math.fabs(d_S2)<0.5:
		# 			compute = False



		# call('mv AP.BCL.3Hz.Popul.Index.'+index+'.AF_Model.' + sys.argv[1] + '.block.'+str(j) + ' Block_Model_'+str(j) + '/', shell=True)


	
	# i = selected_para[ii];
	# index = str(i)
	# run_Normal(para_set[Cell_ID], AF_set[AF_ID], block_set[0],files[0],i)
	# call('mv OneD_output.dat AP.BCL.3Hz.Popul.Index.'+index+'.AF_Model.' + sys.argv[1] + '.block.'+str(0), shell=True);
	# call('mv AP.BCL.3Hz.Popul.Index.'+index+'.AF_Model.' + sys.argv[1] + '.block.'+str(0) + ' Block_Model_'+str(0) + '/', shell=True)
	
	# for j in range(1,8):
	# 	run_Normal(para_set[Cell_ID], AF_set[AF_ID], block_set[j],files[j], i)
	# 	call('mv OneD_output.dat AP.BCL.3Hz.Popul.Index.'+index+'.AF_Model.' + sys.argv[1] + '.block.'+str(j), shell=True);
	# 	call('mv AP.BCL.3Hz.Popul.Index.'+index+'.AF_Model.' + sys.argv[1] + '.block.'+str(j) + ' Block_Model_'+str(j) + '/', shell=True)
# call(['mv', 'log.dat', 'IKurB_1Hz.log.dat'])
for i in range(8):
	# files[i].close()
	singlecell_file[i].close()

singlecell_file_pre.close()
