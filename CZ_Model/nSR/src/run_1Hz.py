

from subprocess import call
import math
# S1 = 500
# S2 = 250
import sys


import numpy as np


def run(Para):
	# global mut
	#call(["./main","BCL", str(S1), "S2", str(S2), "Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./SAN', str(P1), str(P2)]) #"Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./NCZ_Model_New_Ito', '250', '14', '250', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	# call(['./NCZ_Model_New_Ito', '500', '14', '500', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	P=['./NCZ_Model', '1000', '14', '1000', 'WT', 'Normal', '0', '0']
	for i in range(len(Para)):
		P.append(str(Para[i]))

	call(P, stdout=f)



# run(1000, "SimDrug",100,100, 0.1,0.01)


call(['make'])
# para_set = np.loadtxt('../pre_optim.dat')
para_set = np.loadtxt('para_log.dat')
# call(['rm', 'log.dat'])


f = open("AP.log.dat.1Hz", "w+")
# f2 = open("para_log.dat", "w+")

BCL = 1000
# run(BCL,"Normal", 0,0,0,0);
# Mode="SimDrug"

# sys.stdout = open('file', 'w')

for i in range(len(para_set)):

	# # with open("log.dat", "a+") as f:
	# P1 = para_set[i,1]
	# P2 = para_set[i,2]
	# P3 = para_set[i,3]
	# P4 = para_set[i,4]
	# P5 = para_set[i,5]
	# P6 = para_set[i,6]
	# print >> f2, "%.0f " % para_set[i,0], P1, P2, P3, P4, P5, P6
	# print >> f2, "%.0f " % para_set[i,0], P1, P2, P3, P4, P5, P6
	# print >> f2, "%.0f " % para_set[i,0], P1, P2, P3, P4, P5, P6
	# print para_set[i,0], (P1, P2, P3, P4, P5, P6)
	# run(BCL,Mode, P1, P2, P3, P4, P5, P6);
	# call(['mv', 'Drug_CNZ_out.dat',"%.0f" % para_set[i,0]+'.dat'])
	run(para_set[i])
	index = str(i)
	call('mv Drug_CNZ_out.dat AP.BCL.1000.Popul.Index.'+index, shell=True);
# call(['mv', 'log.dat', 'IKurB_4Hz.log.dat'])
