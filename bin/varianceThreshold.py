#!/usr/local/bin/python

import numpy as np
import sklearn as skl
import os 
import sys
from sklearn.feature_selection import VarianceThreshold

if len(sys.argv) != 3:
	print('\nUsage: python '+ sys.argv[0] + ' allele_matrix_file variance_threshold_(0.0 ... 0.9)\n')
	print('  python '+sys.argv[0]+' ./pgProfiles/alleleMatrix 0.5\n')
	sys.exit()

allele_matrix_file = sys.argv[1]
Threshold = sys.argv[2]
Threshold = float(Threshold)

#np.set_printoptions(threshold=30000)
myfile = open(allele_matrix_file)
arr = myfile.readlines()
myfile.close()
arr = np.asarray(arr)
arrsplit = []

### Split arr by tab ###
ii = 0
while ii < len(arr):
	arrsplit.append(arr[ii].split('\t'))
	ii += 1
strain_num = len(arrsplit[0])


### Save A1 - A487 Allele ###
### Dict to save all alleles ###
AllAllele = {}
AllAllele_selected = {}
for x in range(0,strain_num): # x = 0~487
	AllAllele['A'+ repr(x)] = []

### Append the column data to dict
for idx, val in enumerate(arrsplit):
	for strainIdx in range(0,strain_num):
		AllAllele['A'+repr(strainIdx)].append(arrsplit[idx][strainIdx])

### Change data to the form for feature selection ###
integrate = []
for xx in range(1,strain_num):
	integrate.append(AllAllele['A'+repr(xx)])

### Use Variance threshold to select data ###
sele = VarianceThreshold(Threshold*(1-Threshold))#threshold=(.3 * (1 - .3))
selout = sele.fit_transform(integrate) #Get output
seloutsup = sele.get_support(indices=True) #Get index

print ('variance_'+str(Threshold))
for i in seloutsup:
        print (AllAllele['A0'][i])

