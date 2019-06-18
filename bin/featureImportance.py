#!/usr/local/bin/python

import numpy as np
import sklearn as skl
from sklearn.ensemble import ExtraTreesClassifier      # For Feature Importance
import sys
import os

if len(sys.argv) != 3:
  print('\nUsage: python '+sys.argv[0]+' in/out_path'+' top_threshold_(1,2,3,4,5)\n')
  print('  python '+sys.argv[0]+' ./topSchemeSel 5\n')
  sys.exit()

in_path = sys.argv[1]
top_threshold = sys.argv[2]

Group_num = 1
nowpath = in_path
status = True

while status == True:
	pathlocal = nowpath + "/alleleMatrix_grouped/alleleMatrix_grouped" + repr(Group_num)
	if os.path.exists(pathlocal):
		myfile = open(pathlocal)
		arr = myfile.readlines()
		myfile.close()
		arr = np.asarray(arr)
		arrsplit = []
		### Split arr by tab ###
		ii = 0
		while ii < len(arr):
			arrsplit.append(arr[ii].split('\t'))
			ii += 1
		alleleamount = len(arrsplit[0])

		### Save A~ - A~ Allele ###
		### Dict to save all alleles ###
		AllAllele = {}
		AllAllele_selected = {}
		for x in range(0,alleleamount):
			AllAllele['A'+ repr(x)] = []
			AllAllele_selected['A'+ repr(x)] = []

		### Append the column data to dict
		for idx, val in enumerate(arrsplit):
			for allelenum in range(0,alleleamount):
				AllAllele['A'+repr(allelenum)].append(arrsplit[idx][allelenum])

		### Change data to the form for feature selection ###
		integrate = []
		for xx in range(1,alleleamount):
			integrate.append(AllAllele['A'+repr(xx)])

		# Convert 'str' to 'int'
		for idx,value in enumerate(integrate):
			for iidx,iivalue in enumerate(integrate[idx]):
				integrate[idx][int(iidx)] = int(integrate[idx][int(iidx)])

		SAMPLES = []
		for i in range(1,alleleamount):
			SAMPLES.append(i)

		locusend = len(arr)

		classifacation = []
		for i in range(0,alleleamount-1):
			classifacation.append(integrate[i][locusend-1])

		trial1 = np.array(SAMPLES)
		trial2 = np.array(integrate)
		trial3 = np.array(classifacation)

		X = trial2[:,0:locusend-1]
		Y = trial3

		## Feature importance 
		# feature extraction
		estimator = locusend-1
		FIselected = []
		model = ExtraTreesClassifier(n_estimators=estimator)#n_estimators=26980,random_state=0,max_features=10000
		model.fit(X, Y)
		importances = model.feature_importances_

		indices = np.argsort(importances)[::-1]
		shapeX = X.shape[1]

		##### Top x feature output #####
		top_userdef_feature = []
		u = int(top_threshold)
		for f in range(0,u):
		    top_userdef_feature.append(indices[f])
		    top_userdef_feature.sort()

		##### OUTPUT FILE #####
		### Ready to output file ###
		myfile = open(in_path+'/alleleMatrix_new')
		arr = myfile.readlines()
		myfile.close()
		arrsplit = []
		### Split arr by tab ###
		ii = 0
		while ii < len(arr):
			arrsplit.append(arr[ii].split('\t'))
			ii += 1
		alleleamount2 = len(arrsplit[0])
		### Save A1 - A34 Allele ###
		AllAllele = {}
		AllAllele_selected = {}
		for x in range(0,alleleamount2):
			AllAllele['A'+ repr(x)] = []

		### Append the column data to dict
		for idx, val in enumerate(arrsplit):
			for allelenum in range(0,alleleamount2):
				AllAllele['A'+repr(allelenum)].append(arrsplit[idx][allelenum])

		rawdata = []
		f = open(in_path+'/alleleMatrix_new')
		rawdata = f.readlines()
		rawdata = np.asarray(rawdata)
		f.close()

		selrawdata = rawdata[top_userdef_feature]
		Lselrawdata = len(selrawdata)
		
		newpath = '/top' + str(u) + '_locus'
		nowpath = in_path
		pathlocal = nowpath + newpath
		
		if not os.path.exists(pathlocal):
 		   os.makedirs(pathlocal)
#		os.chdir(pathlocal)

		f = open(pathlocal+'/selected_locus.group' + repr(Group_num),'w')
		for i in top_userdef_feature:
			f.write(AllAllele['A0'][i] + '\n')
		f.close()
#		os.chdir(nowpath)

		Group_num += 1
		status = True

	else:
		status = False
		pass

a = []
Group_num = 1
nowpath = in_path
status = True
#os.chdir(nowpath + '/Top' + str(u) + '_var05scheme_selected_locus')

while status == True:
	pathlocal = nowpath + '/top' + str(u) + '_locus' + '/selected_locus.group' + repr(Group_num)

	if os.path.exists(pathlocal):
		f = open(pathlocal)
		for i in f:
			a.append(i)
		Group_num += 1
		status = True
	else:
		status = False
		pass

a = np.unique(a)
f = open(in_path+'/top'+str(u)+'.scheme','w')
for i in a:
  f.write(str(i))
f.close()
#os.chdir(nowpath)



