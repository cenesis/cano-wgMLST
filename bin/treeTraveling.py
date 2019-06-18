#!/usr/local/bin/python

import ete3
from ete3 import Tree
import numpy as np
import pandas as pd
import json
from pandas import read_csv
import os
import sys

if len(sys.argv) != 5:
  print('\nUsage: '+'python '+sys.argv[0]+' '+'newick_tree_file ' + 'userScheme_file '+'alleleMatrix_file '+'out_path\n')
  print('  python '+sys.argv[0]+' '+'./profileCMP/new_tree '+'./profileCMP/core.scheme '+'./profileCMP/alleleMatrix '+'./topSchemeSel\n')
  sys.exit()

newick_tree_file = sys.argv[1]
userScheme_file = sys.argv[2]
alleleMatrix_file = sys.argv[3]
out_path = sys.argv[4]

f = open(newick_tree_file)
t = f.readlines()
f.close()

t[0] = t[0].strip()
t[0] = t[0].replace(";", "")
t[0]=str("(")+t[0]
#t[0]=t[0]+str(");")
t[0]=t[0]+str(":0,A00000000:0):0;")

t = Tree(t[0])
subtree = []

#get all leaf names - species names
species=t.get_leaf_names()

#iterate all nodes
subtrees=[] #will store pairs of subtrees objects with at least 3 species
i=0
sort_group_array = []

#iterate through all the nodes of the tree
for n in t.iter_descendants():
  i += 1
  #skip if node contain not 2 descendants, so node is leaf (final node)
  if len( n.get_children() ) < 2:
    continue 
  #store 2 subtrees as d1 and d2
  d1,d2 = n.get_children()

# check if both subtrees contain at least 3 species (3 leaves)
  if len( d1.get_leaf_names() ) < 3 and len( d2.get_leaf_names() )< 3:
    continue #if not, skip

# Add both subtrees to subtrees list, print iteration, and show node at which tree was split
  subtrees.append( (d1,d2) )

# Sort descendant
  sortd1 = sorted(d1.get_leaf_names())
  sortd2 = sorted(d2.get_leaf_names())
  sort_group_array.append(sortd1)
  sort_group_array.append(sortd2)

# Get the allele profile A000000001 -> int(001) -> 1 as index
length_of_allelename = len(sort_group_array[0][0])  
for idx,val in enumerate(sort_group_array):
  for idx2,val2 in enumerate(sort_group_array[idx]):
    sort_group_array[idx][idx2] = int(sort_group_array[idx][idx2][length_of_allelename-3:length_of_allelename])
sort_group_array = np.asarray(sort_group_array)

# Output subtree leaf list index
f = open(out_path+'/sort_group_array','w')
for i in sort_group_array:
  f.write(str(i)+'\n')
f.close()

##### Grouping part #####
Group1 = []
Group2 = []

allele_dataframe = read_csv(out_path+'/sort_group_array',sep = '\n',header = None)

# Seperate group1 or group2
for idx,val in enumerate(allele_dataframe[0]):
  if idx %2 == 0:
    Group1.append(allele_dataframe[0][idx])
  else:
    Group2.append(allele_dataframe[0][idx])

# Change the data type to array
for idx,val in enumerate(Group1):
  Group1[idx] = json.loads(Group1[idx])
  Group1[idx] = np.asarray(Group1[idx])

for idx,val in enumerate(Group2):
  Group2[idx] = json.loads(Group2[idx])
  Group2[idx] = np.asarray(Group2[idx])
# Get amount of the group element
LGroup1 = range(0,len(Group1))
LGroup2 = range(0,len(Group2))

### Variance threshold scheme ###
f = open(userScheme_file)
scheme = f.readlines()
f.close()
scheme = np.asarray(scheme)
scheme = scheme[1::]

index = [] 
for i in scheme :
  i = int(i[3::]) - 1 
  index.append(i)
index = np.asarray(index)

allele_dataframe = read_csv(alleleMatrix_file,sep = '\t',header = None)
allele_matrix = allele_dataframe.values
allele_matrix = allele_matrix[index]
Rallele_index = len(allele_matrix)
Callele_index = len(allele_matrix[0])
class_labels_int = []

### Classify subgroups ###
for idx in LGroup1:
  class_labels = np.zeros(Callele_index-1)
  class_labels[Group1[idx]-1] = 1
  class_labels[Group2[idx]-1] = 2
  class_labels = class_labels.astype(np.int64)
  class_labels = class_labels.astype(np.object,copy = False)
  class_labels = np.insert(class_labels,0,'class')

  divided = []
  for ind,val in enumerate(class_labels[1:]):
    if val != 0:
      divided.append(ind+1)
    else:
      continue

  divided = np.asarray(divided)
  divided = np.insert(divided,0,0)
  allele_matrix = np.insert(allele_matrix, Rallele_index, class_labels ,axis = 0)
  allele_matrix_new = pd.DataFrame(allele_matrix[0:Rallele_index,:])
  allele_matrix_new.to_csv(out_path+'/alleleMatrix_new', sep = '\t', header = None, index = False)
  allele_matrix_grouped = pd.DataFrame(allele_matrix[0:Rallele_index+1,divided])
  allele_matrix_grouped.to_csv(out_path+'/alleleMatrix_grouped' + repr((idx+1)), sep = '\t', header = None, index = False)

newpath = '/alleleMatrix_grouped'
nowpath = out_path
pathlocal = nowpath + newpath
os.makedirs(pathlocal)
for i in range(1,idx+2):
  os.rename(nowpath + '/alleleMatrix_grouped' + repr(i), pathlocal + '/alleleMatrix_grouped' + repr(i))

