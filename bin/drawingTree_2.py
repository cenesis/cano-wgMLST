#!/usr/local/bin/python

import ete3
from ete3 import Tree, TreeStyle, TextFace, NodeStyle, CircleFace
import os 
import sys

if len(sys.argv) < 2:
  print('\nUsage: xvfb-run python '+sys.argv[0]+' in/out_path'+' [target_name]'+' [scheme_name]\n')
  print('  xvfb-run python '+sys.argv[0]+' '+'./profileCMP '+'Salmonella_enterica '+'core_top5\n')
  sys.exit()

in_path = sys.argv[1]
if len(sys.argv) >=3: target_name = sys.argv[2]
if len(sys.argv) >=4: scheme_name = sys.argv[3]

f = open(in_path+"/new_tree")
t0 = f.readlines()
f.close()

my_dict = {}

file = open(in_path+"/contigfiles.list", "r") 
for line in file: 
  line = line.strip()
# print line
  st1, st2 = line.split("\t")
  st1 = st1.replace(".fa", "")
  st2 = st2.replace(".fa", "")
  my_dict[st1] = st2
# print st1, st2, my_dict[st1]
  t0[0] = t0[0].replace(st1, st2)
file.close()
my_dict["A00000000"] = "Root"


f = open(in_path+"/new_tree")
t = f.readlines()
f.close()

#t[0] = t[0].strip()
#t[0] = t[0].replace(";", "")
#t[0]=str("(")+t[0]
#t[0]=t[0]+str(");")
#t[0]=t[0]+str(":0,A00000000:0):0;")
#print t[0]
#exit()

t = Tree(t[0])

ts = TreeStyle()
ts.show_leaf_name = False
ts.show_branch_length = False
ts.show_branch_support = True

ts.scale = 2
#ts.min_leaf_separation = 3
ts.branch_vertical_margin = 12

ts.legend_position = 4
#ts.legend.add_face(CircleFace(3, "red"), column=0)
mark = TextFace("Outbreak", fsize=10, fgcolor="red")
mark.margin_top = 10
mark.margin_right = 10
mark.margin_left = 5
mark.margin_bottom = 10
#ts.legend.add_face(mark, column=1)

mark2 = TextFace("X", fsize=10, fgcolor="white")
# Set some attributes
mark2.margin_top = 0
mark2.margin_right = 1
mark2.margin_left = 1
mark2.margin_bottom = 0
mark2.opacity = 1 # from 0 to 1
mark2.border.width = 0
mark2.background.color = "white"
#ts.legend.add_face(mark2, column=0)

mark3 = TextFace("Selected branches", fsize=10, fgcolor="white")
mark3.margin_top = 2
mark3.margin_right = 20
mark3.margin_left = 5
mark3.margin_bottom = 2
ts.legend.add_face(mark3, column=1)

ts.margin_left = 20
ts.margin_right = 20
ts.margin_top = 10
ts.margin_bottom = 10

if len(sys.argv) >=3:
  title = TextFace(target_name, fsize=16, fgcolor="SteelBlue", fstyle="italic", bold=True)
  title.margin_top = 10
  title.margin_right = 10
  title.margin_left = 10
  title.margin_bottom = 10
  ts.title.add_face(title, column=0)

if len(sys.argv) >=4:
  title = TextFace(scheme_name, fsize=14, fgcolor="SteelBlue", fstyle="normal", bold=False)
  title.margin_top = 10
  title.margin_right = 10
  title.margin_left = 10
  title.margin_bottom = 10
  ts.title.add_face(title, column=1)

i=0
for n in t.iter_descendants(): #iterate through all the nodes of the tree
  if n.is_leaf():
      cl0 = TextFace(my_dict[n.name], fsize=12, fgcolor="black")
      n.add_face(cl0, 0, "branch-right")

  #skip if node contain not 2 descendants, so node is leaf (final node)
  if len( n.get_children() )<2: continue 

  #store 2 subtrees as d1 and d2
  d1,d2=n.get_children()
  #check if both subtrees contain at least 3 species (3 leaves)
  if len( d1.get_leaf_names() )<3 and len( d2.get_leaf_names() )<3: continue #if not, skip

  i+=1

  cl = TextFace(i, fsize=8, fgcolor="black")
  # Set some attributes
  cl.margin_top = 2
  cl.margin_right = 2
  cl.margin_left = 2
  cl.margin_bottom = 2
  cl.opacity = 1 # from 0 to 1
#  cl.inner_border.width = 1 # 1 pixel border
#  cl.inner_border.type = 1  # dashed line
  cl.border.width = 1
  cl.background.color = "LightGreen"     

t.render(in_path+"/tree.pdf", w=150, units="mm", tree_style=ts)
t0 = Tree(t0[0])
t0.write(format=2, outfile=in_path+"/tree.newick")

