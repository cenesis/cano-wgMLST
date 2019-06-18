#!/usr/local/bin/python

import ete3
from ete3 import Tree, faces, TreeStyle, TextFace, NodeStyle, CircleFace
from ete3 import ClusterTree, RectFace, AttrFace, ProfileFace
from ete3.treeview.faces import add_face_to_node
import pandas as pd
import numpy as np
import colorsys
from IPython.display import display, Image
import os 
import sys

if len(sys.argv) < 2:
	print('\nUsage: /usr/bin/xvfb-run python '+sys.argv[0]+' in/out_path\n')
	print('  /usr/bin/xvfb-run python '+sys.argv[0]+' '+'./topSchemeSel\n')
	sys.exit()

in_path = sys.argv[1]

nameFace = AttrFace("name", fsize=12) #Set leaf node attribute

def setup_heatmap(tree, tree_style, header, center_value=0, nameMap ={}, nameLabel = '', color_up=0.7, color_down=0.2, color_center="white"):
	DEFAULT_COLOR_SATURATION = 0.5
	BASE_LIGHTNESS = 0.7
	def gradient_color(value, max_value, saturation=0.5, hue=0.1):    
		def rgb2hex(rgb):
			return '#%02x%02x%02x' % rgb
		def hls2hex(h, l, s):
			return rgb2hex( tuple(map(lambda x: int(x*255), colorsys.hls_to_rgb(h, l, s))))
	
		lightness = 1 - (value * BASE_LIGHTNESS) / max_value
		return hls2hex(hue, lightness, DEFAULT_COLOR_SATURATION)


	# Calculate max gradient value from the ClusterTree matrix
	maxv = abs(center_value - tree.arraytable._matrix_max)
	minv = abs(center_value - tree.arraytable._matrix_min)
	if center_value <= tree.arraytable._matrix_min:
		MAX_VALUE = minv + maxv
	else:
		MAX_VALUE = max(maxv, minv)
		
	# Add heatmap colors to tree
	cols_add_before_heat = 0
	if nameMap:
		cols_add_before_heat = 1
	for lf in tree:
		if nameMap:
			longNameFace = faces.TextFace(nameMap.get(lf.name, lf.name))
			lf.add_face(longNameFace, column=0, position="aligned")
			
		for i, value in enumerate(getattr(lf, "profile", [])):
			if value > center_value:
				color = gradient_color(abs(center_value - value), MAX_VALUE, hue=color_up)
			elif value < center_value:
				color = gradient_color(abs(center_value - value), MAX_VALUE, hue=color_down)
			else:
				color = center_value
				#color = "LightGreen"
			lf.add_face(RectFace(25, 25, "white", color), position="aligned", column=i+cols_add_before_heat)
			# Uncomment to add numeric values to the matrix
			#lf.add_face(TextFace("%0.2f "%value, fsize=5), position="aligned", column=i)
		lf.add_face(nameFace, column=i+cols_add_before_heat+1, position="aligned")
		
	if nameMap and nameLabel:
		nameF = TextFace(nameLabel, fsize=10)
		#nameF.rotation = -90
		tree_style.aligned_header.add_face(nameF, column=0)
	# Add header 
	for i, name in enumerate(header):
		nameF = TextFace(name, fsize=10)
		nameF.rotation = -90
		tree_style.aligned_header.add_face(nameF, column=i+cols_add_before_heat)



data = pd.read_table(in_path+"/profiles.csv", header=0, index_col=0)
data.index.name = "#Names"
data_mat = data.to_csv(None, sep="\t", float_format="%d")
header = list(data.columns.values)

f = open(in_path+"/wgMLST_tree.newick")
nkTree = f.readlines()
f.close()

t_str = nkTree[0]

t = ClusterTree(t_str, data_mat)

ts = TreeStyle()

ts.margin_left = 20
ts.margin_right = 20
ts.margin_top = 20
ts.margin_bottom = 10

ts.scale = 2
ts.min_leaf_separation = 0
ts.branch_vertical_margin = 0

ts.show_leaf_name = True
ts.show_branch_length = False
ts.show_branch_support = True

setup_heatmap(t, ts, header, center_value=4.5, color_up=0.9, color_down=0.56, color_center="white")

t.render(file_name=in_path+"/heatmap.pdf", tree_style=ts)

