import pytest
import sys
import numpy as np
import swan_vis as swan
import pandas as pd
from swan_vis.utils import *

class TestPlotting(object):
	def test_plotting(self):
		# set up testing swangraph
		sg = swan.SwanGraph()
		sg.datasets = ['annotation', 'a', 'b']
		gene1_loc_df = pd.DataFrame({'chrom': [1, 1, 1, 1, 1, 1, 1],
			'coord': [5, 10, 15, 20, 25, 30, 35],
			'strand': ['+', '+', '+', '+', '+', '+', '+'],
			'vertex_id': [0, 1, 2, 3, 4, 5, 6],
			'annotation': [True, True, True, True, True, True, False],
			'a': [True, True, False, True, True, True, True],
			'b': [False, False, True, True, True, True, False]})
		gene1_t_df = pd.DataFrame({
			'gname': ['GENE01', 'GENE01', 'GENE01', 'GENE01'],
			'gid': ['ENSG01', 'ENSG01', 'ENSG01', 'ENSG01'],
			'tid': ['ENST01', 'ENST02', 'ENST03', 'ENST04'],
			'path': [[0,1,2,3,4,5], [0,3,4,5], [0,1,4,6], [2,3,4,5]],
			'annotation': [True, True, True, False],
			'a': [False, True, True, False], 
			'b': [False, False, False, True]}) 
		gene1_edge_df = pd.DataFrame({
			'v1': [0,1,2,3,4,0,4,1],
			'v2': [1,2,3,4,5,3,6,4],
			'edge_type': ['exon','intron','exon','intron','exon','exon','exon','intron'], 
			'annotation': [True, True, True, True, True, True, False, False],
			'a': [True, False, False, True, True, True, True, True],
			'b': [False, False, True, True, True, False, False, False]
			})
		gene1_edge_df['edge_id'] = gene1_edge_df.apply(lambda x: 
			(x.v1, x.v2), axis=1)

		sg.loc_df = gene1_loc_df
		sg.edge_df = gene1_edge_df
		sg.t_df = gene1_t_df

		sg.loc_df = create_dupe_index(sg.loc_df, 'vertex_id')
		sg.loc_df = set_dupe_index(sg.loc_df, 'vertex_id')
		sg.edge_df = create_dupe_index(sg.edge_df, 'edge_id')
		sg.edge_df = set_dupe_index(sg.edge_df, 'edge_id')
		sg.t_df = create_dupe_index(sg.t_df, 'tid')
		sg.t_df = set_dupe_index(sg.t_df, 'tid')
		sg.get_loc_types()

		# testing
		gene1_tids = gene1_t_df.tid.tolist()
		gene1_locs = gene1_loc_df.vertex_id.tolist()
		gene1_edges = gene1_edge_df.edge_id.tolist()

		# 0th plot - gene summary graph of ENSG01
		sg.datasets = ['annotation', 'a']
		sg = plot0(sg, gene1_tids, gene1_locs, gene1_edges)

		gene2_loc_df = pd.DataFrame({'chrom': [4, 4, 4, 4, 4, 4, 4],
			'coord': [35, 30, 25, 20, 15, 10, 5],
			'strand': ['-', '-', '-', '-', '-', '-', '-'],
			# 'vertex_id': [0, 1, 2, 3,  4,  5,  6],
			'vertex_id': [7, 8, 9, 10, 11, 12, 13],
			'annotation': [True, True, True, True, True, True, False],
			'a': [True, True, False, True, True, True, True],
			'b': [False, False, True, True, True, True, False]})
		gene2_t_df = pd.DataFrame({
			'gname': ['GENE02', 'GENE02', 'GENE02', 'GENE02'],
			'gid': ['ENSG02', 'ENSG02', 'ENSG02', 'ENSG02'],
			'tid': ['ENST05', 'ENST06', 'ENST07', 'ENST08'],
			'path': [[7,8,9,10,11,12], [7,10,11,12], [7,8,11,13], [9,10,11,12]],
			'annotation': [True, True, True, False],
			'a': [False, True, True, False], 
			'b': [False, False, False, True]}) 
		gene2_edge_df = pd.DataFrame({
			'v1': [7,8,9,10,11,7,11,8],
			'v2': [8,9,10,11,12,10,13,11],
			'edge_type': ['exon','intron','exon','intron','exon','exon','exon','intron'], 
			'annotation': [True, True, True, True, True, True, False, False],
			'a': [True, False, False, True, True, True, True, True],
			'b': [False, False, True, True, True, False, False, False]
			})
		gene2_edge_df['edge_id'] = gene2_edge_df.apply(lambda x: 
			(x.v1, x.v2), axis=1)

		sg.loc_df = pd.concat([gene1_loc_df, gene2_loc_df])
		sg.edge_df = pd.concat([gene1_edge_df, gene2_edge_df])
		sg.t_df = pd.concat([gene1_t_df, gene2_t_df])

		sg.loc_df = create_dupe_index(sg.loc_df, 'vertex_id')
		sg.loc_df = set_dupe_index(sg.loc_df, 'vertex_id')
		sg.edge_df = create_dupe_index(sg.edge_df, 'edge_id')
		sg.edge_df = set_dupe_index(sg.edge_df, 'edge_id')
		sg.t_df = create_dupe_index(sg.t_df, 'tid')
		sg.t_df = set_dupe_index(sg.t_df, 'tid')
		sg.get_loc_types()

		# testing
		gene2_tids = gene2_t_df.tid.tolist()
		gene2_locs = gene2_loc_df.vertex_id.tolist()
		gene2_edges = gene2_edge_df.edge_id.tolist()

		# remake the same plot and force it to update
		sg.datasets = ['annotation', 'a', 'b']
		sg = plot0_5(sg, gene2_tids, gene1_locs, gene1_edges)

		# first plot - gene summary graph of ENSG01
		sg = plot1(sg, gene1_tids, gene1_locs, gene1_edges)

		# plot a transcript through the same gene
		sg = plot2(sg, gene1_tids, gene1_locs, gene1_edges)

		# make sure we are doing the right thing after plotting ENST01
		# after plotting it as a browser image
		sg = plot3(sg, gene1_tids, gene1_locs, gene1_edges)

		# plot the same transcript with indicate_novel
		sg = plot4(sg, gene1_tids, gene1_locs, gene1_edges)

		# plot a different transcript but change the indicate opt
		sg = plot5(sg, gene1_tids, gene1_locs, gene1_edges)

		# plot a new gene and use indicate_dataset
		sg = plot6(sg, gene2_tids, gene1_locs, gene1_edges)

		# plot a transcript from the other gene using browser
		sg = plot7(sg, gene1_tids, gene1_locs, gene1_edges)

		# plot a transcript from the other gene using
		# indicate_dataset b
		sg = plot8(sg, gene2_tids, gene1_locs, gene1_edges)

		# plot a gene using indicate_dataset

		# add a new dataset to the graph and plot the same gene has not been 
		# in plotted_graph.py yet

		# add a new dataset to the graph and plot a transcript from a different gene

def plot0(sg, tids, locs, edges):
	sg.plot_graph('ENSG01')
	check_subset(sg, tids, edges, locs)
	node_colors_ctrl = [(0, 'light_blue'), (1, 'yellow'),
			   (2, 'light_blue'), (3, 'yellow'), 
			   (4, 'yellow'), (5, 'orange'), (6, 'orange')]
	curve_ctrl = [((0,1),None), ((1,2),None),
			   ((2,3),None), ((3,4),None),
			   ((4,5),None), ((0,3),'curved'),
			   ((4,6),'curved'),((1,4),'curved')]
	edge_color_ctrl = [((0,1),'green'), ((1,2),'pink'),
			   ((2,3),'green'), ((3,4),'pink'),
			   ((4,5),'green'), ((0,3),'green'),
			   ((4,6),'green'),((1,4),'pink')]
	check_nodes(sg.pg.loc_df, node_colors_ctrl)
	check_edges(sg.pg.edge_df, edge_color_ctrl, curve_ctrl)
	return sg

def plot0_5(sg, tids, locs, edges):
	sg.plot_graph('ENSG02', indicate_dataset='a')
	swan.save_fig('scratch/ensg02_dataset')
	check_subset(sg, tids, edges, locs)
	edgecolor_ctrl = [(0, 'k'), (1, 'k'),
			   (2, None), (3, 'k'), 
			   (4, 'k'), (5, 'k'), (6, 'k')]
	node_colors_ctrl = [(0, 'light_blue'), (1, 'yellow'),
			   (2, 'light_blue'), (3, 'yellow'), 
			   (4, 'yellow'), (5, 'orange'), (6, 'orange')]
	curve_ctrl = [((0,1),None), ((1,2),None),
			   ((2,3),None), ((3,4),None),
			   ((4,5),None), ((0,3),'curved'),
			   ((4,6),'curved'),((1,4),'curved')]
	edge_color_ctrl = [((0,1),'green'), ((1,2),'pink'),
			   ((2,3),'green'), ((3,4),'pink'),
			   ((4,5),'green'), ((0,3),'green'),
			   ((4,6),'green'),((1,4),'pink')]
	style_ctrl = [((0,1),'dashed'), ((1,2),None),
			   ((2,3),None), ((3,4),'dashed'),
			   ((4,5),'dashed'), ((0,3),'dashed'),
			   ((4,6),'dashed'),((1,4),'dashed')]
	check_nodes(sg.pg.loc_df, node_colors_ctrl, edgecolor_ctrl)
	check_edges(sg.pg.edge_df, edge_color_ctrl, curve_ctrl, style_ctrl)
	return sg


def plot1(sg, tids, locs, edges):
	sg.plot_graph('ENSG01')
	check_subset(sg, tids, edges, locs)
	node_colors_ctrl = [(0, 'light_blue'), (1, 'yellow'),
			   (2, 'light_blue'), (3, 'yellow'), 
			   (4, 'yellow'), (5, 'orange'), (6, 'orange')]
	curve_ctrl = [((0,1),None), ((1,2),None),
			   ((2,3),None), ((3,4),None),
			   ((4,5),None), ((0,3),'curved'),
			   ((4,6),'curved'),((1,4),'curved')]
	edge_color_ctrl = [((0,1),'green'), ((1,2),'pink'),
			   ((2,3),'green'), ((3,4),'pink'),
			   ((4,5),'green'), ((0,3),'green'),
			   ((4,6),'green'),((1,4),'pink')]
	check_nodes(sg.pg.loc_df, node_colors_ctrl)
	check_edges(sg.pg.edge_df, edge_color_ctrl, curve_ctrl)
	return sg

def plot2(sg, tids, locs, edges):
	sg.plot_transcript_path('ENST01')
	swan.save_fig('scratch/enst01.png')
	check_subset(sg, tids, edges, locs)
	node_colors_ctrl = [(0, 'light_blue'), (1, 'yellow'),
			   (2, 'yellow'), (3, 'yellow'), 
			   (4, 'yellow'), (5, 'orange'), (6, 'gray_orange')]
	curve_ctrl = [((0,1),None), ((1,2),None),
			   ((2,3),None), ((3,4),None),
			   ((4,5),None), ((0,3),'curved'),
			   ((4,6),'curved'),((1,4),'curved')]
	edge_color_ctrl = [((0,1),'green'), ((1,2),'pink'),
			   ((2,3),'green'), ((3,4),'pink'),
			   ((4,5),'green'), ((0,3),'gray_green'),
			   ((4,6),'gray_green'),((1,4),'gray_pink')]
	used_edges = [(0,1),(1,2),(2,3),(3,4),(4,5)]
	check_nodes(sg.pg.loc_df, node_colors_ctrl)
	check_edges(sg.pg.edge_df, edge_color_ctrl, curve_ctrl)
	return sg

def plot3(sg, tids, locs, edges):
	sg.plot_transcript_path('ENST01', browser=True)
	sg.plot_transcript_path('ENST01')
	check_subset(sg, tids, edges, locs)
	node_colors_ctrl = [(0, 'light_blue'), (1, 'yellow'),
			   (2, 'yellow'), (3, 'yellow'), 
			   (4, 'yellow'), (5, 'orange'), (6, 'gray_orange')]
	curve_ctrl = [((0,1),None), ((1,2),None),
			   ((2,3),None), ((3,4),None),
			   ((4,5),None), ((0,3),'curved'),
			   ((4,6),'curved'),((1,4),'curved')]
	edge_color_ctrl = [((0,1),'green'), ((1,2),'pink'),
			   ((2,3),'green'), ((3,4),'pink'),
			   ((4,5),'green'), ((0,3),'gray_green'),
			   ((4,6),'gray_green'),((1,4),'gray_pink')]
	used_edges = [(0,1),(1,2),(2,3),(3,4),(4,5)]
	check_nodes(sg.pg.loc_df, node_colors_ctrl)
	check_edges(sg.pg.edge_df, edge_color_ctrl, curve_ctrl)
	return sg

def plot4(sg, tids, locs, edges):
	sg.plot_transcript_path('ENST01', indicate_novel=True)
	swan.save_fig('scratch/ENST01_novel.png')
	check_subset(sg, tids, edges, locs)
	node_colors_ctrl = [(0, 'light_blue'), (1, 'yellow'),
			   (2, 'yellow'), (3, 'yellow'), 
			   (4, 'yellow'), (5, 'orange'), (6, 'gray_orange')]
	edgecolor_ctrl = [(0, None), (1, None),
			   (2, None), (3, None), 
			   (4, None), (5, None), (6, '#999999')]
	curve_ctrl = [((0,1),None), ((1,2),None),
			   ((2,3),None), ((3,4),None),
			   ((4,5),None), ((0,3),'curved'),
			   ((4,6),'curved'),((1,4),'curved')]
	edge_color_ctrl = [((0,1),'green'), ((1,2),'pink'),
			   ((2,3),'green'), ((3,4),'pink'),
			   ((4,5),'green'), ((0,3),'gray_green'),
			   ((4,6),'gray_green'),((1,4),'gray_pink')]
	style_ctrl = [((0,1),None), ((1,2),None),
			   ((2,3),None), ((3,4),None),
			   ((4,5),None), ((0,3),None),
			   ((4,6),'dashed'),((1,4),'dashed')]
	used_edges = [(0,1),(1,2),(2,3),(3,4),(4,5)]
	check_nodes(sg.pg.loc_df, node_colors_ctrl, edgecolor_ctrl)
	check_edges(sg.pg.edge_df, edge_color_ctrl, curve_ctrl, style_ctrl)
	return sg

def plot5(sg, tids, locs, edges):
	sg.plot_transcript_path('ENST04')
	swan.save_fig('scratch/enst04.png')
	check_subset(sg, tids, edges, locs)
	node_colors_ctrl = [(0, 'gray_light_blue'), (1, 'gray_yellow'),
			   (2, 'light_blue'), (3, 'yellow'), 
			   (4, 'yellow'), (5, 'orange'), (6, 'gray_orange')]
	curve_ctrl = [((0,1),None), ((1,2),None),
			   ((2,3),None), ((3,4),None),
			   ((4,5),None), ((0,3),'curved'),
			   ((4,6),'curved'),((1,4),'curved')]
	edge_color_ctrl = [((0,1),'gray_green'), ((1,2),'gray_pink'),
			   ((2,3),'green'), ((3,4),'pink'),
			   ((4,5),'green'), ((0,3),'gray_green'),
			   ((4,6),'gray_green'),((1,4),'gray_pink')]
	used_edges = [(2,3),(3,4),(4,5)]
	check_nodes(sg.pg.loc_df, node_colors_ctrl)
	check_edges(sg.pg.edge_df, edge_color_ctrl, curve_ctrl)
	return sg

def plot6(sg, tids, locs, edges):
	sg.plot_graph('ENSG02', indicate_dataset='a')
	swan.save_fig('scratch/ensg02_dataset')
	check_subset(sg, tids, edges, locs)
	edgecolor_ctrl = [(0, 'k'), (1, 'k'),
			   (2, None), (3, 'k'), 
			   (4, 'k'), (5, 'k'), (6, 'k')]
	node_colors_ctrl = [(0, 'light_blue'), (1, 'yellow'),
			   (2, 'light_blue'), (3, 'yellow'), 
			   (4, 'yellow'), (5, 'orange'), (6, 'orange')]
	curve_ctrl = [((0,1),None), ((1,2),None),
			   ((2,3),None), ((3,4),None),
			   ((4,5),None), ((0,3),'curved'),
			   ((4,6),'curved'),((1,4),'curved')]
	edge_color_ctrl = [((0,1),'green'), ((1,2),'pink'),
			   ((2,3),'green'), ((3,4),'pink'),
			   ((4,5),'green'), ((0,3),'green'),
			   ((4,6),'green'),((1,4),'pink')]
	style_ctrl = [((0,1),'dashed'), ((1,2),None),
			   ((2,3),None), ((3,4),'dashed'),
			   ((4,5),'dashed'), ((0,3),'dashed'),
			   ((4,6),'dashed'),((1,4),'dashed')]
	check_nodes(sg.pg.loc_df, node_colors_ctrl, edgecolor_ctrl)
	check_edges(sg.pg.edge_df, edge_color_ctrl, curve_ctrl, style_ctrl)
	return sg

def plot7(sg, tids, locs, edges):
	sg.plot_transcript_path('ENST02', browser=True)
	check_subset(sg, tids, edges, locs)
	return sg

def plot8(sg, tids, locs, edges):
	sg.plot_transcript_path('ENST08', indicate_dataset='b')
	# TODO

def check_subset(sg, t_ctrl, edge_ctrl, loc_ctrl):

	# todo: need a better way of doing this because update_ids runs
	# when plotting. ie check loc coords or something
	t_test = sg.pg.t_df.tid.tolist()
	loc_test = sg.pg.loc_df.vertex_id.tolist()
	edge_test = sg.pg.edge_df.edge_id.tolist()
	check_pairs(t_ctrl, t_test)
	check_pairs(edge_ctrl, edge_test)
	check_pairs(loc_ctrl, loc_test)

def check_edges(edge_df, color_ctrl, curve_ctrl, style_ctrl=None, used_edges_ctrl=None):
	# color
	color_test = get_colors(edge_df)
	check_pairs(color_ctrl, color_test)

	# curve
	curve_test = get_edge_curve(edge_df)
	check_pairs(curve_ctrl, curve_test)

	# style
	if style_ctrl:
		style_test = get_edge_style(edge_df)
		check_pairs(style_ctrl, style_test)
	else:
		style_test = edge_df.line.tolist()
		assert len(set(style_test)) == 1
		assert style_test[0] == None

	# order
	if used_edges_ctrl:
		n_used_edges = len(used_edges)
		used_edges_test = edge_df.tail(n_used_edges).edge_id.tolist()
		check_pairs(used_edges_ctrl, used_edges_test)

def check_nodes(loc_df, color_ctrl, edgecolor_ctrl=None):
	# color
	color_test = get_colors(loc_df)
	check_pairs(color_ctrl, color_test)

	# edgecolor
	if edgecolor_ctrl:
		edgecolor_test = get_node_edgecolor(loc_df)
		check_pairs(edgecolor_ctrl, edgecolor_test)
	else:
		edgecolor_test = loc_df.edgecolor.tolist()
		assert len(set(edgecolor_test)) == 1
		assert edgecolor_test[0] == None

def get_colors(df):

		# plotting colors
		yellow = '#F0E442'
		light_blue = '#56B4E9'
		orange = '#E69F00'
		pink = '#CC79A7'
		green = '#009E73'

		gray_light_blue = '#c1ddec'
		gray_yellow = '#f8f6db'
		gray_orange = '#ebdec3'
		gray_pink = '#efdae5'
		gray_green = '#cdf5ea'

		color_dict = {yellow: 'yellow', light_blue: 'light_blue', 
			orange: 'orange', pink: 'pink', green: 'green',
			gray_light_blue: 'gray_light_blue', gray_yellow: 'gray_yellow',
			gray_orange: 'gray_orange', gray_pink: 'gray_pink',
			gray_green: 'gray_green'}

		df['hr_color'] = df.apply(lambda x: color_dict[x.color], axis=1)

		if 'vertex_id' in df.columns:
			id_color_pairs = df.apply(lambda x: (x.vertex_id, x.hr_color), axis=1)
		elif 'edge_id' in df.columns:
			id_color_pairs = df.apply(lambda x: (x.edge_id, x.hr_color), axis=1)

		df.drop('hr_color', axis=1, inplace=True)
		return id_color_pairs

def get_edge_curve(edge_df):
	edge_df['hr_curve'] = edge_df.apply(lambda x:
		'curved' if x.curve else None, axis=1)
	id_curve_pairs = edge_df.apply(lambda x: (x.edge_id, x.hr_curve), axis=1)
	return id_curve_pairs

def get_edge_style(edge_df):
	id_style_pairs = edge_df.apply(lambda x:
		(x.edge_id, x.line), axis=1)
	return id_style_pairs

def get_node_edgecolor(loc_df):
	id_edgecolor_pairs = loc_df.apply(lambda x:
		(x.vertex_id, x.edgecolor), axis=1)
	return id_edgecolor_pairs

def check_pairs(control, test):
	print('control')
	print(control)
	print('test')
	print(test)
	for t in test:
		assert t in control
	assert len(test) == len(control)