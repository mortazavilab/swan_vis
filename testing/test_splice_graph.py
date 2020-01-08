import pytest
import sys
import numpy as np
sys.path.append('../utils/')
sys.path.append('../../refactor_splice_graph/')
import SpliceGraph as sw
import PlottedGraph as pg
from utils import *
from plotting_tools import * 
import networkx as nx
import math

class TestSpliceGraph(object):
	
	###########################################################################
	######################### Merging tests ###################################
	###########################################################################


	# TODO add_dataset test
	# def test_add_dataset(self):
	# 	a = get_dummy_sg()

	# TODO add_annotation test

	# TODO add_dataset (with abundance) test

	# test for merge_dfs, the following is a summative assesment of 
	# many methods, individual tests for which are also present in this file
	def test_merge_dfs(self):

		a,b = get_dummy_merge_sgs()

		a.merge_dfs(b, 'b')

		print(a.loc_df)
		print(a.edge_df)
		print(a.t_df)

		# loc_df
		print('loc_df')
		test = a.loc_df.apply(
			lambda x: (x.vertex_id, x.a, x.b), axis=1)
		control = [(0,True,True),(1,True,True),(2,True,False),
				   (3,True,False),(4,False,True),(5,False,True),
				   (6,False,True)]
		check_pairs(control, test)

		# edge_df
		print('edge_df')
		test = a.edge_df.apply(
			lambda x: (x.edge_id, x.a, x.b), axis=1)
		control = [((0,1),True,True),((1,2),True, False),((0,2),True,False),
			       ((2,3),True,False),((0,4),False,True),((1,5),False,True),
				   ((4,5),False,True)]
		check_pairs(control, test)

		# t_df
		print('t_df')
		test = a.t_df.apply(
				lambda x: (x.tid, x.a, x.b), axis=1)
		control = [(0,True,True),(1,True,False),(2,False,True),
				   (3,True,False),(4,False,True)]
		check_pairs(control, test)

	# test merge loc_dfs
	def test_merge_loc_dfs(self):
		a,b = get_dummy_merge_sgs()

		a.merge_loc_dfs(b, 'b')

		test = a.loc_df.apply(
			lambda x: (x.vertex_id_a, x.vertex_id_b), axis=1)
		test = [nan_to_neg(pair) for pair in test]
		control = [(0,1),(1,2),(2,np.nan),(3,np.nan),(np.nan,3),
				   (np.nan,4),(np.nan,5)]
		control = [nan_to_neg(pair) for pair in control]
		check_pairs(control, test)

	def test_get_merged_id_map(self):
		a,b = get_dummy_merge_sgs()
		a.merge_loc_dfs(b, 'b')
		id_map = a.get_merged_id_map()

		test = list(id_map.items())
		control = [(1,0),(2,1),(3,4),(4,5),(5,6)]
		check_pairs(control, test)

	def test_update_loc_df_ids_merge(self):
		a,b = get_dummy_merge_sgs()
		a.merge_loc_dfs(b, 'b')
		id_map = a.get_merged_id_map()
		a.update_loc_df_ids_merge(id_map)

		test = a.loc_df.apply(
			lambda x: (x.vertex_id, x.a, x.b), axis=1)
		control = [(0,True,True),(1,True,True),(2,True,False),
				   (3,True,False),(4,False,True),(5,False,True),
				   (6,False,True)]
		check_pairs(control, test)

	def test_update_edge_df_ids_merge(self):
		a,b = get_dummy_merge_sgs()
		a.merge_loc_dfs(b, 'b')
		id_map = a.get_merged_id_map()
		a.update_loc_df_ids_merge(id_map)
		b.update_edge_df_ids_merge(id_map)

		test = b.edge_df.edge_id.tolist()
		control = [(0,1),(0,4),(1,5),(4,5)]
		check_pairs(control, test)

	def test_merge_edge_dfs(self):
		a,b = get_dummy_merge_sgs()
		a.merge_loc_dfs(b, 'b')
		id_map = a.get_merged_id_map()
		a.update_loc_df_ids_merge(id_map)
		b.update_edge_df_ids_merge(id_map)
		a.merge_edge_dfs(b, 'b')

		test = a.edge_df.apply(
			lambda x: (x.edge_id, x.a, x.b), axis=1)
		control = [((0,1),True,True),((1,2),True, False),((0,2),True,False),
			       ((2,3),True,False),((0,4),False,True),((1,5),False,True),
				   ((4,5),False,True)]
		check_pairs(control, test)

	def test_update_t_df_paths_merge(self):
		a,b = get_dummy_merge_sgs()
		a.merge_loc_dfs(b, 'b')
		id_map = a.get_merged_id_map()
		a.update_loc_df_ids_merge(id_map)
		b.update_edge_df_ids_merge(id_map)
		a.merge_edge_dfs(b, 'b')
		b.update_t_df_paths_merge(id_map)

		test = b.t_df.path.tolist()
		test = [tuple(path) for path in test]
		control = [(0,1),(0,1,5),(0,1,4)]
		check_pairs(control, test)

	def test_merge_t_dfs(self):
		a,b = get_dummy_merge_sgs()
		a.merge_loc_dfs(b, 'b')
		id_map = a.get_merged_id_map()
		a.update_loc_df_ids_merge(id_map)
		b.update_edge_df_ids_merge(id_map)
		a.merge_edge_dfs(b, 'b')
		b.update_t_df_paths_merge(id_map)
		a.merge_t_dfs(b, 'b')

		test = a.t_df.apply(
				lambda x: (x.tid, x.a, x.b), axis=1)
		control = [(0,True,True),(1,True,False),(2,False,True),
				   (3,True,False),(4,False,True)]
		check_pairs(control, test)

	############################################################################
	####################### Creation utility tests #############################
	############################################################################
	
	# TODO test create_dfs_gtf (will need dummy input file)

	# TODO test create_dfs_db (will need dummy input file)

	def test_get_loc_types(self):
		a = get_dummy_sg('no_locs')
		a.get_loc_types()

		control_tss = [0,1]
		control_alt_tss = [0,1]
		control_tes = [1,2]
		control_alt_tes = [1,2]
		control_internal = [1]

		test_tss = a.loc_df[a.loc_df.apply(
			lambda x: True if x.TSS else False, axis=1)].vertex_id.tolist()
		test_alt_tss = a.loc_df[a.loc_df.apply(
			lambda x: True if x.alt_TSS else False, axis=1)].vertex_id.tolist()
		test_tes = a.loc_df[a.loc_df.apply(
			lambda x: True if x.TES else False, axis=1)].vertex_id.tolist()
		test_alt_tes = a.loc_df[a.loc_df.apply(
			lambda x: True if x.alt_TES else False, axis=1)].vertex_id.tolist()
		test_internal = a.loc_df[a.loc_df.apply(
			lambda x: True if x.internal else False, axis=1)].vertex_id.tolist()

		print('tss')
		check_pairs(control_tss, test_tss)
		print('alt_tss')
		check_pairs(control_alt_tss, test_alt_tss)
		print('tes')
		check_pairs(control_tes, test_tes)
		print('alt_tes')
		check_pairs(control_alt_tes, test_alt_tes)
		print('internal')
		check_pairs(control_internal, test_internal)

	# This stuff lives in test_graph.py now because these methods live there now
	# def test_update_ids(self):
	# 	a = get_dummy_sg()
	# 	a.update_ids()

	# 	# loc_df
	# 	print('loc_df')
	# 	test = a.loc_df.apply(
	# 		lambda x: (x.vertex_id, x.coord), axis=1)
	# 	control = [(0,1),(1,2),(2,3)]
	# 	check_pairs(control, test)

	# 	# edge_df
	# 	print('edge_df')
	# 	test = a.edge_df.edge_id.tolist()
	# 	control = [(0,1),(0,2),(2,1)]
	# 	check_pairs(control, test)

	# 	# t_df
	# 	print('t_df')
	# 	test = a.t_df.path.tolist()
	# 	test = [tuple(path) for path in test]
	# 	control = [(0,2,1),(2,1),(0,2)]
	# 	check_pairs(control, test)

	# def test_update_loc_df_ids(self):
	# 	a = get_dummy_sg()
	# 	id_map = a.get_ordered_id_map()
	# 	a.update_loc_df_ids(id_map)

	# 	test = a.loc_df.apply(
	# 		lambda x: (x.vertex_id, x.coord), axis=1)
	# 	control = [(0,1),(1,2),(2,3)]
	# 	check_pairs(control, test)

	# def test_update_edge_df_ids(self):
	# 	a = get_dummy_sg()
	# 	id_map = a.get_ordered_id_map()
	# 	a.update_edge_df_ids(id_map)

	# 	test = a.edge_df.edge_id.tolist()
	# 	control = [(0,1),(0,2),(2,1)]
	# 	check_pairs(control, test)

	# def test_update_t_df_paths(self):
	# 	a = get_dummy_sg()
	# 	id_map = a.get_ordered_id_map()
	# 	a.update_t_df_paths(id_map)

	# 	test = a.t_df.path.tolist()
	# 	test = [tuple(path) for path in test]
	# 	control = [(0,2,1),(2,1),(0,2)]
	# 	check_pairs(control, test)		

	# # TODO extend this to work on multiple genes,
	# # will require a more robust test
	# def test_get_ordered_id_map(self):
	# 	a = get_dummy_sg()
	# 	id_map = a.get_ordered_id_map()

	# 	test = list(id_map.items())
	# 	control = [(0,0),(1,2),(2,1)]
	# 	check_pairs(control, test)

	############################################################################
	###################### Ordering transcripts tests ##########################
	############################################################################

	def test_order_transcripts(self):
		a = get_dummy_sg()
		a.t_df['a_counts'] =[0,0,12]
		a.t_df['b_counts'] = [1,0,14]
		a.counts = ['a_counts', 'b_counts']

		# order by expression level
		a.order_transcripts(order='expression')
		print(a.t_df.head())
		assert a.t_df.tid.tolist() == [0,2,1]

		# order by transcript id
		a.order_transcripts()
		print(a.t_df.head())
		assert a.t_df.tid.tolist() == [0,1,2]

		# order by coordinate of tss
		a.order_transcripts(order='tss')
		print(a.t_df.head())
		assert a.t_df.tid.tolist() == [0,2,1]

		# order by coordinate of tes
		a.order_transcripts()
		a.order_transcripts(order='tes')

		print()
		print(a.loc_df)
		print(a.edge_df)
		print(a.t_df)
		assert a.t_df.tid.tolist() == [1,2,0]

	############################################################################
	####################### Adding counts/tpms test ############################
	############################################################################

	def test_add_abundance(self):
		a = get_dummy_sg(special='dataset')
		a.add_abundance('input_files/test_abundance.tsv',
			['count_1a', 'count_2a'], 'dataset_a')

		assert 'dataset_a_counts' in a.t_df.columns
		assert 'dataset_a_tpm' in a.t_df.columns
		assert a.get_count_cols() == ['dataset_a_counts']
		assert a.get_tpm_cols() == ['dataset_a_tpm']

def nan_to_neg(pair):
	new_tuple = []
	for i in pair:
		if math.isnan(i):
			new_tuple.append(-1)
		else: 
			new_tuple.append(i)
	return tuple(new_tuple)


def check_pairs(control, test):
	print('control')
	print(control)
	for t in test:
		print(t)
		assert t in control
	assert len(test) == len(control)

def get_dummy_sg(special=None):
	a = sw.SpliceGraph()

	loc_df = pd.DataFrame({'chrom': [1,1,1],
		'coord': [1,3,2],
		'strand': ['+', '+', '+'],
		'vertex_id': [0,1,2]})
	loc_df = create_dupe_index(loc_df, 'vertex_id')
	loc_df = set_dupe_index(loc_df, 'vertex_id')

	edge_df = pd.DataFrame({'edge_id': [(0,2),(0,1),(1,2)],
			  'v1': [0,0,1],
			  'v2': [2,1,2],
			  'edge_type': ['exon', 'exon', 'intron'], 
			  'strand': ['+','+','+']})
	edge_df = create_dupe_index(edge_df, 'edge_id')
	edge_df = set_dupe_index(edge_df, 'edge_id')

	t_df = pd.DataFrame({'tid':[2,1,0],
			  'gid':[0,0,0],
			  'gname':['0','0','0'],
			  'path':[[0,1,2],[1,2],[0,1]],
			  'counts_a': [0,0,12],
			  'counts_b': [1,0,14]})
	t_df = create_dupe_index(t_df, 'tid')
	t_df = set_dupe_index(t_df, 'tid')

	if special == 'intron':
		edge_df.loc[(0,1), 'edge_type'] = 'intron'


	a.loc_df = loc_df
	a.edge_df = edge_df
	a.t_df = t_df 

	if special == 'no_locs':	
		pass
	else:
		a.get_loc_types()

	if special == 'dataset':
		a.datasets = ['dataset_a']
		a.loc_df['dataset_a'] = True
		a.edge_df['dataset_a'] = True
		a.t_df['dataset_a'] = True

	return a

def get_dummy_merge_sgs():
		a = sw.SpliceGraph()
		b = sw.SpliceGraph()

		a.loc_df = pd.DataFrame({'chrom': [1,1,1,2],
								 'coord': [1,2,3,1],
								 'strand': ['+', '+', '+', '+'],
								 'vertex_id': [0,1,2,3]})
		b.loc_df = pd.DataFrame({'chrom': [1,1,1,2,3],
								 'coord': [1,2,4,1,1],
								 'strand': ['+', '+', '+', '-', '+'],
								 'vertex_id': [1,2,3,4,5]})
		a.loc_df = create_dupe_index(a.loc_df, 'vertex_id')
		b.loc_df = create_dupe_index(b.loc_df, 'vertex_id')
		a.loc_df = set_dupe_index(a.loc_df, 'vertex_id')
		b.loc_df = set_dupe_index(b.loc_df, 'vertex_id')

		a.edge_df = pd.DataFrame({'edge_id': [(0,1),(1,2),(0,2),(2,3)],
								  'v1': [0,1,0,2],
								  'v2': [1,2,2,3],
								  'edge_type': ['exon', 'intron', 'exon', 'exon'], 
								  'strand': ['+','+','+','+']})
		b.edge_df = pd.DataFrame({'edge_id': [(1,2),(1,3),(2,4),(3,4)],
								  'v1': [1,1,2,3],
								  'v2': [2,3,4,4],
								  'strand': ['+','+','+','+'],
								  'edge_type': ['exon', 'exon', 'intron', 'intron']})
		a.edge_df = create_dupe_index(a.edge_df, 'edge_id')
		b.edge_df = create_dupe_index(b.edge_df, 'edge_id')
		a.edge_df = set_dupe_index(a.edge_df, 'edge_id')
		b.edge_df = set_dupe_index(b.edge_df, 'edge_id')

		a.t_df = pd.DataFrame({'tid': [0,1,3],
								'gid': [0,0,0],
								'gname': ['0','0','0'],
								'path': [[0,1],[0,1,2,3],[0,2]]})
		b.t_df = pd.DataFrame({'tid': [0,2,4],
								'gid': [0,0,0],
								'gname': ['0','0','0'],
								'path': [[1,2],[1,2,4],[1,2,3]]})
		a.t_df = create_dupe_index(a.t_df, 'tid')
		b.t_df = create_dupe_index(b.t_df, 'tid')
		a.t_df = set_dupe_index(a.t_df, 'tid')
		b.t_df = set_dupe_index(b.t_df, 'tid')

		# add 'dataset a' to a
		a.datasets = ['a']
		a.loc_df['a'] = True
		a.edge_df['a'] = True
		a.t_df['a'] = True

		a.get_loc_types()
		b.get_loc_types()


		return a, b




