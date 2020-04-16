import pytest
import sys
import numpy as np
import swan_vis as swan
import networkx as nx
import math
import pandas as pd

# tests the functions in superclass Graph.py
class TestGraph(object):

	##########################################################################
	################# Related to checking contents of Graph ##################
	##########################################################################

	# test for check_datasets
	def test_check_datasets(self):

		g = swan.SwanGraph()
		g.datasets = ['a','b']

		# 1: test for a dataset that's not there
		with pytest.raises(Exception) as excinfo:
			g.check_datasets(['c'])
		assert "c not present in graph" in str(excinfo.value)

		# 2: test for a dataset that is there, input in string form
		result = g.check_datasets('a')
		assert result == None

		# 3: test for multiple datasets that are in the graph
		result = g.check_datasets(['a', 'b'])
		assert result == None

		# 4: test for a dataset that's not there in a list 
		# of datasets that are there
		with pytest.raises(Exception) as excinfo:
			g.check_datasets(['a', 'b', 'c'])
		assert 'c not present in graph' in str(excinfo.value)

	# test for check_abundances
	def test_check_abundances(self):
		g = swan.SwanGraph()
		g.datasets = ['a', 'b', 'c']
		g.counts = ['a_counts', 'b_counts']

		# 1: test for abundance that isn't there
		# also makes sure it works with character input (ie not list)
		with pytest.raises(Exception) as excinfo:
			g.check_abundances('d')
		assert 'Abundance for dataset d' in str(excinfo.value)

		# 2: test for a abundance from a dataset that is there but
		# does not have abundance info
		# also makes sure it works with len(list) = 1 info
		with pytest.raises(Exception) as excinfo:
			g.check_abundances(['c'])			
		assert 'Abundance for dataset c' in str(excinfo.value)

		# 3: test for multiple datasets that are in the graph
		result = g.check_abundances(['a', 'b'])
		assert result == None

		# 4: test for a dataset that's not there in a list 
		# of datasets that are there
		with pytest.raises(Exception) as excinfo:
			g.check_abundances(['a', 'b', 'd'])
		assert 'Abundance for dataset d' in str(excinfo.value)

	##########################################################################
	####################### Related to creating Graph ########################
	##########################################################################

	# test update_ids
	def test_update_ids(self):
		a = get_dummy_sg()
		a.update_ids()

		# loc_df
		print('loc_df')
		test = a.loc_df.apply(
			lambda x: (x.vertex_id, x.coord), axis=1)
		control = [(0,1),(1,2),(2,3)]
		check_pairs(control, test)

		# edge_df
		print('edge_df')
		test = a.edge_df.edge_id.tolist()
		control = [(0,1),(0,2),(2,1)]
		check_pairs(control, test)

		# t_df
		print('t_df')
		test = a.t_df.path.tolist()
		test = [tuple(path) for path in test]
		control = [(0,2,1),(2,1),(0,2)]
		check_pairs(control, test)

	# test update_loc_df_ids
	def test_update_loc_df_ids(self):
		a = get_dummy_sg()
		id_map = a.get_ordered_id_map()
		a.dfs_to_dicts()
		a.update_loc_df_ids(id_map)
		a.dicts_to_dfs()

		test = a.loc_df.apply(
			lambda x: (x.vertex_id, x.coord), axis=1)
		control = [(0,1),(1,2),(2,3)]
		check_pairs(control, test)

	# tests update_edge_df_ids
	def test_update_edge_df_ids(self):
		a = get_dummy_sg()
		id_map = a.get_ordered_id_map()
		a.dfs_to_dicts()
		a.update_edge_df_ids(id_map)
		a.dicts_to_dfs()

		test = a.edge_df.edge_id.tolist()
		control = [(0,1),(0,2),(2,1)]
		check_pairs(control, test)

	# tests update_t_df_paths
	def test_update_t_df_paths(self):
		a = get_dummy_sg()
		id_map = a.get_ordered_id_map()
		a.dfs_to_dicts()
		a.update_t_df_paths(id_map)
		a.dicts_to_dfs()

		test = a.t_df.path.tolist()
		test = [tuple(path) for path in test]
		control = [(0,2,1),(2,1),(0,2)]
		check_pairs(control, test)		

	# TODO extend this to work on multiple genes,
	# will require a more robust test
	# tests get_ordered_id_map
	def test_get_ordered_id_map(self):
		a = get_dummy_sg()
		id_map = a.get_ordered_id_map()

		test = list(id_map.items())
		control = [(0,0),(1,2),(2,1)]
		check_pairs(control, test)

	# tests create_graph_from_dfs
	def test_create_graph_from_dfs(self):
		### TODO ###
		pass
		### TODO ###

	# tests oder_edge_df
	def test_order_edge_df(self):
		a = get_dummy_sg()
		a.order_edge_df()

		v1 = a.edge_df.v1.tolist()
		assert v1 == sorted(v1)

	##########################################################################
	############################# Other utilities ############################
	##########################################################################

	# tests is_empty
	def test_is_empty(self):
		a = swan.SwanGraph()

		# 1: graph is empty, no datasets (including annotation) 
		print(a.datasets)
		assert a.is_empty() == True

		a.datasets = ['a']

		# 2: graph has one dataset
		print(a.datasets)
		assert a.is_empty() == False

		a.datasets = ['a', 'b']

		# 3: graph has more than one dataset
		print(a.datasets)
		assert a.is_empty() == False

	# tests get_dataset_cols
	def test_get_dataset_cols(self):
		a = swan.SwanGraph()

		# 1: graph is empty, no datasets (including annotation) 
		print(a.datasets)
		assert a.get_dataset_cols() == []

		a.datasets = ['a']

		# 2: graph has one dataset
		print(a.datasets)
		assert a.get_dataset_cols() == ['a']

		a.datasets = ['a', 'b']

		# 3: graph has more than one dataset
		print(a.datasets)
		assert a.get_dataset_cols() == ['a', 'b']

	# tests get_count_cols
	def test_get_count_cols(self):
		a = swan.SwanGraph()

		# 1: empty graph, no datasets argument
		assert a.get_count_cols() == []

		# 2: one dataset in graph, no datasets argument
		a.counts = ['a_counts']
		assert a.get_count_cols() == ['a_counts']

		# 3: > one dataset in graph, no datasets argument
		a.counts = ['a_counts', 'b_counts']
		assert a.get_count_cols() == ['a_counts', 'b_counts']

		# 4: empty graph, datasets argument
		# also tests ability to handle char input
		a = swan.SwanGraph()
		with pytest.raises(Exception) as excinfo:
			a.get_count_cols('a')
		assert 'Abundance for dataset a' in str(excinfo.value)

		# 5: graph with one dataset and datasets argument already in graph
		a.counts = ['a_counts']
		assert a.get_count_cols('a') == ['a_counts']

		# 6: graph with one dataset and datasets argument not in graph
		with pytest.raises(Exception) as excinfo:
			a.get_count_cols('b')
		assert 'Abundance for dataset b' in str(excinfo.value)

		# 7: graph with more than one dataset and datasets argument in graph
		a.counts = ['a_counts', 'b_counts']
		assert a.get_count_cols(['a', 'b']) == ['a_counts', 'b_counts']

		# 8: graph with more than one dataset and not all datasets argument in graph
		with pytest.raises(Exception) as excinfo:
			a.get_count_cols(['a', 'b', 'c'])
		assert 'Abundance for dataset c'

	# tests get_tpm_cols
	def test_get_tpm_cols(self):
		a = swan.SwanGraph()

		# 1: empty graph, no datasets argument
		assert a.get_tpm_cols() == []

		# 2: one dataset in graph, no datasets argument
		a.tpm = ['a_tpm']
		assert a.get_tpm_cols() == ['a_tpm']

		# 3: > one dataset in graph, no datasets argument
		a.tpm = ['a_tpm', 'b_tpm']
		assert a.get_tpm_cols() == ['a_tpm', 'b_tpm']

		# 4: empty graph, datasets argument
		# also tests ability to handle char input
		a = swan.SwanGraph()
		with pytest.raises(Exception) as excinfo:
			a.get_tpm_cols('a')
		assert 'Abundance for dataset a' in str(excinfo.value)

		# 5: graph with one dataset and datasets argument already in graph
		a.tpm = ['a_tpm']
		a.counts = ['a_counts']
		assert a.get_tpm_cols('a') == ['a_tpm']

		# 6: graph with one dataset and datasets argument not in graph
		with pytest.raises(Exception) as excinfo:
			a.get_tpm_cols('b')
		assert 'Abundance for dataset b' in str(excinfo.value)

		# 7: graph with more than one dataset and datasets argument in graph
		a.tpm = ['a_tpm', 'b_tpm']
		a.counts = ['a_counts', 'b_counts']
		assert a.get_tpm_cols(['a', 'b']) == ['a_tpm', 'b_tpm']

		# 8: graph with more than one dataset and not all datasets argument in graph
		with pytest.raises(Exception) as excinfo:
			a.get_tpm_cols(['a', 'b', 'c'])
		assert 'Abundance for dataset c'

	# tests get_strand_from_tid
	def test_get_strand_from_tid(self):
		a = get_dummy_sg()
		a.loc_df.loc[1, 'strand'] = '-'

		print(a.loc_df.head())

		assert a.get_strand_from_tid(0) == '+'
		assert a.get_strand_from_tid(1) == '-'
		assert a.get_strand_from_tid(2) == '+'

	# tests get_path_from_tid
	def test_get_path_from_tid(self):
		a = get_dummy_sg()

		assert a.get_path_from_tid(0) == [0,1]
		assert a.get_path_from_tid(1) == [1,2]
		assert a.get_path_from_tid(2) == [0,1,2]

	# tests get_gid_from_tid
	def test_get_gid_from_tid(self):
		a = get_dummy_sg()

		a.t_df.loc[1, 'gid'] = 3

		assert a.get_gid_from_tid(0) == 0
		assert a.get_gid_from_tid(1) == 3
		assert a.get_gid_from_tid(2) == 0

	# tests get_gene_min_max
	def test_get_gene_min_max(self):
		a = get_dummy_sg()
		assert a.get_gene_min_max(0) == (1, 3)

	# tests get_transcript_min_max
	def test_get_transcript_min_max(self):
		a = get_dummy_sg()
		assert a.get_transcript_min_max(0) == (1,3)
		assert a.get_transcript_min_max(1) == (2,3)
		assert a.get_transcript_min_max(2) == (1,2)

	# tests subset_on_gene
	def test_subset_on_gene(self):
		gid = 0
		a = swan.SwanGraph()
		a.t_df = pd.DataFrame({'tid': [0,1,2,3,4,5],
							   'gid': [0,0,1,1,2,2],
							   'path': [[0,1,2],
							   			[2,3,4],
							   			[5,6,7],
							   			[6,7,8],
							   			[9,10,11,12],
							   			[9,11,12]]})
		a.loc_df = pd.DataFrame({'vertex_id': [0,1,2,3,4,5,6,7,8,9,10,11,12], 
								 'strand': ['+','+','+','+','+','-','-','-','-','-','-','-','-'],
								 'chrom': [1,1,1,1,1,2,2,2,2,2,2,2,2],
								 'coord': [0,1,2,3,4,5,6,7,8,9,10,11,12]})
		a.edge_df = pd.DataFrame({'edge_id': [(0,1),(1,2),(2,3),(3,4),(5,6),(6,7),
											  (7,8),(9,10),(10,11),(11,12)]})
		a.t_df = swan.create_dupe_index(a.t_df, 'tid')
		a.t_df = swan.set_dupe_index(a.t_df, 'tid')
		a.loc_df = swan.create_dupe_index(a.loc_df, 'vertex_id')
		a.loc_df = swan.set_dupe_index(a.loc_df, 'vertex_id')
		a.edge_df = swan.create_dupe_index(a.edge_df, 'edge_id')
		a.edge_df = swan.set_dupe_index(a.edge_df, 'edge_id')

		# check subsetting for gene 0 
		a = swan.subset_on_gene(a, 0)

		test = a.t_df['tid'].tolist() 
		control = [0,1]
		check_pairs(control, test)

		test = a.loc_df['vertex_id'].tolist()
		control = [0,1,2,3,4]
		check_pairs(control, test)

		test = a.edge_df['edge_id'].tolist()
		control = [(0,1),(1,2),(2,3),(3,4)]
		check_pairs(control, test)

##########################################################################
############################# Pytest utilities ###########################
##########################################################################

def check_pairs(control, test):
	print('control')
	print(control)
	for t in test:
		print(t)
		assert t in control
	assert len(test) == len(control)

def get_dummy_sg(special=None):
	a = swan.SwanGraph()

	loc_df = pd.DataFrame({'chrom': [1,1,1],
		'coord': [1,3,2],
		'strand': ['+', '+', '+'],
		'vertex_id': [0,1,2]})
	loc_df = swan.create_dupe_index(loc_df, 'vertex_id')
	loc_df = swan.set_dupe_index(loc_df, 'vertex_id')

	edge_df = pd.DataFrame({'edge_id': [(0,2),(0,1),(1,2)],
			  'v1': [0,0,1],
			  'v2': [2,1,2],
			  'edge_type': ['exon', 'exon', 'intron'], 
			  'strand': ['+','+','+']})
	edge_df = swan.create_dupe_index(edge_df, 'edge_id')
	edge_df = swan.set_dupe_index(edge_df, 'edge_id')

	t_df = pd.DataFrame({'tid':[2,1,0],
			  'gid':[0,0,0],
			  'gname':['0','0','0'],
			  'path':[[0,1,2],[1,2],[0,1]],
			  'counts_a': [0,0,12],
			  'counts_b': [1,0,14]})
	t_df = swan.create_dupe_index(t_df, 'tid')
	t_df = swan.set_dupe_index(t_df, 'tid')

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






