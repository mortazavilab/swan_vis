import pytest
import sys
import numpy as np
sys.path.append('../utils/')
sys.path.append('../../refactor_splice_graph/')
import SpliceGraph as sg
import PlottedGraph as pg
from utils import *
from plotting_tools import * 
import networkx as nx

class TestSpliceGraph(object):
	def test_get_loc_types(self):
		t_df = pd.DataFrame({'tid':[0,1,2,4,3],
							 'gid':[0,0,0,0,0],
							 'gname':['0','0','0','0','0'],
							 'path':[[0,1,2],[0,2],[0,3],[0,1,3],[0,1]]})
		loc_df = pd.DataFrame({'vertex_id':[0,1,2,3]})
		loc_df = create_dupe_index(loc_df, 'vertex_id')
		loc_df = set_dupe_index(loc_df, 'vertex_id')
		loc_df = sg.get_loc_types(loc_df, t_df)

		c_tss = [0]
		c_alt_tss = []
		c_tes = [1,2,3]
		c_alt_tes = [1,2,3]
		c_internal = [1]

		tss = loc_df[loc_df.apply(lambda x: True if x.TSS else False, axis=1)].vertex_id.tolist()
		alt_tss = loc_df[loc_df.apply(lambda x: True if x.alt_TSS else False, axis=1)].vertex_id.tolist()
		tes = loc_df[loc_df.apply(lambda x: True if x.TES else False, axis=1)].vertex_id.tolist()
		alt_tes = loc_df[loc_df.apply(lambda x: True if x.alt_TES else False, axis=1)].vertex_id.tolist()
		internal = loc_df[loc_df.apply(lambda x: True if x.internal else False, axis=1)].vertex_id.tolist()

		print(tss)
		assert c_tss == tss
		print(alt_tss)
		assert c_alt_tss == alt_tss
		print(tes)
		assert c_tes == tes
		print(alt_tes)
		assert c_alt_tes == alt_tes
		print(internal)
		assert c_internal == internal

	def test_get_dataset_fields(self):
		G = nx.Graph()
		G.add_node(0)
		attrs = {0: {'coord': 3, 'internal': True, 'TES': False,
				  'dataset_1': True, 'dataset_lizard': True,
				  'dataset_-1': False}}
		nx.set_node_attributes(G, attrs)
		print(G.nodes(data=True))
		fields = sg.get_dataset_fields(G)
		print('fields: ')
		print(fields)
		c_fields = ['dataset_1', 'dataset_lizard', 'dataset_-1']

		assert set(fields) == set(c_fields)

	def test_add_dataset(self):
		loc_df = pd.DataFrame({'chrom': [1,2,3],
			'coord': [1,1,1],
			'strand': ['+', '-', '+'],
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

		t_df = pd.DataFrame({'tid':[0,1,3],
				  'gid':[0,0,0],
				  'gname':['0','0','0'],
				  'path':[[0,1,2],[0,2],[0,1]]})
		t_df = create_dupe_index(t_df, 'tid')
		t_df = set_dupe_index(t_df, 'tid')

		loc_df = sg.get_loc_types(loc_df, t_df)

		splice_graph = sg.SpliceGraph(loc_df=loc_df, edge_df=edge_df, t_df=t_df)
		splice_graph.add_dataset('annotation')
		print(splice_graph.loc_df.head())
		print(splice_graph.edge_df.head())
		print(splice_graph.t_df.head())
		print(splice_graph.G.nodes(data=True))
		print(splice_graph.G.edges(data=True))
		
		# df cols added correctly
		assert all(d == True for d in splice_graph.loc_df.dataset_annotation)
		assert all(d == True for d in splice_graph.edge_df.dataset_annotation)
		assert all(d == True for d in splice_graph.t_df.dataset_annotation)

		# node/edge designations added correctly
		assert all(i[1] == True for i in splice_graph.G.nodes(data='dataset_annotation'))
		assert all(i[2] == True for i in splice_graph.G.edges(data='dataset_annotation'))

	def test_remove_dataset(self):
		loc_df = pd.DataFrame({'chrom': [1,2,3],
			'coord': [1,1,1],
			'strand': ['+', '-', '+'],
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

		t_df = pd.DataFrame({'tid':[0,1,3],
				  'gid':[0,0,0],
				  'gname':['0','0','0'],
				  'path':[[0,1,2],[0,2],[0,1]]})
		t_df = create_dupe_index(t_df, 'tid')
		t_df = set_dupe_index(t_df, 'tid')

		loc_df = sg.get_loc_types(loc_df, t_df)

		splice_graph = sg.SpliceGraph(loc_df=loc_df, edge_df=edge_df, t_df=t_df)
		splice_graph.add_dataset('annotation')

		splice_graph.remove_dataset('annotation')

		print(splice_graph.loc_df.head())
		print(splice_graph.edge_df.head())
		print(splice_graph.t_df.head())

		assert 'dataset_annotation' not in splice_graph.loc_df.columns
		assert 'dataset_annotation' not in splice_graph.edge_df.columns
		assert 'dataset_annotation' not in splice_graph.t_df.columns

	def test_add_abundance_dataset(self):
		loc_df = pd.DataFrame({'chrom': [1,2,3],
			'coord': [1,1,1],
			'strand': ['+', '-', '+'],
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

		t_df = pd.DataFrame({'tid':[0,1,2],
				  'gid':[0,0,0],
				  'gname':['0','0','0'],
				  'path':[[0,1,2],[0,2],[0,1]]})
		t_df = create_dupe_index(t_df, 'tid')
		t_df = set_dupe_index(t_df, 'tid')

		loc_df = sg.get_loc_types(loc_df, t_df)

		splice_graph = sg.SpliceGraph(loc_df=loc_df, edge_df=edge_df, t_df=t_df)

		# add dataset test
		file = 'input_files/test_abundance.tsv'
		splice_graph.add_abundance_dataset(file, ['count_1a', 'count_1b'], '1')
		splice_graph.add_abundance_dataset(file, ['count_2a', 'count_2b'], '2')

		print(splice_graph.t_df)

		# dataset 1 pairs
		dataset_1_pairs = splice_graph.t_df.apply(lambda x: (x.tid, x.counts_1), axis=1)
		control_pairs = [(0,4),(1,8),(2,12)]
		check_pairs(control_pairs, dataset_1_pairs)

		# dataset 2 pairs
		dataset_2_pairs = splice_graph.t_df.apply(lambda x: (x.tid, x.counts_2), axis=1)
		control_pairs = [(0,6),(1,10),(2,14)]
		check_pairs(control_pairs, dataset_2_pairs)

def check_pairs(control, test):
	print('control')
	print(control)
	for t in test:
		print(t)
		assert t in control





