import pytest
import sys
import numpy as np
sys.path.append('../utils/')
sys.path.append('../../refactor_splice_graph/')
import SpliceGraph as sg
import PlottedGraph as pg
from utils import *
from plotting_tools import * 

class TestMergeSGs(object):
	def test_merge_loc_dfs(self):
		""" 
		
		"""
		a = pd.DataFrame({'chrom': [1,2,3],
			'coord': [1,1,1],
			'strand': ['+', '-', '+'],
			'vertex_id': [0,1,2]})
		b = pd.DataFrame({'chrom': [1,2,4],
			'coord': [1,1,1],
			'strand': ['+', '-', '+'],
			'vertex_id': [0,1,2]})
		loc_df = sg.merge_loc_dfs(a,b,'dataset_a','dataset_b')
		id_map = sg.get_vertex_id_map(loc_df, 'dataset_a', 'dataset_b')
		loc_df = sg.assign_new_vertex_ids(loc_df, id_map)

		print(loc_df)
		assert len(loc_df.index) == 4

		# make sure that locations are being assigned to correct datasets
		loc_dataset_pair_a = [(loc_id,d_a)
			for loc_id,d_a in zip(loc_df.vertex_id.tolist(),loc_df.dataset_a.tolist())]
		print('loc dataset a pair')
		print(loc_dataset_pair_a)
		control = [(0,True), (1,True), (2,True), (3,False)]
		check_pairs(control, loc_dataset_pair_a)

		loc_dataset_pair_b = [(loc_id,d_b)
			for loc_id,d_b in zip(loc_df.vertex_id.tolist(),loc_df.dataset_b.tolist())]
		print('loc dataset b pair')
		print(loc_dataset_pair_b)
		control = [(0,True), (1,True), (2,False), (3,True)]
		check_pairs(control, loc_dataset_pair_b)	

	def test_get_vertex_id_map(self):
		df = pd.DataFrame({'chrom': [1,2,3,4],
			'coord': [1,1,1,1],
			'strand': ['+', '-', '+', '+'],
			'vertex_id_a': [0,1,2,np.nan],
			'vertex_id_b': [3,4,np.nan,5],
			'dataset_a': [True,True,True,False],
			'dataset_b': [True,True,False,True]})
		b_map = sg.get_vertex_id_map(df, 'dataset_a', 'dataset_b')
		print(b_map)
		assert b_map == {3:0, 4:1, 5:3}

	def test_assign_new_vertex_ids(self):
		df = pd.DataFrame({'chrom': [1,2,3,4],
			'coord': [1,1,1,1],
			'strand': ['+', '-', '+', '+'],
			'vertex_id_a': [0,1,2,np.nan],
			'vertex_id_b': [3,4,np.nan,5],
			'present_in': [['a','b'],['a','b'],['a'],['b']]})
		b_map = {5:3}

		df = sg.assign_new_vertex_ids(df, b_map)
		test = df.vertex_id.tolist()
		print(df)
		print('vertex ids: {}'.format(test))
		assert test == [0,1,2,3]

	def test_assign_new_edge_ids(self):
		loc_df = pd.DataFrame({'chrom': [1,2,3,4], 
							   'coord': [1,1,1,1],
							   'strand': ['+', '+', '+', '+'],
							   'vertex_id_a': [0,1,2,np.nan],
							   'vertex_id_b': [3,4,np.nan,5],
							   'present_in': [['a','b'],['a','b'],['a'],['b']],
							   'vertex_id': [0,1,2,3]})
		b = pd.DataFrame({'edge_id': [(3,5),(3,4),(4,5)],
						  'v1': [3,3,4],
						  'v2': [5,4,5],
						  'edge_type': ['exon', 'exon', 'intron'],
						  'strand': ['+','+','+']})
		id_map = {3:0, 4:1, 5:3}
		b = sg.assign_new_edge_ids(b, id_map)
		print('b after replacing entries')
		print(b)

		assert b.v1.tolist() == [0,0,1]
		assert b.v2.tolist() == [3,1,3]
		assert b.edge_id.tolist() == [(0,3),(0,1),(1,3)]

	def test_merge_edge_dfs(self):
		loc_df = pd.DataFrame({'chrom': [1,2,3,4], 
							   'coord': [1,1,1,1],
							   'strand': ['+', '+', '+', '+'],
							   'vertex_id_a': [0,1,2,np.nan],
							   'vertex_id_b': [3,4,np.nan,5],
							   'present_in': [['a','b'],['a','b'],['a'],['b']],
							   'vertex_id': [0,1,2,3]})
		a = pd.DataFrame({'edge_id': [(0,2),(0,1),(1,2)],
						  'v1': [0,0,1],
						  'v2': [2,1,2],
						  'edge_type': ['exon', 'exon', 'intron'], 
						  'strand': ['+','+','+']})
		b = pd.DataFrame({'edge_id': [(0,3),(0,1),(1,3)],
						  'v1': [0,0,1],
						  'v2': [3,1,3],
						  'edge_type': ['exon', 'exon', 'intron'],
						  'strand': ['+','+','+']})
		edge_df = sg.merge_edge_dfs(a, b, 'dataset_a', 'dataset_b')

		# test that correct edge ids assigned
		assert set(edge_df.edge_id.tolist()) == set([(0,2),(0,1),(1,2),(0,3),(1,3)])

		# check that edge ids were indicated as present in the right dfs
		id_dataset_a_pairs = [(e_id, d_a)
				for e_id,d_a in zip(edge_df.edge_id.tolist(),edge_df.dataset_a.tolist())]
		print('id dataset a pairs')
		print(id_dataset_a_pairs)
		control_pairs = [((0,2),True),((0,1),True),((1,2),True),((0,3),False),((1,3),False)]
		check_pairs(control_pairs, id_dataset_a_pairs)

		# check that edge ids were indicated as present in the right dfs
		print(edge_df)
		print(edge_df.edge_id.tolist())
		print(edge_df.dataset_b.tolist())
		id_dataset_b_pairs = [(e_id, d_b)
				for e_id,d_b in zip(edge_df.edge_id.tolist(),edge_df.dataset_b.tolist())]
		print('id dataset b pairs')
		print(id_dataset_b_pairs)
		control_pairs = [((0,2),False),((0,1),True),((1,2),False),((0,3),True),((1,3),True)]
		check_pairs(control_pairs, id_dataset_b_pairs)

	def test_assign_new_paths(self):
		id_map = {3:0, 4:1, 5:3}
		b = pd.DataFrame({'tid':[4,2,3],
				  'gid':[0,0,0],
				  'gname':['0','0','0'],
				  'path':[[3,4,5],[3,5],[3,4]]})
		b = sg.assign_new_paths(b, id_map)

		# check that ids in paths are being correctly reassigned
		print(b)
		assert set(tuple(i) for i in b.path.tolist()) == set([(0,1,3),(0,3),(0,1)])

	def test_merge_t_dfs(self):
		id_map = {3:0, 4:1, 5:3}
		a = pd.DataFrame({'tid':[0,1,3],
						  'gid':[0,0,0],
						  'gname':['0','0','0'],
						  'path':[[0,1,2],[0,2],[0,1]]})
		b = pd.DataFrame({'tid':[4,2,3],
				  'gid':[0,0,0],
				  'gname':['0','0','0'],
				  'path':[[0,1,3],[0,3],[0,1]]})
		t_df = sg.merge_t_dfs(a,b,'dataset_a','dataset_b')

		# this will make sure that tids have the same path?
		paths = [tuple(path) for path in t_df.path.tolist()]
		tid_path_pairs = [(tid,path) for tid,path in zip(t_df.tid.tolist(),paths)]
		print('tid_path_pairs')
		print(tid_path_pairs)
		control_pairs = [(0,(0,1,2)),(1,(0,2)),(2,(0,3)),(3,(0,1)),(4,(0,1,3))]
		check_pairs(control_pairs, tid_path_pairs)

		# test that these were correctly identified as from each dataset
		tid_dataset_a_pairs = [(tid,d_a)
			for tid,d_a in zip(t_df.tid.tolist(),t_df.dataset_a.tolist())]
		print('tid dataset a pairs')
		print(tid_dataset_a_pairs)
		control_pairs = [(0,True),(1,True),(2,False),(4,False),(3,True)]
		check_pairs(control_pairs, tid_dataset_a_pairs)

		tid_dataset_b_pairs = [(tid,d_b)
			for tid,d_b in zip(t_df.tid.tolist(),t_df.dataset_b.tolist())]
		print('tid dataset b pairs')
		print(tid_dataset_b_pairs)
		control_pairs = [(0,False),(1,False),(2,True),(4,True),(3,True)]
		check_pairs(control_pairs, tid_dataset_b_pairs)

def check_pairs(control, test):
	print('control')
	print(control)
	for t in test:
		print(t)
		assert t in control






