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
		merged = sg.merge_loc_dfs(a,b)

		print(merged)
		# TODO this test case isn't good enough
		assert len(merged.index) == 4

	def test_assign_locs_present(self):
		df = pd.DataFrame({'chrom': [1,2,3,4],
			'coord': [1,1,1,1],
			'strand': ['+', '-', '+','+'],
			'vertex_id_a': [0,1,2,np.nan],
			'vertex_id_b': [0,1,np.nan,2],
			'present_in_a':[True,True,True,False],
			'present_in_b':[True,True,False,True]})

		df['present_in'] = df.apply(lambda x: sg.present_in(x), axis=1)

		print(df)

		p = df.loc[df.vertex_id_a == 0, 'present_in'].tolist()[0]
		print(p)
		assert p == ['a', 'b']

		p = df.loc[df.vertex_id_a == 1, 'present_in'].tolist()[0]
		print(p)
		assert p == ['a', 'b']

		p = df.loc[df.vertex_id_a == 2, 'present_in'].tolist()[0]
		print(p)
		assert p == ['a']

		p = df.loc[df.vertex_id_b == 2, 'present_in'].tolist()[0]
		print(p)
		assert p == ['b']

	def test_get_vertex_id_map(self):
		df = pd.DataFrame({'chrom': [1,2,3,4],
			'coord': [1,1,1,1],
			'strand': ['+', '-', '+', '+'],
			'vertex_id_a': [0,1,2,np.nan],
			'vertex_id_b': [3,4,np.nan,5],
			'present_in': [['a','b'],['a','b'],['a'],['b']]})
		b_map = sg.get_vertex_id_map(df)
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
		edge_df = sg.merge_edge_dfs(a, b)

		# test that correct edge ids assigned
		assert set(edge_df.edge_id.tolist()) == set([(0,2),(0,1),(1,2),(0,3),(1,3)])

		# check that edge ids were indicated as present in the right dfs
		id_present_pairs = [(e_id, tuple(p))
				for e_id,p in zip(edge_df.edge_id.tolist(),edge_df.present_in.tolist())]
		print(id_present_pairs)
		control_pairs = [((0,2),('a',)),((0,1),('a','b')),((1,2),('a',)),((0,3),('b',)),((1,3),('b',))]
		check_pairs(control_pairs, id_present_pairs)

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
		t_df = sg.merge_t_dfs(a,b)

		# this will make sure that tids have the same path?
		paths = [tuple(path) for path in t_df.path.tolist()]
		tid_path_pairs = [(tid,path) for tid,path in zip(t_df.tid.tolist(),paths)]
		print('tid_path_pairs')
		print(tid_path_pairs)
		control_pairs = [(0,(0,1,2)),(1,(0,2)),(2,(0,3)),(3,(0,1)),(4,(0,1,3))]
		check_pairs(control_pairs, tid_path_pairs)

		# test that these were correctly identified as from each dataset
		tid_present_pairs = [(tid,tuple(p))
			for tid,p in zip(t_df.tid.tolist(),t_df.present_in.tolist())]
		print(tid_present_pairs)
		control_pairs = [(0,('a',)),(1,('a',)),(2,('b',)),(4,('b',)),(3,('a','b'))]
		check_pairs(control_pairs, tid_present_pairs)

def check_pairs(control, test):
	print('control')
	print(control)
	for t in test:
		print(t)
		assert t in control






