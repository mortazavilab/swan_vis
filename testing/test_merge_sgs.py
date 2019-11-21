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
		assert len(merged.index) == 4

	def test_assign_locs_present(self):
		df = pd.DataFrame({'chrom': [1,2,3,4],
			'coord': [1,1,1,1],
			'strand': ['+', '-', '+','+'],
			'vertex_id_a': [0,1,2,np.nan],
			'vertex_id_b': [0,1,np.nan,2]})

		df['present_in'] = df.apply(lambda x: sg.loc_present_in(x), axis=1)

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
		assert b_map == {5:3}

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
		b = sg.assign_new_edge_ids(loc_df, b)
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
		b = pd.DataFrame({'edge_id': [(3,5),(3,4),(4,5)],
						  'v1': [3,3,4],
						  'v2': [5,4,5],
						  'edge_type': ['exon', 'exon', 'intron'],
						  'strand': ['+','+','+']})
		print('a')
		print(a)
		print('b')
		print(b)
		edge_df = sg.merge_edge_dfs(loc_df, a, b)
		print('edge_df after merging')
		print(edge_df)
		assert 1==2





