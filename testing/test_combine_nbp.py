
import pytest
import sys
import numpy as np
sys.path.append('../utils/')
sys.path.append('../../refactor_splice_graph/')
from SpliceGraph import SpliceGraph
from PlottedGraph import PlottedGraph
from utils import *
import networkx as nx
import math

class TestCombineNBP(object):

	# graph where not all nodes must be combined
	def test_not_all_nodes_combined(self):
		sg = SpliceGraph()
		sg.loc_df = pd.DataFrame({'vertex_id':[0,1,2,3,4,5],
								  'coord': [0,1,2,3,4,5],
								  'strand': ['+','+','+','+','+','+'],
								  'chrom': [1,1,1,1,1,1]})
		sg.edge_df = pd.DataFrame({'edge_id':[(0,1),(1,2),(2,3),(2,4),(3,5),(4,5)],
								   'v1': [0,1,2,2,3,4],
								   'v2': [1,2,3,4,5,5],
								   'strand':['+','+','+','+','+','+'],
								   'edge_type': ['exon', 'intron', 'exon', 'exon', 'intron', 'intron']})
		sg.t_df = pd.DataFrame({'tid': [0,1],
								'gname': ['0','0'],
								'gid': [0,0],
								'path': [[0,1,2,3,5],[0,1,2,4,5]]})
		sg.loc_df = create_dupe_index(sg.loc_df, 'vertex_id')
		sg.edge_df = create_dupe_index(sg.edge_df, 'edge_id')
		sg.t_df = create_dupe_index(sg.t_df, 'tid')
		sg.loc_df = set_dupe_index(sg.loc_df, 'vertex_id')
		sg.edge_df = set_dupe_index(sg.edge_df, 'edge_id')                                         
		sg.t_df = set_dupe_index(sg.t_df, 'tid')
		sg.get_loc_types()
		sg.create_graph_from_dfs()

		sg.pg = PlottedGraph(sg, True, False, False)

		# first make sure that only the correct edges/nodes were 
		# kept in loc_df and edge_df
		assert set(sg.pg.loc_df.vertex_id.tolist()) == set(['c0',2,3,4,5])
		assert set(sg.pg.edge_df.edge_id.tolist()) == set([('c0',2),(2,3),(2,4),(3,5),(4,5)])

		# make sure that the correct nbps were assigned to each node
		assert tuple(sg.pg.loc_df.loc['c0', 'agg_path']) == (0,1)

		# make sure the paths were updated correctly in t_df
		test = set([tuple(path) for path in sg.pg.t_df.path.tolist()])
		control = [('c0',2,3,5), ('c0',2,4,5)]
		check_pairs(control, test)

	# graph where an NBP is halted because it reaches an alt TSS
	def test_TSS(self):
		sg = SpliceGraph()
		sg.loc_df = pd.DataFrame({'vertex_id':[0,1,2,3,4,5],
								  'coord': [0,1,2,3,4,5],
								  'strand': ['+','+','+','+','+','+'],
								  'chrom': [1,1,1,1,1,1]})
		sg.edge_df = pd.DataFrame({'edge_id': [(0,1),(1,2),(2,3),(3,4),(4,5)],
								   'v1': [0,1,2,3,4],
								   'v2': [1,2,3,4,5],
								   'strand': ['+','+','+','+','+'],
								   'edge_type': ['exon','intron','exon','intron','exon']})
		sg.t_df = pd.DataFrame({'tid': [0,1],
								'gname': ['0','0'],
								'gid': [0,0],
								'path': [[0,1,2,3,4,5],[2,3,4,5]]})
		sg.loc_df = create_dupe_index(sg.loc_df, 'vertex_id')
		sg.edge_df = create_dupe_index(sg.edge_df, 'edge_id')
		sg.t_df = create_dupe_index(sg.t_df, 'tid')
		sg.loc_df = set_dupe_index(sg.loc_df, 'vertex_id')
		sg.edge_df = set_dupe_index(sg.edge_df, 'edge_id')                                         
		sg.t_df = set_dupe_index(sg.t_df, 'tid')
		sg.get_loc_types()
		sg.create_graph_from_dfs()

		sg.pg = PlottedGraph(sg, True, False, False)

		# first make sure that only the correct edges/nodes were 
		# kept in loc_df and edge_df
		assert set(sg.pg.loc_df.vertex_id.tolist()) == set(['c0','c1'])
		assert set(sg.pg.edge_df.edge_id.tolist()) == set([('c0','c1')])

		# make sure that the correct nbps were assigned to each node
		assert tuple(sg.pg.loc_df.loc['c0', 'agg_path']) == (0,1)
		assert tuple(sg.pg.loc_df.loc['c1', 'agg_path']) == (2,3,4,5)

		# make sure the paths were updated correctly in t_df
		test = set([tuple(path) for path in sg.pg.t_df.path.tolist()])
		control = [('c0','c1'), ('c1',)]
		check_pairs(control, test)

	# graph where an NBP is halted because it reaches an alt TES
	def test_TES(self):
		sg = SpliceGraph()
		sg.loc_df = pd.DataFrame({'vertex_id':[0,1,2,3,4,5],
								  'coord': [0,1,2,3,4,5],
								  'strand': ['+','+','+','+','+','+'],
								  'chrom': [1,1,1,1,1,1]})
		sg.edge_df = pd.DataFrame({'edge_id': [(0,1),(1,2),(2,3),(3,4),(4,5)],
								   'v1': [0,1,2,3,4],
								   'v2': [1,2,3,4,5],
								   'strand': ['+','+','+','+','+'],
								   'edge_type': ['exon','intron','exon','intron','exon']})
		sg.t_df = pd.DataFrame({'tid': [0,1],
								'gname': ['0','0'],
								'gid': [0,0],
								'path': [[0,1,2,3,4,5],[0,1,2,3]]})
		sg.loc_df = create_dupe_index(sg.loc_df, 'vertex_id')
		sg.edge_df = create_dupe_index(sg.edge_df, 'edge_id')
		sg.t_df = create_dupe_index(sg.t_df, 'tid')
		sg.loc_df = set_dupe_index(sg.loc_df, 'vertex_id')
		sg.edge_df = set_dupe_index(sg.edge_df, 'edge_id')                                         
		sg.t_df = set_dupe_index(sg.t_df, 'tid')
		sg.get_loc_types()
		sg.create_graph_from_dfs()

		sg.pg = PlottedGraph(sg, True, False, False)

		# first make sure that only the correct edges/nodes were 
		# kept in loc_df and edge_df
		assert set(sg.pg.loc_df.vertex_id.tolist()) == set(['c0','c1'])
		assert set(sg.pg.edge_df.edge_id.tolist()) == set([('c0','c1')])

		# make sure that the correct nbps were assigned to each node
		assert tuple(sg.pg.loc_df.loc['c0', 'agg_path']) == (0,1,2,3)
		assert tuple(sg.pg.loc_df.loc['c1', 'agg_path']) == (4,5)

		# make sure the paths were updated correctly in t_df
		test = set([tuple(path) for path in sg.pg.t_df.path.tolist()])
		control = [('c0','c1'), ('c0',)]
		check_pairs(control, test)


	# more complicated example with several things being tested
	# nbp interrupted by TSS
	# nbp interrupted by TES
	# all nodes lumped together
	def test_combine(self):
		sg = SpliceGraph()
		sg.loc_df = pd.DataFrame({'vertex_id':[0,1,2,3,4,5,6,7,8,9,10,11],
								  'coord': [0,1,2,3,4,5,6,7,8,9,10,11],
								  'strand': ['+','+','+','+','+','+','+','+','+','+','+','+'],
								  'chrom': [1,1,1,1,1,1,1,1,1,1,1,1]})
		sg.edge_df = pd.DataFrame({'edge_id': [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,10),(10,11),(8,9),(9,10)],
								   'v1': [0,1,2,3,4,5,6,7,10,8,9],
								   'v2': [1,2,3,4,5,6,7,10,11,9,10],
								   'strand': ['+','+','+','+','+','+','+','+','+','+','+'],
								   'edge_type': ['exon','intron','exon','intron','exon','intron','exon','intron','exon','exon','intron']})
		sg.t_df = pd.DataFrame({'tid': [0,1,2],
								'gname': ['0','0','0'],
								'gid': [0,0,0],
								'path': [[0,1,2,3,4,5,6,7],[2,3,4,5,6,7,10,11],[8,9,10,11]]})
		sg.loc_df = create_dupe_index(sg.loc_df, 'vertex_id')
		sg.edge_df = create_dupe_index(sg.edge_df, 'edge_id')
		sg.t_df = create_dupe_index(sg.t_df, 'tid')
		sg.loc_df = set_dupe_index(sg.loc_df, 'vertex_id')
		sg.edge_df = set_dupe_index(sg.edge_df, 'edge_id')                                         
		sg.t_df = set_dupe_index(sg.t_df, 'tid')
		sg.get_loc_types()
		sg.create_graph_from_dfs()
		sg.pg = PlottedGraph(sg, True, False, False)

		# first make sure that only the correct edges/nodes were 
		# kept in loc_df and edge_df
		print(set(sg.pg.loc_df.vertex_id.tolist()))
		assert set(sg.pg.loc_df.vertex_id.tolist()) == set(['c0','c1','c2','c3'])
		print(set(sg.pg.edge_df.edge_id.tolist()))
		assert set(sg.pg.edge_df.edge_id.tolist()) == set([('c0','c1'),('c1','c3'),('c2','c3')])

		# make sure that the correct nbps were assigned to each node
		assert tuple(sg.pg.loc_df.loc['c0', 'agg_path']) == (0,1)
		assert tuple(sg.pg.loc_df.loc['c1', 'agg_path']) == (2,3,4,5,6,7)
		assert tuple(sg.pg.loc_df.loc['c2', 'agg_path']) == (8,9)
		assert tuple(sg.pg.loc_df.loc['c3', 'agg_path']) == (10,11)

		# make sure the paths were updated correctly in t_df
		test = set([tuple(path) for path in sg.pg.t_df.path.tolist()])
		control = [('c0','c1'), ('c1','c3'), ('c2','c3')]
		check_pairs(control, test)




def check_pairs(control, test):
	print('control')
	print(control)
	for t in test:
		print(t)
		assert t in control
	assert len(test) == len(control)


