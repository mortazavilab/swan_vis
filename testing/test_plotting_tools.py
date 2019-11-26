import pytest
import sys
import pandas as pd
import numpy as np
sys.path.append('../utils/')
sys.path.append('../../refactor_splice_graph/')
import SpliceGraph as sg
import PlottedGraph as pg
from utils import *
from plotting_tools import * 
import networkx as nx

class TestPlottingTools(object):
	def test_get_gene_min_max(self):
		loc_df = pd.DataFrame({'vertex_id': [0,1,2,3,4],
							   'coord': [1,2,3,4,5]})
		loc_df = create_dupe_index(loc_df, 'vertex_id')
		loc_df = set_dupe_index(loc_df, 'vertex_id')

		t_df = pd.DataFrame({'tid': [0,1,2],
							 'gid': [0,0,1],
							 'path': [[0,1], [0,2], [3,4]]})
		t_df = create_dupe_index(t_df, 'tid')
		t_df = set_dupe_index(t_df, 'tid')

		print(loc_df)
		print(t_df)

		assert sg.get_gene_min_max(loc_df, t_df, 0) == (1,3)
		assert sg.get_gene_min_max(loc_df, t_df, 1) == (4,5)

	def test_get_coord_map(self):
		loc_df = pd.DataFrame({'vertex_id': [0,1,2,3,4],
							   'coord': [1,2,3,2,3],
							   'chrom': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
							   'strand': ['+', '+', '+', '-', '-']})
		loc_df = create_dupe_index(loc_df, 'vertex_id')
		loc_df = set_dupe_index(loc_df, 'vertex_id')

		# gene 0
		coord_map = get_coord_map(1,3, loc_df, 'chr1', '+')
		assert len(coord_map.index) == 3

		# gene 1
		coord_map = get_coord_map(2,3, loc_df, 'chr1', '-')
		assert len(coord_map.index) == 2


