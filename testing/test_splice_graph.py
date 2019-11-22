import pytest
import sys
import numpy as np
sys.path.append('../utils/')
sys.path.append('../../refactor_splice_graph/')
import SpliceGraph as sg
import PlottedGraph as pg
from utils import *
from plotting_tools import * 

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