import pytest
import sys
import numpy as np
import swan_vis as swan

class TestAddDatasets(object):
	def test_add_annotation(self):

		# adding an annotation to an empty graph
		print('testing for correct novelty assignment when adding annotation')
		sg = swan.SwanGraph()
		sg.add_annotation('input_files/annot.gtf')
		control = len(sg.t_df.index)
		test = len(sg.t_df.loc[sg.t_df.novelty == 'Known'].index)
		assert control == test

		# adding an annotation to a graph that already has data 
		# but preexisting data does not have novelty categories
		print('testing for correct novelty assignment when adding annotation '
			  'to graph with preexisting data that does not contain novelty info')
		sg = swan.SwanGraph()
		sg.add_dataset('a', 'input_files/annot_2.gtf')
		sg.add_annotation('input_files/annot.gtf')
		control = [('ENST01', 'Known'), ('ENST03', 'Known'), ('ENST07', 'Known'),
				   ('ENST02', 'Undefined'), ('ENST04', 'Undefined'), ('ENST08', 'Undefined')]
		test = sg.t_df.apply(lambda x: (x.tid, x.novelty), axis=1)
		check_pairs(control, test)

		# adding an annotation to a graph that already has data 
		# and preexisting data has novelty categories
		print('testing for correct novelty assignment when adding annotation '
			  'to graph with preexisting data that does contain novelty info')
		sg = swan.SwanGraph()
		sg.add_dataset('a', 'input_files/annot_2.gtf')
		sg.t_df['novelty'] = ['ISM', 'NNC', 'NIC', 'NIC']
		sg.add_annotation('input_files/annot.gtf')
		control = [('ENST01', 'Known'), ('ENST03', 'Known'), ('ENST07', 'Known'),
				   ('ENST02', 'NNC'), ('ENST04', 'NIC'), ('ENST08', 'NIC')]
		test = sg.t_df.apply(lambda x: (x.tid, x.novelty), axis=1)
		check_pairs(control, test)

	def test_weird_gtf(self):
		sg = swan.SwanGraph()
		sg.add_dataset('test', 'input_files/weird_gtf_entries.gtf')
		print(sg.t_df)

		# check each transcript 
		tid = 'ENST00000002165.11_3'
		path = sg.t_df.loc[tid, 'path']
		print(path)
		coords = sg.loc_df.loc[path, 'coord'].tolist()
		ctrl_coords = [143832857, 143832548, 143828561,
					   143828374, 143825389, 143825050,
					   143823702, 143823492, 143823259,
					   143823069, 143818634, 143818526, 
					   143816984, 143815949]
		check_pairs(ctrl_coords, coords)
		tid = 'ENST00000514436.1'
		path = sg.t_df.loc[tid, 'path']
		print(path)
		coords = sg.loc_df.loc[path, 'coord'].tolist()
		ctrl_coords = [326096, 326569, 327348, 328112]
		check_pairs(ctrl_coords, coords)

def check_pairs(control, test):
	print('control')
	print(control)
	print('test')
	print(test)
	for t in test:
		assert t in control
	assert len(test) == len(control)