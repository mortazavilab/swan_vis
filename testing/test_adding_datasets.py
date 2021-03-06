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

	def test_minus_strand_unordered_gtf(self):
		sg = process_gtf()
		print(sg.t_df)

		# check transcript 
		tid = 'ENST00000002165.11_3'
		path = sg.t_df.loc[tid, 'path']
		print(path)
		coords = sg.loc_df.loc[path, 'coord'].tolist()
		ctrl_coords = [143832857, 143832548, 143828561,
					   143828374, 143825389, 143825050,
					   143823702, 143823492, 143823259,
					   143823069, 143818634, 143818526, 
					   143816984, 143815949]
		assert ctrl_coords == coords
		edge_ids = [(path[i],path[i+1]) for i in range(len(path)-1)]
		for eid in edge_ids:
			print(eid)
			assert eid in sg.edge_df.index.tolist()

	def test_plus_strand_unordered_gtf(self):
		sg = process_gtf()
		print(sg.t_df)

		# check transcript
		tid = 'ENST00000514436.1'
		path = sg.t_df.loc[tid, 'path']
		print(path)
		coords = sg.loc_df.loc[path, 'coord'].tolist()
		ctrl_coords = [326096, 326569, 327348, 328112]
		assert ctrl_coords == coords
		edge_ids = [(path[i],path[i+1]) for i in range(len(path)-1)]
		for eid in edge_ids:
			print(eid)
			assert eid in sg.edge_df.index.tolist()

	def test_unordered_gtf_numeric_coord_sort(self):
		sg = process_gtf()
		print(sg.t_df)

		# check transcript
		tid = 'ENST00000620134.4'
		path = sg.t_df.loc[tid, 'path']
		print(path)
		coords = sg.loc_df.loc[path, 'coord'].tolist()
		ctrl_coords = [138860, 138642, 138334, 138150,
			130591, 130522, 119255, 119126, 117375, 117301,
			112775, 112622, 110606, 110525, 100509, 100372, 
			98301, 98145, 93325, 93219, 92725, 92596, 89902,
			89713, 88890, 88698, 86870, 85820, 84595, 84274,]
		assert ctrl_coords == coords
		edge_ids = [(path[i],path[i+1]) for i in range(len(path)-1)]
		for eid in edge_ids:
			print(eid)
			assert eid in sg.edge_df.index.tolist()

	def test_weird_gtf_weird_chrom(self):
		sg = process_gtf()
		print(sg.t_df)

		# check transcript
		tid = 'TALONT000206784'
		path = sg.t_df.loc[tid, 'path']
		print(path)
		coords = sg.loc_df.loc[path, 'coord'].tolist()
		ctrl_coords = [115043, 114986, 112850, 112792, 62949, 61468]
		assert ctrl_coords == coords
		edge_ids = [(path[i],path[i+1]) for i in range(len(path)-1)]
		for eid in edge_ids:
			print(eid)
			assert eid in sg.edge_df.index.tolist()

	# test to make sure adding a gtf without gene names
	# should use the gene_id instead
	def test_no_gene_name_gtf(self):
		sg = swan.SwanGraph()
		sg.add_dataset('test', 'input_files/Canx.gtf')

		gnames = sg.t_df.gname.tolist()
		gids = sg.t_df.gid.tolist()

		assert gnames == gids
	
	def test_no_gene_name_db(self):
		sg = swan.SwanGraph()
		sg.add_dataset('test', 'input_files/chr11_and_Tcf3_no_gname.db')

		gnames = sg.t_df.gname.tolist()
		gids = sg.t_df.gid.tolist()

		assert gnames == gids

def process_gtf():	
	sg = swan.SwanGraph()
	sg.add_dataset('test', 'input_files/weird_gtf_entries.gtf')

	return sg 

def check_pairs(control, test):
	print('control')
	print(control)
	print('test')
	print(test)
	for t in test:
		assert t in control
	assert len(test) == len(control)