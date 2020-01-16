import pytest
import sys
import numpy as np
sys.path.append('../utils/')
sys.path.append('../../refactor_splice_graph/')
from SpliceGraph import SpliceGraph
from utils import *


# i suppose this is also testing the df init because 
# everything needs to be here? And in the right order -\(' ')/-
class TestMergeSGs(object):
	def test_merge_sgs(self):
		a_gtf = 'input_files/annot.gtf'
		b_gtf = 'input_files/annot_2.gtf'
		sg = SpliceGraph()
		sg.add_dataset('a', gtf=a_gtf)
		sg.add_dataset('b', gtf=b_gtf)
		sg.order_transcripts()

		print(sg.loc_df.head())
		print(sg.edge_df.head())
		print(sg.t_df.head())

		# print(sg.loc_df[['chrom', 'coord', 'strand', 'a', 'b']])
		# print(sg.edge_df[['edge_type', 'a', 'b']])
		# print(sg.t_df[['path', 'a', 'b']])

		# check that the format of dfs are ok
		assert sg.loc_df.index.names == ['vertex_id']
		control = ['coord', 'chrom', 'strand', 'a', 'b', 'vertex_id',
				   'internal', 'TSS', 'TES', 'alt_TSS', 'alt_TES'] 
		test = sg.loc_df.columns.tolist()
		check_pairs(control, test)

		assert sg.edge_df.index.names == ['edge_id']
		control = ['v1', 'v2', 'edge_type', 'strand', 'a', 'b', 'edge_id'] 
		test = sg.edge_df.columns.tolist()
		check_pairs(control, test)

		assert sg.t_df.index.names == ['tid']
		control = ['tid', 'gid', 'gname', 'path', 'a', 'b']
		test = sg.t_df.columns.tolist()
		check_pairs(control, test)

		# test that loc_df merging happened correctly
		# query chr, coord, strand, a and b columns
		chrs = sg.loc_df['chrom'].tolist()
		control = [1,1,1,1,1,1,1,7,7,7,7,1,1,1,1,1,4,4]
		print('test chrs: ')
		print(chrs)
		print('control chrs: ')
		print(chrs)
		assert control == control

		coords = sg.loc_df['coord'].tolist()
		control = [1,90,100,500,600,900,1000,1,10,15,
				   20,2000,1500,1000,900,800,4000,1000]
		print('test coords: ')
		print(coords)
		print('control coords: ')
		print(control)
		assert coords == control

		strand = sg.loc_df['strand'].tolist()
		control = ['+','+','+','+','+','+','+','+','+',
				   '+','+','-','-','-','-','-','-','-']
		print('test strands: ')
		print(strand)
		print('control strands: ')
		print(control)
		assert strand == control

		a = sg.loc_df['a'].tolist()
		control = [bool(i) for i in [1,0,1,1,1,1,1,0,0,0,0,1,1,1,1,0,1,1]]
		print('test a presence: ')
		print(a)
		print('control a presence: ')
		print(control)
		assert a == control

		b = sg.loc_df['b'].tolist()
		control = [bool(i) for i in [1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0]]
		print('test b presence: ')
		print(b)
		print('control b presence: ')
		print(control)
		assert b == control

		# test that edge_df merging and id mapping happened correctly
		# query edge_id, edge_type, a and b columns
		edge_id = sg.edge_df['edge_id'].tolist()
		control = [(0,1),(0,2),(1,3),(2,3),(3,4),(4,5),(5,6),
				   (7,8),(8,9),(9,10),(11,12),(12,13),(13,14),(13,15),(16,17)]
		print('test edge_ids: ')
		print(edge_id)
		print('control edge_ids')
		print(control)
		assert edge_id == control

		edge_type = sg.edge_df['edge_type'].tolist()
		control = ['exon','exon','intron','intron','exon','intron','exon',
				   'exon','intron','exon','exon','intron','exon','exon','exon']
		print('test edge_types: ')
		print(edge_type)
		print('control edge_types: ')
		print(control)
		assert edge_type == control

		a = sg.edge_df['a'].tolist()
		control = [bool(i) for i in [0,1,0,1,1,1,1,0,0,0,1,1,1,0,1]]
		print('test a presence: ')
		print(a)
		print('control a presence: ')
		print(control)
		assert a == control

		b = sg.edge_df['b'].tolist()
		control = [bool(i) for i in [1,1,1,1,1,1,1,1,1,1,1,1,0,1,0]]
		print('test b presence: ')
		print(b)
		print('control b presence: ')
		print(control)
		assert b == control

		# test that t_df merging and id mapping happened correctly
		# query tid, path, a and b columns
		tid = sg.t_df['tid'].tolist()
		control = ['ENST01','ENST02','ENST03','ENST04','ENST07','ENST08']
		print('test tids: ')
		print(tid)
		print('control tids: ')
		print(control)
		assert tid == control

		paths = [tuple(path) for path in sg.t_df['path'].tolist()]
		control = [(0,2,3,4,5,6),(0,1,3,4,5,6),(11,12,13,14),(11,12,13,15),
				   (16,17),(7,8,9,10)]
		print('test paths: ')
		print(paths)
		print('control paths: ')
		print(control)
		assert paths == control

		a = sg.t_df['a'].tolist()
		control = [bool(i) for i in [1,0,1,0,1,0]]
		print('test a presence: ')
		print(a)
		print('control a presence: ')
		print(control)
		assert a == control

		b = sg.t_df['b'].tolist()
		control = [bool(i) for i in [1,1,0,1,0,1]]
		print('test b presence: ')
		print(b)
		print('control b presence: ')
		print(control)
		assert b == control

def check_pairs(control, test):
	print()
	print('control')
	print(control)
	print('test')
	print(test)
	for t in test:
		assert t in control
	assert len(test) == len(control)






