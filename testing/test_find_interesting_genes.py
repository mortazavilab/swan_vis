import pytest
import sys
import numpy as np
import swan_vis as swan
import pandas as pd

class TestFindInterestingGenes(object):

	# 
	def test_num_novel_known_isoforms(self):
		sg = swan.SwanGraph()
		sg.t_df = pd.DataFrame({'tid': [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],
							    'gid': [0,0,0,0,1,1,1,2,2,2,2,2,2,2,3,3],
							    'annotation': [True,True,True,True,False,True,True,False,False,False,True,False,False,True,True,True],
							    'a': [False,True,True,False,True,False,False,True,True,False,False,True,False,False,False,False],
							    'b': [False,False,True,True,True,False,True,False,True,True,True,False,True,False,False,False]})
		sg.t_df = swan.create_dupe_index(sg.t_df, 'tid')
		sg.t_df = swan.set_dupe_index(sg.t_df, 'tid')
		sg.datasets = ['annotation', 'a', 'b']

		genes, g_df = sg.find_genes_with_novel_isoforms()

		# check that the genes were returned in the right order
		print(genes)
		control = [2,1,0]
		assert genes == control

		# check that gene 3 didn't get included
		g_df_genes = g_df.index.tolist()
		print(g_df_genes)
		control = [0,1,2]
		check_pairs(control, g_df_genes)

		# check that sum of known/novel models is correct
		known_models = g_df.loc[0, 'known']
		print('gene 0 num known:')
		print(known_models)
		known_control = 3
		novel_models = g_df.loc[0, 'novel']
		print('gene 0 num novel:')
		print(novel_models)
		novel_control = 0
		assert known_models == known_control
		assert novel_models == novel_control

		known_models = g_df.loc[1, 'known']
		print('gene 1 num known:')
		print(known_models)
		known_control = 1
		novel_models = g_df.loc[1, 'novel']
		print('gene 1 num novel:')
		print(novel_models)
		novel_control = 1
		assert known_models == known_control
		assert novel_models == novel_control

		known_models = g_df.loc[2, 'known']
		print('gene 2 num known:')
		print(known_models)
		known_control = 1
		novel_models = g_df.loc[2, 'novel']
		print('gene 2 num novel:')
		print(novel_models)
		novel_control = 5
		assert known_models == known_control
		assert novel_models == novel_control

		# make sure that when we're trying this w/o an annotation,
		# we raise the error
		sg.t_df.drop('annotation', axis=1, inplace=True)
		sg.datasets.remove('annotation')
		with pytest.raises(Exception) as excinfo:
			genes, g_df = sg.find_genes_with_novel_isoforms()
		assert 'No annotation data' in str(excinfo.value)


def check_pairs(control, test):
	print('control')
	print(control)
	for t in test:
		print(t)
		assert t in control
	assert len(test) == len(control)

