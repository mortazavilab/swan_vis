import sys
import os
from collections import defaultdict
lib_path = '/'.join(os.path.abspath(__file__).split('/')[0:-2])
sys.path.append(lib_path+'/utils/')
sys.path.append(lib_path)
from SpliceGraph import SpliceGraph
from utils import *
import cProfile
import time
import numpy as np

gtf = '/Users/fairliereese/mortazavi_lab/ref/gencode.v24/gencode.v24.annotation.txt'
# gtf = 'input_files/gencode_200k.gtf'
# gtf = 'input_files/hippocampus_mapt.gtf'
# gtf = 'input_files/test_combine_2.gtf'
# gtf = 'input_files/annot.gtf'
sg = SpliceGraph()

# print('dana')
# cProfile.run('sg.create_dfs_gtf_dana(gtf)')
# print('fairlie')
# cProfile.run('sg.create_dfs_gtf(gtf)')

# print('dana')
# sg.create_dfs_gtf(gtf)
sg.add_annotation(gtf)

# start_time = time.time()

# loc_df = sg.loc_df
# t_df = sg.t_df

# def label_node_types(locs, vertex_ids, node_type):
# 	for vertex_id in vertex_ids:
# 		locs[vertex_id][node_type] = True
# 	return locs 

# # create a dictionary to hold loc info to speed stuff up
# locs = loc_df.loc[:, ['vertex_id', 'coord']].copy(deep=True)
# locs['internal'] = False
# locs['TSS'] = False
# locs['TES'] = False
# locs['alt_TSS'] = False
# locs['alt_TES'] = False
# locs.drop(['coord', 'vertex_id'], axis=1, inplace=True)
# locs = locs.to_dict('index')

# # label each TSS and TES
# paths = t_df.path.tolist()
# tss = np.unique([path[0] for path in paths])
# locs = label_node_types(locs, tss, 'TSS')
# tes = np.unique([path[-1] for path in paths])
# locs = label_node_types(locs, tes, 'TES')
# internal = np.unique([n for path in paths for n in path[1:-1]])
# locs = label_node_types(locs, internal, 'internal')

# # also create a dictionary of gid: [path1, ... pathn] to speed up
# path_dict = defaultdict(list)
# for tid, entry in t_df.iterrows():
# 	path_dict[entry.gid].append(entry.path)

# # label alt TES/TSS
# for gid in path_dict.keys():
# 	paths = path_dict[gid]
# 	if len(paths) > 1:
# 		tss = np.unique([path[0] for path in paths])
# 		tes = np.unique([path[-1] for path in paths])
# 		if len(tss) > 1:  
# 			locs = label_node_types(locs, tss, 'alt_TSS')
# 		if len(tes) > 1:
# 			locs = label_node_types(locs, tes, 'alt_TES')

# print('time to get loc types')
# print(time.time()-start_time)

