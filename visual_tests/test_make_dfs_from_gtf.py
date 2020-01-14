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

# gtf = '/Users/fairliereese/mortazavi_lab/ref/gencode.v24/gencode.v24.annotation.txt'
gtf = 'input_files/gencode_200k.gtf'
# gtf = 'input_files/hippocampus_mapt.gtf'
# gtf = 'input_files/test_combine_2.gtf'
# gtf = 'input_files/annot.gtf'
sg = SpliceGraph()

# print('dana')
# cProfile.run('sg.create_dfs_gtf_dana(gtf)')
# print('fairlie')
# cProfile.run('sg.create_dfs_gtf(gtf)')

# print('dana')
sg.create_dfs_gtf(gtf)
# sg.add_annotation(gtf)

# def label_node_types(locs, vertex_ids, node_type):
# 	for vertex_id in vertex_ids:
# 		locs[vertex_id][node_type] = True
# 	return locs 

# start_time = time.time()

# loc_df = sg.loc_df
# t_df = sg.t_df

# locs = sg.loc_df.loc[:, ['vertex_id', 'coord']].copy(deep=True)
# locs['internal'] = False
# locs['TSS'] = False
# locs['TES'] = False
# locs['alt_TSS'] = False
# locs['alt_TES'] = False
# locs.drop(['coord', 'vertex_id'], axis=1, inplace=True)
# locs = locs.to_dict('index')

# t_df.reset_index(drop=True, inplace=True)
# t_df.set_index('gid', inplace=True)

# # label each TSS and TES
# paths = t_df.path.tolist()
# tss = np.unique([path[0] for path in paths])
# locs = label_node_types(locs, tss, 'TSS')
# # loc_df[tss]['TSS'] = True
# tes = np.unique([path[-1] for path in paths])
# # loc_df.loc[tes]['TES'] = True
# locs = label_node_types(locs, tes, 'TES')
# internal = np.unique([n for path in paths for n in path[1:-1]])
# # loc_df.loc[internal]['internal'] = True
# locs = label_node_types(locs, internal, 'internal')

# for gid in t_df.index.unique().tolist():

# 	gene_paths = t_df.loc[gid, 'path']
# 	if type(gene_paths) != list:
# 		gene_paths = gene_paths.tolist()
# 	else: 
# 		gene_paths = [gene_paths]

# 	if len(gene_paths) > 1:
# 		tss = [path[0] for path in gene_paths] # could combine this and next line
# 		tes = [path[0] for path in gene_paths] # and return a tuple
# 		if len(set(tss)) > 1:  # probably don't need to use set, each entry should be unique already
# 			# loc_df.loc[tss, 'alt_TSS'] = True
# 			locs = label_node_types(locs, tss, 'alt_TSS')
# 		if len(set(tes)) > 1:
# 			locs = label_node_types(locs, tes, 'alt_TES')
# 			# loc_df.loc[tes, 'alt_TES'] = True

# # create df from locs dict 
# locs = pd.DataFrame.from_dict(locs, orient='index')

# print(locs.head())
# print(loc_df.head())

# # append old loc_df with node types loc_df
# loc_df = pd.concat([loc_df, locs], axis=1)

print('time to get loc types')
print(time.time()-start_time)

print(loc_df.head())


# print(sg.loc_df)
# print(sg.edge_df)
# print(sg.t_df)

# print(transcripts)
# print(exons)