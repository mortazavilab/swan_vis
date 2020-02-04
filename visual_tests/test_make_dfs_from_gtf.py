import sys
import os
from collections import defaultdict
import swan as sw
import cProfile
import time
import numpy as np

gtf = '/Users/fairliereese/mortazavi_lab/ref/gencode.v24/gencode.v24.annotation.txt'
# gtf = 'input_files/gencode_200k.gtf'
# gtf = 'input_files/hippocampus_mapt.gtf'
# gtf = 'input_files/test_combine_2.gtf'
# gtf = 'input_files/annot.gtf'
sg = sw.SwanGraph()

# print('dana')
# cProfile.run('sg.create_dfs_gtf_dana(gtf)')
# print('fairlie')
# cProfile.run('sg.create_dfs_gtf(gtf)')

# print('dana')
# sg.create_dfs_gtf(gtf)
sg.add_annotation(gtf)

# start_time = time.time()
# print(sg.loc_df['coord'])
# print(sg.loc_df.loc[:, ['coord', 'strand', 'vertex_id']])
# print(sg.edge_df.loc[:, ['v1', 'v2']])
# print(sg.t_df['path'])
# sg.update_ids()
# print('after updating ids')
# print('hello world')
# print(sg.loc_df.loc[:, ['coord', 'strand', 'vertex_id']])
# print(sg.edge_df.loc[:, ['v1', 'v2']])
# print(sg.t_df['path'])
# print('time to get ordered id map')
# print(time.time()-start_time)


