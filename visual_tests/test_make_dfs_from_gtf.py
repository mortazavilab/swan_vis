import sys
import os
from collections import defaultdict
lib_path = '/'.join(os.path.abspath(__file__).split('/')[0:-2])
sys.path.append(lib_path+'/utils/')
sys.path.append(lib_path)
from SpliceGraph import SpliceGraph
from utils import *
import cProfile

# gtf = '/Users/fairliereese/mortazavi_lab/ref/gencode.v24/gencode.v24.annotation.txt'
# gtf = 'input_files/gencode_200k.gtf'
# gtf = 'input_files/hippocampus_mapt.gtf'
gtf = 'input_files/test_combine_2.gtf'
# gtf = 'input_files/annot.gtf'
sg = SpliceGraph()

# print('dana')
# cProfile.run('sg.create_dfs_gtf_dana(gtf)')
# print('fairlie')
# cProfile.run('sg.create_dfs_gtf(gtf)')

# print('dana')
# sg.create_dfs_gtf_dana(gtf)
sg.add_annotation(gtf)

print(sg.loc_df.head())
print(sg.edge_df.head())
print(sg.t_df.head())

loc_dict = sg.loc_df.to_dict('index')
print(loc_dict)
sg.loc_df = pd.DataFrame.from_dict(loc_dict, orient='index')
# loc_df.rename({'index': 'vertex_id'}, inplace=True)
sg.loc_df.reset_index(drop=True, inplace=True)
# sg.loc_df = create_dupe_index(sg.loc_df, 'vertex_id')
# sg.loc_df = set_dupe_index(sg.loc_df, 'vertex_id')
sg.loc_df.index.names = ['vertex_id']
print(sg.loc_df)



# print(sg.loc_df)
# print(sg.edge_df)
# print(sg.t_df)

# print(transcripts)
# print(exons)