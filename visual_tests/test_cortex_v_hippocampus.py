import sys
import os
from collections import defaultdict
lib_path = '/'.join(os.path.abspath(__file__).split('/')[0:-2])
sys.path.append(lib_path+'/utils/')
sys.path.append(lib_path)
from SpliceGraph import SpliceGraph
# from utils import *
# from plotting_tools import * 
# from report_tools import *


# annot = sg.SpliceGraph(gtf='input_files/annot_mapt.gtf')
# cortex = sg.SpliceGraph(gtf='input_files/cortex_mapt.gtf')
# hippo = sg.SpliceGraph(gtf='input_files/hippocampus_mapt.gtf')

# merged = sg.merge_graphs(annot, cortex, 'annotation', 'Cortex')
# merged = sg.add_graph(merged, hippo, 'Hippo')

# G = merged.G
# loc_df = merged.loc_df
# edge_df = merged.edge_df
# t_df = merged.t_df



# # add cortex abundances to graph 
# file = 'input_files/mouse_brain_talon_abundance_filtered.tsv'
# merged.add_abundance_dataset(file, ['PB82', 'PB84'], 'Cortex')

# # add hippocampus abundances to graph
# file = 'input_files/mouse_brain_talon_abundance_filtered.tsv'
# merged.add_abundance_dataset(file, ['PB83', 'PB85'], 'Hippo')

# # report ???
# args = defaultdict()
# args['color_edges'] = True
# args['color_nodes'] = True
# args['color_alt_nodes'] = True
# args['indicate_dataset'] = False
# args['combine'] = False
# args['indicate_novel'] = True

# gen_report(merged, args, 'figures/mapt', browser=True, order='tss')
# gen_report(merged, args, 'figures/mapt', order='tss')

ab_file = 'input_files/mouse_brain_talon_abundance_filtered.tsv'

sg = SpliceGraph()
sg.add_annotation(gtf='input_files/annot_mapt.gtf')

# sg.add_dataset('cortex',
# 				gtf='input_files/cortex_mapt.gtf',
# 			    counts_file=ab_file,
# 			    count_cols=['PB82', 'PB84'])
# sg.add_dataset('hippocampus',
# 			   gtf='input_files/hippocampus_mapt.gtf',
# 			   counts_file=ab_file,
# 			   count_cols=['PB83', 'PB85'])

sg.plot_graph()
sg.save_fig('figures/test_mapt.png')

# print(sg.loc_df.head())
# print(sg.edge_df.head())
# print(sg.t_df.head())

# print()

# print(sg.pg.loc_df.head())
# print(sg.pg.edge_df.head())
# print(sg.pg.t_df.head())



