import sys
import os
from collections import defaultdict
lib_path = '/'.join(os.path.abspath(__file__).split('/')[0:-2])
sys.path.append(lib_path+'/utils/')
sys.path.append(lib_path)
from SpliceGraph import SpliceGraph
from utils import *
import time
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


# with abundances, takes longer to run
sg.add_dataset('cortex',
				gtf='input_files/cortex_mapt.gtf',
			    counts_file=ab_file,
			    count_cols=['PB82', 'PB84'])
sg.add_dataset('hippocampus',
			   gtf='input_files/hippocampus_mapt.gtf',
			   counts_file=ab_file,
			   count_cols=['PB83', 'PB85'])

# # without abundances, much faster to run
# sg.add_dataset('cortex',
# 				gtf='input_files/cortex_mapt.gtf')

# sg.add_dataset('hippocampus',
# 			   gtf='input_files/hippocampus_mapt.gtf')

# # testing
# sg.plot_graph(indicate_dataset='cortex')
# sg.save_fig('figures/test_mapt_cortex.png')
# sg.plot_graph(indicate_dataset='hippocampus')
# sg.save_fig('figures/test_mapt_hippocampus.png') 
# sg.plot_graph(indicate_novel=True)
# sg.save_fig('figures/test_mapt_novel.png')

# sg.plot_transcript_path('ENSMUST00000106992.9')
# sg.save_fig('figures/test_mapt_ENSMUST00000106992.9.png')
# sg.plot_transcript_path('ENSMUST00000106992.9', browser=True)
# sg.save_fig('figures/test_mapt_ENSMUST00000106992_browser.png')

# print(sg.t_df.head())

# # plot each transcript test
# sg.plot_each_transcript('ENSMUSG00000018411.17','figures/test_mapt')

# # plot each transcript browser test
# sg.plot_each_transcript('ENSMUSG00000018411.17', 'figures/test_mapt', browser=True)

# report test
# gid = 'ENSMUSG00000018411.17'
# counts_cols = sg.get_count_cols(['cortex', 'hippocampus'])
# report_tids = sg.t_df.loc[(sg.t_df.gid==gid)&(sg.t_df[counts_cols].sum() > 0),
# 			  	'tid'].tolist()
# sg.t_df.to_csv('test_t_df.csv')	
# print(report_tids)
# sg.gen_report('ENSMUSG00000018411.17', 'figures/test_mapt', order='expression', tpm=True)
# sg.gen_report('ENSMUSG00000018411.17', 'figures/test_mapt', order='expression')
# sg.gen_report('ENSMUSG00000018411.17', 'figures/test_mapt', order='expression', browser=True)




# print(sg.loc_df[['annotation', 'cortex', 'hippocampus']])
# print(sg.loc_df.head())
# print(sg.edge_df.head())
# print(sg.t_df.head())

# print()

# print(sg.pg.loc_df.head())
# print(sg.pg.edge_df.head())
# print(sg.pg.t_df.head())



