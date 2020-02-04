import sys
import os
from collections import defaultdict
import swan as sw
import cProfile
import time
import numpy as np


# # creating the sg
# annot_gtf = 'input_files/gencode.vM21.annotation.gtf'
# hippocampus_gtf = 'input_files/hippocampus_talon.gtf'
# cortex_gtf = 'input_files/cortex_talon.gtf'
# ab_file = 'input_files/mouse_brain_talon_abundance_filtered.tsv'

# sg = sw.SwanGraph()
# sg.add_annotation(annot_gtf)
# sg.add_dataset('hippocampus', gtf=hippocampus_gtf,
# 			   counts_file=ab_file,
# 			   count_cols=['PB83', 'PB85'])
# sg.add_dataset('cortex', gtf=cortex_gtf,
# 			    counts_file=ab_file,
# 			    count_cols=['PB82', 'PB84'])
# sg.save_graph('input_files/cortex_hippocampus')

# loading the sg
sg = sw.SwanGraph()
sg.load_graph('input_files/cortex_hippocampus.p')


# finding "interesting" genes
# genes = sg.find_interesting_genes(how='isoform_switching')
# print(genes)

# plot the summary graph 
sg.plot_graph('ENSMUSG00000018411.17', indicate_dataset='cortex')



