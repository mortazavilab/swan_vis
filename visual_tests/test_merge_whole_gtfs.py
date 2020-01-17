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


# # creating the sg
# annot_gtf = 'input_files/gencode.vM21.annotation.gtf'
# hippocampus_gtf = 'input_files/hippocampus_talon.gtf'
# cortex_gtf = 'input_files/cortex_talon.gtf'
# ab_file = 'input_files/mouse_brain_talon_abundance_filtered.tsv'

# sg = SpliceGraph()
# sg.add_annotation(annot_gtf)
# sg.add_dataset('hippocampus', gtf=hippocampus_gtf,
# 			   counts_file=ab_file,
# 			   count_cols=['PB83', 'PB85'])
# sg.add_dataset('cortex', gtf=cortex_gtf,
# 			    counts_file=ab_file,
# 			    count_cols=['PB82', 'PB84'])
# sg.save_graph('input_files/cortex_hippocampus')

# loading the sg
sg = SpliceGraph()
sg.load_graph('input_files/cortex_hippocampus.p')
genes = sg.find_interesting_genes(how='num_novel_isoforms')
print(genes)


