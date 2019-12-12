import sys
import os
from collections import defaultdict
lib_path = '/'.join(os.path.abspath(__file__).split('/')[0:-2])
sys.path.append(lib_path+'/utils/')
sys.path.append(lib_path)
import SpliceGraph as sg
import PlottedGraph as pg
from utils import *
from plotting_tools import * 
from report_tools import *


# mapt
annot = sg.SpliceGraph(gtf='input_files/gaby/annot_mapt.gtf')
wt = sg.SpliceGraph(gtf='input_files/gaby/WT_mapt.gtf')
xFAD = sg.SpliceGraph(gtf='input_files/gaby/5xFAD_mapt.gtf')

merged = sg.merge_graphs(annot, wt, 'annotation', 'WT')
merged = sg.add_graph(merged, xFAD, '5xFAD')

G = merged.G
loc_df = merged.loc_df
edge_df = merged.edge_df
t_df = merged.t_df

# add cortex abundances to graph 
file = 'input_files/gaby/Cortex_talon_abundance_filtered.tsv'
merged.add_abundance_dataset(file, ['PB82', 'PB84'], 'WT')

# add hippocampus abundances to graph
file = 'input_files/gaby/Cortex_talon_abundance_filtered.tsv'
merged.add_abundance_dataset(file, ['PB130', 'PB131'], '5xFAD')

# report ???
args = defaultdict()
args['color_edges'] = True
args['color_nodes'] = True
args['color_alt_nodes'] = True
args['indicate_dataset'] = False
args['combine'] = False
args['indicate_novel'] = True

gen_report(merged, args, 'gaby/figures/mapt', order='expression')
# gen_report(merged, args, 'gaby/figures/mapt', order='expression')


# cst7
annot = sg.SpliceGraph(gtf='input_files/gaby/annot_cst7.gtf')
wt = sg.SpliceGraph(gtf='input_files/gaby/WT_cst7.gtf')
xFAD = sg.SpliceGraph(gtf='input_files/gaby/5xFAD_cst7.gtf')

merged = sg.merge_graphs(annot, wt, 'annotation', 'WT')
merged = sg.add_graph(merged, xFAD, '5xFAD')

G = merged.G
loc_df = merged.loc_df
edge_df = merged.edge_df
t_df = merged.t_df

# add cortex abundances to graph 
file = 'input_files/gaby/Cortex_talon_abundance_filtered.tsv'
merged.add_abundance_dataset(file, ['PB82', 'PB84'], 'WT')

# add hippocampus abundances to graph
file = 'input_files/gaby/Cortex_talon_abundance_filtered.tsv'
merged.add_abundance_dataset(file, ['PB130', 'PB131'], '5xFAD')

# report ???
args = defaultdict()
args['color_edges'] = True
args['color_nodes'] = True
args['color_alt_nodes'] = True
args['indicate_dataset'] = False
args['combine'] = False
args['indicate_novel'] = True

gen_report(merged, args, 'gaby/figures/cst7', order='expression')
# gen_report(merged, args, 'gaby/figures/cst7', order='expression')

# trem2
annot = sg.SpliceGraph(gtf='input_files/gaby/annot_trem2.gtf')
wt = sg.SpliceGraph(gtf='input_files/gaby/WT_trem2.gtf')
xFAD = sg.SpliceGraph(gtf='input_files/gaby/5xFAD_trem2.gtf')

merged = sg.merge_graphs(annot, wt, 'annotation', 'WT')
merged = sg.add_graph(merged, xFAD, '5xFAD')

G = merged.G
loc_df = merged.loc_df
edge_df = merged.edge_df
t_df = merged.t_df

# add cortex abundances to graph 
file = 'input_files/gaby/Cortex_talon_abundance_filtered.tsv'
merged.add_abundance_dataset(file, ['PB82', 'PB84'], 'WT')

# add hippocampus abundances to graph
file = 'input_files/gaby/Cortex_talon_abundance_filtered.tsv'
merged.add_abundance_dataset(file, ['PB130', 'PB131'], '5xFAD')

# report ???
args = defaultdict()
args['color_edges'] = True
args['color_nodes'] = True
args['color_alt_nodes'] = True
args['indicate_dataset'] = False
args['combine'] = False
args['indicate_novel'] = True

gen_report(merged, args, 'gaby/figures/trem2', order='expression')
# gen_report(merged, args, 'gaby/figures/trem2', order='expression')




