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


annot = sg.SpliceGraph(gtf='input_files/annot_mapt.gtf')
cortex = sg.SpliceGraph(gtf='input_files/cortex_mapt.gtf')
hippo = sg.SpliceGraph(gtf='input_files/hippocampus_mapt.gtf')

merged = sg.add_graph(annot, cortex, 'Cortex')
merged = sg.add_graph(merged, hippo, 'Hippo')

loc_df = merged.loc_df
edge_df = merged.edge_df
t_df = merged.t_df


# # testing
# print(merged.loc_df.head())
# print(merged.edge_df.head())
# print(merged.t_df.head())

# print(merged.loc_df.tail())
# print(merged.edge_df.tail())
# print(merged.t_df.tail())

# # this makes sense because we filtered based on reproducibility across replicates so everything
# # in one dataset must be present in the other! (wait, is this right?)

# print('loc_df')
# print(len(loc_df.loc[loc_df.dataset_cortex == False].index))
# print(len(loc_df.loc[loc_df.dataset_hippo == False].index))
# print(loc_df.loc[loc_df.dataset_cortex == False])
# print(loc_df.loc[loc_df.dataset_hippo == False])

# print('edge_df')
# print(len(edge_df.loc[edge_df.dataset_cortex == False].index))
# print(len(edge_df.loc[edge_df.dataset_hippo == False].index))
# print(edge_df.loc[edge_df.dataset_cortex == False])
# print(edge_df.loc[edge_df.dataset_hippo == False])

# print('t_df')
# print(len(t_df.loc[t_df.dataset_cortex == False].index))
# print(len(t_df.loc[t_df.dataset_hippo == False].index))
# print(t_df.loc[t_df.dataset_cortex == False])
# print(t_df.loc[t_df.dataset_hippo == False])


# add cortex abundances to graph 
file = 'input_files/mouse_brain_talon_abundance_filtered.tsv'
merged.add_abundance_dataset(file, ['PB82', 'PB84'], 'Cortex')

# add hippocampus abundances to graph
file = 'input_files/mouse_brain_talon_abundance_filtered.tsv'
merged.add_abundance_dataset(file, ['PB83', 'PB85'], 'Hippo')

# # testing
# print(merged.t_df.columns)
# print(merged.t_df[['tid','dataset_hippo_counts','dataset_cortex_counts']].head())

# report ???
args = defaultdict()
args['color_edges'] = True
args['color_nodes'] = True
args['color_alt_nodes'] = True
args['indicate_dataset'] = False
args['combine'] = False

# gen_report(merged, args, browser=True)
gen_report(merged, args, 'figures/mapt')




