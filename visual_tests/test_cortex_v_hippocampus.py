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


annot = sg.SpliceGraph(gtf='input_files/annot_mapt.gtf')
cortex = sg.SpliceGraph(gtf='input_files/cortex_mapt.gtf')
hippo = sg.SpliceGraph(gtf='input_files/hipp_mapt.gtf')

merged = sg.add_graph(annot, cortex, 'cortex')
merged = sg.add_graph(merged, hippo, 'hippo')

loc_df = merged.loc_df
edge_df = merged.edge_df
t_df = merged.t_df

print(merged.loc_df.head())
print(merged.edge_df.head())
print(merged.t_df.head())

print(merged.loc_df.tail())
print(merged.edge_df.tail())
print(merged.t_df.tail())

# this makes sense because we filtered based on reproducibility across replicates so everything
# in one dataset must be present in the other! (wait, is this right?)

print('loc_df')
print(len(loc_df.loc[loc_df.dataset_cortex == False].index))
print(len(loc_df.loc[loc_df.dataset_hippo == False].index))
# print(len(loc_df.loc[loc_df.dataset_annot == False].index))
# print(loc_df.loc[loc_df.dataset_annot == False])
print(loc_df.loc[loc_df.dataset_cortex == False])
print(loc_df.loc[loc_df.dataset_hippo == False])

print('edge_df')
print(len(edge_df.loc[edge_df.dataset_cortex == False].index))
print(len(edge_df.loc[edge_df.dataset_hippo == False].index))

print('t_df')
print(len(t_df.loc[t_df.dataset_cortex == False].index))
print(len(t_df.loc[t_df.dataset_hippo == False].index))

# add abundances to graph 
file = 'input_files/mouse_brain_talon_abundance_filtered.tsv'

# cortex
merged.add_abundance_dataset(file, ['PB82', 'PB84'], 'cortex')

