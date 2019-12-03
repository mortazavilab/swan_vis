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


a = sg.SpliceGraph(gtf='input_files/annot_mapt.gtf')
b = sg.SpliceGraph(gtf='input_files/obs_mapt.gtf')

merged = sg.merge_graphs(a, b, 'annotation', 'obs')

args = defaultdict()
args['color_edges'] = True
args['color_nodes'] = True
args['color_alt_nodes'] = True
args['combine'] = False
args['indicate_dataset'] = 'obs'
args['indicate_novel'] = False

# using indicate_dataset
plotted_graph = pg.PlottedGraph(merged, args)
plot_graph(plotted_graph, args)
save_fig('figures/obs_annot_mapt_dataset.png')

# using indicate_novel
args['indicate_datset'] = False
args['indicate_novel'] = True
plotted_graph = pg.PlottedGraph(merged, args)
plot_graph(plotted_graph, args)
save_fig('figures/obs_annot_mapt_novel.png')