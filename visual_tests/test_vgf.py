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
annot = sg.SpliceGraph(gtf='input_files/annot_vgf.gtf')
novel = sg.SpliceGraph(gtf='input_files/cortex_hippocampus_vgf.gtf')

args = defaultdict()
args['color_edges'] = True
args['color_nodes'] = True
args['color_alt_nodes'] = True
args['indicate_dataset'] = False
args['combine'] = False
args['indicate_novel'] = False

plotted_graph = pg.PlottedGraph(annot, args)
plot_graph(plotted_graph, args)
save_fig('figures/vgf_annotation.png')

merged = sg.merge_graphs(annot, novel, 'annotation', 'Novel')
args['indicate_novel'] = True

plotted_graph = pg.PlottedGraph(merged, args)
plot_graph(plotted_graph, args)
save_fig('figures/vgf_novel.png')
