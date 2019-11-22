import sys
import os
lib_path = '/'.join(os.path.abspath(__file__).split('/')[0:-2])
sys.path.append(lib_path+'/utils/')
sys.path.append(lib_path)
import SpliceGraph as sg
import PlottedGraph as pg
from utils import *
from plotting_tools import * 

splice_graph = sg.SpliceGraph(gtf='input_files/test_combine_2.gtf')

args = defaultdict()
args['color_edges'] = True
args['color_nodes'] = True
args['color_alt_nodes'] = True
args['combine'] = False

plotted_graph = pg.PlottedGraph(splice_graph, args)
plot_graph(plotted_graph, args)
save_fig('figures/test_combine_2_exons_introns_alt_TSS_alt_TES.png')

args['combine'] = True
plotted_graph = pg.PlottedGraph(splice_graph, args)
plot_graph(plotted_graph, args)
save_fig('figures/test_combine_2_exons_introns_alt_TSS_alt_TES_combined.png')


