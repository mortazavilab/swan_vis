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

splice_graph = sg.SpliceGraph(gtf='input_files/annot_mapt.gtf')

args = defaultdict()
args['color_edges'] = True
args['color_nodes'] = True
args['color_alt_nodes'] = True
args['indicate_dataset'] = False
args['combine'] = False

plotted_graph = pg.PlottedGraph(splice_graph, args)
plot_graph(plotted_graph, args)
save_fig('mapt_exons_introns_alt_TSS_alt_TES.png')

args['combine'] = True
plotted_graph = pg.PlottedGraph(splice_graph, args)
plot_graph(plotted_graph, args)
save_fig('mapt_exons_introns_alt_TSS_alt_TES_combined.png')

# each transcript through the combined graph
for tid in plotted_graph.t_df.tid.tolist():
	print(tid)
	oname = 'figures/mapt_combined_{}.png'.format(tid)
	path = plotted_graph.t_df.loc[tid, 'path']
	plot_overlaid_path(plotted_graph, path, args)
	save_fig(oname)

# also plot genome browser style!
for tid in plotted_graph.t_df.tid.tolist():
	print(tid)
	oname = 'figures/mapt_browser_{}.png'.format(tid)
	plot_path_browser(splice_graph, tid, oname)