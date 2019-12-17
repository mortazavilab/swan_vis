import sys
import os
lib_path = '/'.join(os.path.abspath(__file__).split('/')[0:-2])
sys.path.append(lib_path+'/utils/')
sys.path.append(lib_path)
from SpliceGraph import SpliceGraph


# splice_graph = sg.SpliceGraph(gtf='input_files/test_combine_2.gtf')

# args = defaultdict()
# args['color_edges'] = True
# args['color_nodes'] = True
# args['color_alt_nodes'] = True
# args['indicate_dataset'] = False
# args['combine'] = False

# plotted_graph = pg.PlottedGraph(splice_graph, args)
# plot_graph(plotted_graph, args)
# save_fig('figures/test_combine_2_exons_introns_alt_TSS_alt_TES.png')

# args['combine'] = True
# plotted_graph = pg.PlottedGraph(splice_graph, args)
# plot_graph(plotted_graph, args)
# save_fig('figures/test_combine_2_exons_introns_alt_TSS_alt_TES_combined.png')


sg = SpliceGraph()
sg.add_dataset('test',
	gtf='input_files/test_combine_2.gtf')

sg.plot_graph()
sg.save_fig('figures/test_combine_not_combined.png')