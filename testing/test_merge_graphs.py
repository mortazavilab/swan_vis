import sys
import os
lib_path = '/'.join(os.path.abspath(__file__).split('/')[0:-2])
sys.path.append(lib_path+'/utils/')
sys.path.append(lib_path)
import SpliceGraph as sg
import PlottedGraph as pg
from utils import *
from plotting_tools import * 


a = sg.SpliceGraph(gtf='input_files/annot_mapt.gtf')
b = sg.SpliceGraph(gtf='input_files/obs_mapt.gtf')

merged = sg.merge_graphs(a, b, 'annot', 'obs')

args['color_edges'] = True
args['color_nodes'] = True
args['color_alt_nodes'] = True
args['combine'] = False
args['indicate_dataset'] = 'obs'

pg.plot_graph()