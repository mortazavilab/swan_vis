import sys
import os
from collections import defaultdict
lib_path = '/'.join(os.path.abspath(__file__).split('/')[0:-2])
sys.path.append(lib_path+'/utils/')
sys.path.append(lib_path)
import SpliceGraph as sw
import PlottedGraph as pg
from utils import *
from plotting_tools import * 

sg = sw.SpliceGraph()
sg.add_annotation(db='input_files/combine.db')

print(sg.loc_df.head())
print(sg.edge_df.head())
print(sg.t_df.head())