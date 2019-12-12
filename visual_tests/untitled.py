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
cortex = sg.SpliceGraph(gtf='input_files/cortex_vgf.gtf')
