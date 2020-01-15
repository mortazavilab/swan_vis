import sys
import os
from collections import defaultdict
lib_path = '/'.join(os.path.abspath(__file__).split('/')[0:-2])
sys.path.append(lib_path+'/utils/')
sys.path.append(lib_path)
from SpliceGraph import SpliceGraph
from utils import *
import cProfile
import time
import numpy as np

annot_gtf = 'input_files/gencode.vM21.annotation.gtf'
hippocampus_gtf = 'input_files/hippocampus_talon.gtf'

sg = SpliceGraph()
sg.add_annotation(annot_gtf)
sg.add_dataset('hippocampus', gtf=hippocampus_gtf)