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

gtf1 = 'input_files/annot.gtf'
gtf2 = 'input_files/annot_2.gtf'

sg = SpliceGraph()
sg.add_dataset('a', gtf=gtf1)
sg.add_dataset('b', gtf=gtf2)
