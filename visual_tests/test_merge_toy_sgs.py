import sys
import os
from collections import defaultdict
import swan as sw
import cProfile
import time
import numpy as np

gtf1 = 'input_files/annot.gtf'
gtf2 = 'input_files/annot_2.gtf'

sg = sw.SwanGraph()
sg.add_dataset('a', gtf=gtf1)
sg.add_dataset('b', gtf=gtf2)
