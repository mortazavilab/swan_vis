import sys
import os
from collections import defaultdict
import swan as sw

gtf = 'input_files/annot.gtf'

sg = sw.SwanGraph()
sg.add_annotation(gtf)