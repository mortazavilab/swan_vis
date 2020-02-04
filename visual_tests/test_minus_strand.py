import sys
import os
from collections import defaultdict
import swan as sw

sg = sw.SwanGraph()
sg.add_annotation(gtf='input_files/annotation_ENSMUSG00000051951.5.gtf')

# sg.add_dataset('cortex',
# 				gtf='input_files/cortex_ENSMUSG00000051951.5.gtf')

sg.add_dataset('hippocampus',
			   gtf='input_files/hippocampus_ENSMUSG00000051951.5.gtf')

sg.plot_transcript_path('ENSMUST00000159265.1')
sg.save_fig('figures/test_Xkr4_ENSMUST00000159265.1.png')
sg.plot_transcript_path('ENSMUST00000159265.1', browser=True)
sg.save_fig('figures/test_Xkr4_ENSMUST00000159265.1_browser.png')



