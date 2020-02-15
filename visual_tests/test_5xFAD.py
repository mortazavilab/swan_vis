import sys
import os
from collections import defaultdict
import swan_vis as sw
import numpy as np

# # add each dataset to the SpliceGraph
# prefix = '/Users/fairliereese/mortazavi_lab/bin/swan_vis/visual_tests/'
# annot_gtf = '/Users/fairliereese/mortazavi_lab/ref/gencode.vM21/gencode.vM21.annotation.gtf'
# wt_1_gtf = prefix+'input_files/wt_1_filtered_talon.gtf'
# wt_2_gtf = prefix+'input_files/wt_2_filtered_talon.gtf'
# alz_1_gtf = prefix+'input_files/5xFAD_1_filtered_talon.gtf'
# alz_2_gtf = prefix+'input_files/5xFAD_2_filtered_talon.gtf'
# ab_file = prefix+'input_files/wt_5xFAD_filtered_talon_abundance.tsv'
# sg = sw.SwanGraph()
# sg.add_annotation(annot_gtf)
# sg.add_dataset('wt_1', wt_1_gtf,
#   counts_file=ab_file,
#   count_cols='PB132')
# sg.add_dataset('wt_2', wt_2_gtf,
#   counts_file=ab_file,
#   count_cols='PB133')
# sg.add_dataset('5xFAD_1', alz_1_gtf,
#   counts_file=ab_file,
#   count_cols='PB130')
# sg.add_dataset('5xFAD_2', alz_2_gtf,
#   counts_file=ab_file,
#   count_cols='PB131')

# # save the thing so we can just load it at any point
# sg.save_graph('input_files/wt_5xFAD_sg')

# load it back in
sg = sw.SwanGraph()
sg.load_graph('input_files/wt_5xFAD_sg.p')

# # plot a summary graph
# sg.plot_graph('ENSMUSG00000018411.17')
# sg.save_fig('../figures/ENSMUSG00000018411.17_summary.png')

# # plot a transcript path through a gene
# sg.plot_transcript_path('ENSMUST00000106992.9')
# sg.save_fig('../figures/ENSMUST00000106992.9.png')

# # plot each transcript through the gene
# sg.plot_each_transcript('ENSMUSG00000018411.17', '../figures/wt_5xFAD')

# # plot combined summary graph
# sg.plot_graph('ENSMUSG00000018411.17', combine=True)
# sg.save_fig('../figures/ENSMUSG00000018411.17_combined.png')

# # plot with indicate_novel
# sg.plot_graph('ENSMUSG00000018411.17', indicate_novel=True)
# sg.save_fig('../figures/ENSMUSG00000018411.17_novel.png')

# # plot with indicate_dataset
# sg.plot_graph('ENSMUSG00000018411.17', indicate_dataset='wt_1')
# sg.save_fig('../figures/ENSMUSG00000018411.17_wt_1.png')

# # plot transcript path browser
# sg.plot_transcript_path('ENSMUST00000106992.9', browser=True)
# sg.save_fig('../figures/ENSMUST00000106992.9_browser.png')

# generate report
# sg.gen_report('ENSMUSG00000018411.17', '../figures/wt_5xFAD')

# you are here
# # generate report with combine
# sg.gen_report('ENSMUSG00000018411.17', '../figures/wt_5xFAD', combine=True)

# # generate report with indicate_novel
# sg.gen_report('ENSMUSG00000018411.17', '../figures/wt_5xFAD', indicate_novel=True)

# # generate report with indicate_dataset
# sg.gen_report('ENSMUSG00000018411.17', '../figures/wt_5xFAD', indicate_dataset='wt_1')

# # generate report in browser style
# sg.gen_report('ENSMUSG00000018411.17', '../figures/wt_5xFAD', browser=True)

# # generate report with specific datasets
# sg.gen_report('ENSMUSG00000018411.17', 'figures/wt_5xFAD', datasets=['wt_1', '5xFAD_1'])

# # generate report with tpm shown
# sg.gen_report('ENSMUSG00000018411.17', 'figures/wt_5xFAD', tpm=True)

# # generate report with heatmap
# sg.gen_report('ENSMUSG00000018411.17', 'figures/wt_5xFAD', heatmap=True)

# # generate report ordered by tid
# sg.gen_report('ENSMUSG00000018411.17', 'figures/wt_5xFAD', order='tid', heatmap=True)

# # generate report ordered by TSS location
# sg.gen_report('ENSMUSG00000018411.17', 'figures/wt_5xFAD', order='tss', heatmap=True)

# # generate report ordered by TES location
# sg.gen_report('ENSMUSG00000018411.17', 'figures/wt_5xFAD', order='tes', heatmap=True)

# # include unexpressed
# sg.gen_report('ENSMUSG00000018411.17', 'figures/wt_5xFAD', order='tid', tpm=True, include_unexpressed=True)

# make a report of genes heavily involved in isoform switching
genes, t_df = sg.find_differentially_expressed_transcripts(['wt_1', 'wt_2'], ['5xFAD_1', '5xFAD_2'])
t_df.to_csv('input_files/wt_5xFAD_isoform_switches_t_df.csv')
# g_df.to_csv('input_files/wt_5xFAD_isoform_switches_g_df.csv')
print(genes)
sg.gen_report(genes, 'figures/wt_5xFAD_isoform_switches', heatmap=True)



