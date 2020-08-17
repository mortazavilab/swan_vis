### Getting Started

import swan_vis as swan

annot_gtf = 'data/gencode.v29.annotation.gtf'
hep_1_gtf = 'data/hepg2_1_talon.gtf'
hep_2_gtf = 'data/hepg2_2_talon.gtf'
hff_1_gtf = 'data/hffc6_1_talon.gtf'
hff_2_gtf = 'data/hffc6_2_talon.gtf'
hff_3_gtf = 'data/hffc6_3_talon.gtf'
ab_file = 'data/all_talon_abundance_filtered.tsv'
talon_db = 'data/talon.db'

# initialize a new SwanGraph
sg = swan.SwanGraph() 

# add an annotation transcriptome 
sg.add_annotation(annot_gtf)

# add a dataset's transcriptome and abundance information to
# the SwanGraph
sg.add_dataset('HepG2_1', hep_1_gtf,
    counts_file=ab_file,
    count_cols='hepg2_1')
sg.add_dataset('HepG2_2', hep_2_gtf,
    counts_file=ab_file,
    count_cols='hepg2_2')
sg.add_dataset('HFFc6_1', hff_1_gtf,
    counts_file=ab_file,
    count_cols='hffc6_1')
sg.add_dataset('HFFc6_2', hff_2_gtf,
    counts_file=ab_file,
    count_cols='hffc6_2')
sg.add_dataset('HFFc6_3', hff_3_gtf,
    counts_file=ab_file,
    count_cols='hffc6_3')

# save the SwanGraph as a Python pickle file
sg.save_graph('swan')

### Analysis

dataset_groups = [['HepG2_1','HepG2_2'],
				  ['HFFc6_1','HFFc6_2','HFFc6_3']]

# perform a differential gene expression 
# Wald test on the provided two lists of datasets
sg.de_gene_test(dataset_groups)

# return a list of gene ids and their corresponding entries in the 
# differential expression test summary table for a given q value
gene_ids, gene_summary = sg.get_de_genes(q=0.05)

dataset_groups = [['HepG2_1','HepG2_2'],
				  ['HFFc6_1','HFFc6_2','HFFc6_3']]

# perform a differential transcript expression 
# Wald test on the provided two lists of datasets
sg.de_transcript_test(dataset_groups)

# return a list of transcript ids and their corresponding entries in the 
# differential expression test summary table for a given q value
transcript_ids, transcript_summary = sg.get_de_transcripts(q=0.05)

# find genes that exhibit isoform switching with a given input
# q value significance threshold
# returned entries will be those with de gene q-value > 0.05
# and de transcript q-value <= 0.05
is_genes, is_table = sg.find_isoform_switching_genes(q=0.05)

# returns a list of genes, transcripts, and specific edges in 
# the SwanGraph with novel exon skipping events
es_genes, es_transcripts, es_edges = sg.find_es_genes()

# returns a list of genes, transcripts, and specific edges in 
# the SwanGraph with novel intron retention events
ir_genes, ir_transcripts, ir_edges = sg.find_ir_genes()


### Visualization

# plot a gene summary graph
sg.plot_graph('ADRM1')

# save the currently-plotted figure with the given filename and location
swan.save_fig('figures/my_gene_summary.png')

# plot a gene summary graph with novel splice sites and 
# splice junctions plotted as outlined nodes and dashed edges respectively
sg.plot_graph('ADRM1', indicate_novel=True)

# plot a gene summary graph with splice sites and splice junctions
# that are present in the indicated dataset as outlined nodes
# and dashed edges respectively
sg.plot_graph('ADRM1', indicate_dataset='HepG2_1')

# plot a gene summary graph with novel splice sites and junctions
# outlined and dashed respectively
# save the figure with the given prefix
sg.plot_graph('ADRM1', indicate_novel=True, prefix='figures/adrm1')

# plot the path of a specific transcript through its parent gene
# summary graph for a given transcript
sg.plot_transcript_path('TALONT000301953')

# plot the path of a specific transcript through its parent gene summary
# graph for a given transcript
# plot novel splice sites as outlined nodes
# plot novel splice junctions as dashed edges
sg.plot_transcript_path('TALONT000301953', indicate_novel=True)

# plot the path of a specific transcript through its parent gene summary
# graph for a given transcript
# plot splice sites from the given dataset as outlined nodes
# plot splice junctions from the given dataset as dashed edges
sg.plot_transcript_path('TALONT000301953', indicate_dataset='HFFc6_1')

# plot the traditional browser-style representation 
# for a given transcript
sg.plot_transcript_path('TALONT000301953', browser=True)

# generate a report for the given gene 
# save it with the given filepath prefix
# display abundance information in the form of a heatmap
# include differential transcript expression test significance results
# (defaults to significance threshold q >= 0.05)
# display the novelty category associated with the transcript 
# display novel splice sites and junctions 
# as outlined nodes and dashed edges respectively
sg.gen_report('ADRM1',
              prefix='figures/adrm1_paper',
              heatmap=True,
              include_qvals=True,
              novelty=True, 
              indicate_novel=True)

# generate a report for the given gene 
# save it with the given filepath prefix
# combine datasets using the input dataset groupings and group names
# display abundance information in the form of a heatmap
# display novel splice sites and junctions 
# as outlined nodes and dashed edges respectively
sg.gen_report('ADRM1', prefix='figures/adrm1',
              dataset_groups=[['HepG2_1','HepG2_2'],['HFFc6_1','HFFc6_2']],
              dataset_group_names=['HepG2', 'HFFc6'],
              heatmap=True,
              indicate_novel=True)

# generate reports for the given genes
# with the given prefix
# only display the input datasets
# display abundance information as TPM numbers
# order the entries in the report by the genomic coordinate
# of each transcript's transcription start site (TSS)
sg.gen_report(['ADRM1','PSAP'], prefix='figures/multi_gene',
            datasets=['HepG2_1', 'HFFc6_1'],
            tpm=True,
            order='TSS')

# generate reports for the given gene
# with the given prefix
# display abundance information as a heatmap
# include differential transcript expression test significance results
# use significance threshold q >= 0.1
# use the browser-style representation as the transcript visualization
sg.gen_report('ADRM1', prefix='figures/adrm1',
             heatmap=True,
             include_qvals=True, q=0.1,
             browser=True)












