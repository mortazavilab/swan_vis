# Getting started
import swan_vis as swan
sg = swan.SwanGraph()
annot_gtf = '/dfs6/pub/freese/mortazavi_lab/ref/gencode.v29/gencode.v29.annotation.gtf'
data_gtf = 'data/all_talon_observedOnly.gtf'
ab_file = 'data/all_talon_abundance_filtered.tsv'
talon_db = 'data/talon.db'
pass_list = 'data/all_pass_list.csv'
meta = 'data/metadata.tsv'
sg.add_annotation(annot_gtf)
sg.add_transcriptome(data_gtf)
sg.add_abundance(ab_file)
sg.save_graph('swan')
sg = swan.read('swan.p')
sg = swan.SwanGraph()


sg.add_annotation(annot_gtf)
sg.add_transcriptome(talon_db, pass_list=pass_list)
sg.add_abundance(ab_file)
sg.add_metadata(meta)
sg.save_graph('data/swan')

# Analysis
sg = swan.read('data/swan.p')
obs_col = 'cell_line'
obs_conditions = ['hepg2', 'hffc6']

# perform a differential gene expression
# Wald test on the provided metadata column and conditions
test = sg.de_gene_test(obs_col, obs_conditions=obs_conditions)
test.head(2)

# deg - differential gene expression
uns_key = swan.make_uns_key('deg',
                            obs_col=obs_col,
                            obs_conditions=obs_conditions)
test = sg.adata.uns[uns_key]
test.head(2)

# return a table of significantly differentially-expressed genes
# for a given q val + log2fc threshold
de_genes = sg.get_de_genes(obs_col, obs_conditions=obs_conditions,
                           q=0.05, log2fc=1)

obs_col = 'cell_line'
obs_conditions = ['hepg2', 'hffc6']

# perform a differential transcript expression
# Wald test on the provided metadata column and conditions
test = sg.de_transcript_test(obs_col, obs_conditions=obs_conditions)


# In[43]:


# save the SwanGraph as a Python pickle file
sg.save_graph('data/swan')
sg.save_graph('swan')
sg.save_graph('data/swan_files_full')
# sg.save_graph('data/swan_back')


# The output in `test` is a summary table for the differential expression test.

# In[8]:


test.head(2)


# The results of this test are similarly stored in an automatically-generated key in `sg.adata.uns`, and will be saved to the SwanGraph if you save it. You can regenerate this key and access the summary table by running the following code:

# In[9]:


# det - differential transcript expression
uns_key = swan.make_uns_key('det',
                            obs_col=obs_col,
                            obs_conditions=obs_conditions)
test = sg.adata.uns[uns_key]
test.head(2)


# Again, Swan can also automatically subset the test summary table to pull out genes that pass a certain significance threshold.

# In[10]:


# return a table of significantly differentially-expressed genes
# for a given q val + log2fc threshold
de_transcripts = sg.get_de_transcripts(obs_col, obs_conditions=obs_conditions,
                           q=0.05, log2fc=1)


# In[11]:


de_transcripts.head()


# ## <a name="is"></a>Isoform switching / Differential isoform expression testing

# Isoform switching / differential isoform expression (DIE) testing is implemented according to the strategy in [Joglekar et. al., 2021](https://www.nature.com/articles/s41467-020-20343-5). DIE can roughly be described as finding statistically significant changes in isoform expression between two conditions along with a change in percent isoform usage per gene.
#
# Pairwise comparisons can be set up using different columns in the metadata that was added to the SwanGraph with the `obs_col` and `obs_conditions` arguments.

# In[12]:


# look at valid metadata options
sg.adata.obs


# In[13]:


# find genes that exhibit DIE between HFFc6 and HepG2
obs_col = 'cell_line'
obs_conditions = ['hepg2', 'hffc6']
die_table = sg.die_gene_test(obs_col=obs_col,
                             obs_conditions=obs_conditions,
                             verbose=True)


# The resultant table contains an entry for each gene with the p value (`p_val`), adjusted p value (`adj_p_val`), and change in percent isoform usage for the top two isoforms (`dpi`). Exact details on these calculations can be found in [Joglekar et. al., 2021](https://www.nature.com/articles/s41467-020-20343-5).

# In[14]:


die_table.head()


# As with differential expression testing, differential isoform expression testing results are stored automatically in `sg.adata.uns`, and will be saved to the SwanGraph if you save it. You can regenerate this key and access the summary table by running the following code:

# In[15]:


# die_iso - isoform level differential isoform expression test results
uns_key = swan.make_uns_key('die',
                            obs_col=obs_col,
                            obs_conditions=obs_conditions)
test = sg.adata.uns[uns_key]
test.head(2)


# Swan comes with an easy way to filter your DIE test results based on adjusted p value and dpi thresholds.

# In[16]:


test = sg.get_die_genes(obs_col=obs_col, obs_conditions=obs_conditions,
                       p=0.05, dpi=10)
test.head()


# Swan also now automatically tracks transcription start site (TSS) and transcription end site (TES) usage, and find genes that exhibit DIE on the basis of their starts or ends. To do this, use the `kind` argument to `die_gene_test`.

# In[17]:


# find genes that exhibit DIE for TSSs between HFFc6 and HepG2
die_table = sg.die_gene_test(kind='tss',
                             obs_col=obs_col,
                             obs_conditions=obs_conditions,
                             verbose=True)
die_table.head()


# To access the die results on the tss level, use `die_kind='tss'` as input to `make_uns_key()`.

# In[18]:


# die_iso - TSS level differential isoform expression test results
uns_key = swan.make_uns_key('die',
                            obs_col=obs_col,
                            obs_conditions=obs_conditions,
                            die_kind='tss')
test = sg.adata.uns[uns_key]
test.head(2)


# And provide the `kind='tss'` option to `get_die_genes()` when trying to filter your test results.

# In[19]:


test = sg.get_die_genes(kind='tss', obs_col=obs_col,
                        obs_conditions=obs_conditions,
                        p=0.05, dpi=10)
test.head()


# For TESs, use `kind='tes'` as input to `die_genes_test()`, `die_kind='tes'` to `make_uns_key()`, and `kind='tes'` to `get_die_genes()`.

# In[20]:


# find genes that exhibit DIE for TESs between HFFc6 and HepG2
die_table = sg.die_gene_test(kind='tes', obs_col='cell_line', obs_conditions=['hepg2', 'hffc6'])
die_table.head()


# In[21]:


# die_iso - TSS level differential isoform expression test results
uns_key = swan.make_uns_key('die',
                            obs_col=obs_col,
                            obs_conditions=obs_conditions,
                            die_kind='tes')
test = sg.adata.uns[uns_key]
test.head(2)


# In[22]:


test = sg.get_die_genes(kind='tes', obs_col=obs_col,
                        obs_conditions=obs_conditions,
                        p=0.05, dpi=10)
test.head()



sg.adata.obs.head()

col_name = sg.add_multi_groupby(['cell_line', 'replicate'])

print(col_name)
print(sg.adata.obs.head())

obs_col = col_name
obs_conditions = ['hffc6_3', 'hepg2_1']

deg_summary = sg.de_gene_test(obs_col=obs_col,
                              obs_conditions=obs_conditions)
det_summary = sg.de_transcript_test(obs_col=obs_col,
                                    obs_conditions=obs_conditions)
die_summary = sg.die_gene_test(obs_col=obs_col,
                               obs_conditions=obs_conditions)


# returns a DataFrame of genes, transcripts, and specific edges in
# the SwanGraph with novel exon skipping events
es_df = sg.find_es_genes(verbose=True)


es_df.head()


# returns a DataFrame of genes, transcripts, and specific edges in
# the SwanGraph with novel intron retaining events
ir_df = sg.find_ir_genes(verbose=True)




# save the SwanGraph as a Python pickle file
sg.save_graph('data/swan')

ir_df.head()

# visualization

# this is just to display saved images
# from IPython.display import display, Image


# In[26]:


get_ipython().run_line_magic('matplotlib', 'inline')
import swan_vis as swan

# code to download this data is in the Getting started tutorial
sg = swan.read('data/swan.p')


# ## <a name="gene_summary"></a>Gene summary graphs
# Gene summary graphs display a high-level summary of the complexity of splicing within a certain gene locus. All observed splice sites and splice junctions from input datasets, and the annotation (if added) are plotted in full color. Genes can be called to be plotted either using their gene ids or gene names, but we recommend using gene ids as we have encountered redundant gene names during testing.

# In[27]:


# plot a gene summary graph for the given gene
sg.plot_graph('ADRM1')


# Gene summary graphs are a type of graph plot and therefore have plotting options that allow the user to highlight nodes and edges that are either not present in the annotation with the `indicate_novel=True` option.

# For instance, say we want to highlight the nodes and edges that are not seen in the annotation. In this representation, nodes (splice sites) that are not present in the annotation will appear outlined, and intronic or exonic regions (edges) that are not present in the annotation will appear dashed.

# In[28]:


# plot a gene summary graph with novel splice sites and
# splice junctions plotted as outlined nodes and dashed edges respectively
sg.plot_graph('ADRM1', indicate_novel=True)


# <!-- Similarly, you can highlight the nodes that come from a specific dataset. Outlined nodes and dashed edges are those that are present in the queried dataset. -->

# In[29]:


# # plot a gene summary graph without displaying it
# sg.plot_graph('ADRM1', indicate_novel=True, display=False)


# ## <a name="transcript_path"></a>Transcript path graphs
#

# Transcript path graphs display the same structure as gene summary graphs but gray out nodes and edges (splice sites and intronic/exonic regions) that are not present in the given transcript. In this case, the transcript id field is needed to plot the path.

# In[30]:


# plot the path of a specific transcript through its parent gene
# summary graph for a given transcript
sg.plot_transcript_path('TALONT000301961')


# You can use the `indicate_novel=True` option for transcript path Swan graphs too.

# In[31]:


# plot the path of a specific transcript through its parent gene summary
# graph for a given transcript
# plot novel splice sites as outlined nodes
# plot novel splice junctions as dashed edges
sg.plot_transcript_path('TALONT000301961', indicate_novel=True)


# In[32]:


# # plot the path of a specific transcript through its parent gene summary
# # graph for a given transcript
# # plot splice sites from the given dataset as outlined nodes
# # plot splice junctions from the given dataset as dashed edges
# sg.plot_transcript_path('TALONT000301953', indicate_dataset='HFFc6_1')


# For transcripts, there is also a unique option that allows you to generate the genome-browser style representation of a transcript, using the `browser` option.

# In[33]:


# plot the traditional browser-style representation
# for a given transcript
sg.plot_transcript_path('TALONT000301961', browser=True)


# ## <a name="change_colors"></a>Changing colors

# If you are unhappy with the default color scheme for the transcript models (both browser and Swan style), Swan now includes an option to change the colors. Use the `set_plotting_colors()` function using a dictionary that maps the different components of the transcript model (tss, tes, internal, exon, intron, browser) to either a [named Matplotlib colors](https://matplotlib.org/stable/gallery/color/named_colors.html) or hex code.
#
# The user doesn't have to supply colors for all of the components. For instance, say we just want to change the color of the introns to purple instead of pink.

# In[34]:


cmap = {'intron': 'purple'}
sg.set_plotting_colors(cmap)


# In[35]:


# plot the path of a specific transcript through its parent gene
# summary graph for a given transcript
sg.plot_transcript_path('TALONT000301961')


# Of course, you can also change all of the colors too! Beware that Swan will automatically compute a grayed-out version of the color. Try to choose colors that all have similar low "lightness" in HSL space so that the grayed-out components in the transcript path plots are still distinguishable from one another.

# In[36]:


cmap = {'intron': 'rebeccapurple', 'exon': 'chartreuse', 'tss': 'dodgerblue',
        'tes': 'salmon', 'internal': 'goldenrod'}
sg.set_plotting_colors(cmap=cmap)
# plot the path of a specific transcript through its parent gene
# summary graph for a given transcript
sg.plot_transcript_path('TALONT000301961')


# You can also change the color of plotted browser figures at the same time by adding a `browser` color to your colormap.

# In[37]:


cmap = {'browser': 'palevioletred'}
sg.set_plotting_colors(cmap=cmap)
# plot the path of a specific transcript in browser format
sg.plot_transcript_path('TALONT000301961', browser=True)


# If at any point you want to revert to the default color settings, you can run the following.

# In[38]:


sg.set_plotting_colors(default=True)
sg.plot_transcript_path('TALONT000301961')
sg.plot_transcript_path('TALONT000301961', browser=True)


# ## <a name="save_fig"></a>Saving a figure

# Saving a figure in Swan for `plot_graph` and for `plot_transcript_path` can be done in two different ways.

# The first way involves calling `save_fig` after your figure has been generated. This method allows you to give your figure whatever name you want.

# In[39]:


# plot gene summary graph for a given gene
sg.plot_graph('ADRM1')

# save the currently-plotted figure with the given filename and location
swan.save_fig('figures/my_gene_summary.png')


# In[40]:


# display(Image(filename='figures/my_gene_summary.png'))


# The second way only requires one line of code and requires that the user pass the corresponding plotting function a `prefix` for the filename and path. The file will be automatically named according to the settings in the graph.

# In[41]:


# plot a gene summary graph with novel splice sites and junctions
# outlined and dashed respectively
# save the figure with the given prefix
sg.plot_graph('ADRM1', indicate_novel=True, prefix='figures/adrm1')


# In[42]:


# display(Image(filename='figures/adrm1_novel_ENSG00000130706.12_summary.png'))


# As you can see, here the gene name ADRM1 is not used to save the figure because we have encountered conflicing gene names in our internal use of Swan. To avoid these clashes, Swan automatically fetches the gene id associated with the first instance of the gene name it finds and uses it to save the graph, in the interest of not accidentally overwriting a preexisting file.

# ## <a name="swan_report"></a>Swan reports

# Swan reports display all the expressed transcripts in a given gene in a PDF format. There are many ways to customize these reports so here are a few. Unlike the above plotting options, the user must provide a `prefix` argument as there are many files that must be automatically generated to create the report.

# In[43]:


# generate a report for the given gene
# save it with the given filepath prefix
# include differential transcript expression test significance results
# (defaults to significance threshold q >= 0.05)
# display the novelty category associated with the transcript
# display novel splice sites and junctions
#     as outlined nodes and dashed edges respectively
sg.gen_report('ADRM1',
              prefix='figures/adrm1_paper',
              novelty=True,
              indicate_novel=True)


# In[44]:


# display(Image(filename='figures/adrm1_paper_novel_ENSG00000130706.12_report.png'))


# You can also generate a report where differentially-expressed transcripts from prior tests are bolded. Use the `include_qvals=True` option and tell Swan which differential transcript expression test results to display using `qval_obs_col` and `qval_obs_conditions`. Here, I'll show the results of the HFFc6 vs. HepG2 differential transcript expression test that was run in the Analysis tutorial.

# In[45]:


# generate a report for the given gene
# save it with the given filepath prefix
# include differential transcript expression test significance results
# (defaults to significance threshold q >= 0.05)
# display the novelty category associated with the transcript
# display novel splice sites and junctions
#     as outlined nodes and dashed edges respectively
# include the qvals from differential transcript expression test
# differential transcript expression test metadata variable: cell_line
# differential transcript expression test categories to compare: hepg2, hffc6
sg.gen_report('ADRM1',
              prefix='figures/adrm1_paper',
              novelty=True,
              indicate_novel=True,
              include_qvals=True,
              qval_obs_col='cell_line',
              qval_obs_conditions=['hepg2', 'hffc6'])


# In[46]:


# display(Image(filename='figures/adrm1_paper_novel_qval_ENSG00000130706.12_report.png'))


# You can also plot transcripts based on their percent isoform (pi) values using `layer=pi`, which help illustrate the basis on which isoform switches are called in Swan. Here I'm choosing to plot it in a separate color using the `cmap` argument to make it clear that the metric being plotted is different, and overlaying each cell with the pi value using `display_numbers=True`.
#
# For this example, I'll plot a gene that was called as isoform switching by the Swan isoform switching module, NIPAL3.

# In[47]:


# generate a report for the given gene
# save it with the given filepath prefix
# plot the percent isoform (pi) values
# use the magma color way
# display values on top of each cell
# display the novelty category associated with the transcript
# display novel splice sites and junctions
#     as outlined nodes and dashed edges respectively
sg.gen_report('NIPAL3',
              prefix='figures/nipal3',
              layer='pi',
              cmap='magma',
              display_numbers=True,
              novelty=True,
              indicate_novel=True)


# In[48]:


# display(Image(filename='figures/nipal3_novel_ENSG00000001461.16_report.png'))


# In these cases, it can be beneficial to actually look at the data in the original groups that the isoform switching test was performed. To group your input samples by a metadata column that can be found in `sg.adata.obs` in the report, use the `groupby` option. Here, I've also demonstrated that if you have transcript names in your GTF or TALON db, that those can be displayed instead of the transcript IDs using `transcript_name=True`.

# In[49]:


# generate a report for the given gene
# save it with the given filepath prefix
# plot the percent isoform (pi) values
# use the magma color way
# display values on top of each cell
# display the novelty category associated with the transcript
# display novel splice sites and junctions
#     as outlined nodes and dashed edges respectively
# group datasets based on the 'cell_line' metadata column
sg.gen_report('NIPAL3',
              prefix='figures/nipal3',
              layer='pi',
              cmap='magma',
              display_numbers=True,
              novelty=True,
              indicate_novel=True,
              groupby='cell_line',
              transcript_name=True)


# Using this strategy, the basis of the isoform switch is a little clearer. The longer isoform, NIPAL3-204, is proportionally higher-expressed in HFFc6, wherease NIPAL3-202 is proportionally higher-expressed in HepG2.

# In[50]:


# display(Image(filename='figures/nipal3_novel_cell_line_ENSG00000001461.16_report.png'))


# Swan now supports using colors to represent metadata categories which can be useful for more complex sets of samples. Swan cannot automatically resize dataset names for gene reports and therefore I recommend using this strategy when plotting a large number of datasets.

# First, assign colors to different metadata columns in `sg.adata.obs`. You can use hexcodes or [named Matplotlib colors](https://matplotlib.org/stable/gallery/color/named_colors.html).

# In[51]:


sg.set_metadata_colors('cell_line', {'hepg2': 'gold', 'hffc6': '#ba55d3'})


# Then use the `metadata_cols` option to indicate what colored metadata categories you'd like to plot at the top of the gene report. Here I'm also demonstrating the option to plot the browser-style transcript representation using the `browser=True` option.

# In[52]:


# generate a report for the given gene
# save it with the given filepath prefix
# plot the percent isoform (pi) values
# use the magma color way
# display values on top of each cell
# display the novelty category associated with the transcript
# display novel splice sites and junctions
#     as outlined nodes and dashed edges respectively
# group datasets based on the 'cell_line' metadata column
# color cell lines by metadata colors
# plot the genome browser representation of the transcript models
# include the qvals from differential transcript expression test
# differential transcript expression test metadata variable: cell_line
# differential transcript expression test categories to compare: hepg2, hffc6
sg.gen_report('ADRM1',
              prefix='figures/adrm1',
              layer='pi',
              cmap='magma',
              display_numbers=True,
              novelty=True,
              groupby='cell_line',
              transcript_name=True,
              metadata_cols=['cell_line'],
              browser=True,
              include_qvals=True,
              qval_obs_col='cell_line',
              qval_obs_conditions=['hepg2', 'hffc6'])


# In[53]:


# display(Image(filename='figures/adrm1_browser_color_cell_line_ENSG00000130706.12_report.png'))


# You can also include more than one metadata column to color.

# In[54]:


# from lighter to darker blue
sg.set_metadata_colors('replicate', {'1': '#bef4ff',
                                     '2': '#73a8b2',
                                     '3': '#263133'})


# In[55]:


# generate a report for the given gene
# save it with the given filepath prefix
# use the virids color way
# display values on top of each cell
# display the novelty category associated with the transcript
# color cell lines and replicates by metadata colors
# plot the genome browser representation of the transcript models
sg.gen_report('ADRM1',
              prefix='figures/adrm1',
              cmap='viridis',
              display_numbers=True,
              novelty=True,
              transcript_name=True,
              metadata_cols=['cell_line', 'replicate'],
              browser=True)


# In[56]:


# display(Image(filename='figures/adrm1_browser_color_replicate_cell_line_ENSG00000130706.12_report.png'))


# Note that if I try to use the `groupby='cell_line'` option with `metadata_cols=['cell_line', 'replicate']`, Swan will throw an error because there are multiple distinct replicates that belong to each cell line which makes the groupby impossible.

# In[57]:



# In[66]:


# display(Image(filename='figures/srf_ENSMUSG00000015605.6_report.png'))


# You can choose which datasets and the order in which they display in the report using the `datasets` option combined with different metadata. Suppose we want to restrict the datasets in the report to just the HFFc6 datasets, and we want to display the replicates 2 and 3 only, in descending order. Specify which categories from the relevant metadata columns you wish to include using a dictionary with the format `{metadata_column: [metadata_category1, metadata_category2 ...]}`. The output report will include the intersection of datasets that satisfy each condition.
#
# Here I've also demonstrated the `order` option, where you can order the transcripts based on transcript ID, expression level (default), or genomic location of TSS / TES (`[tid', 'expression', 'tss', 'tes']` respectively).

# In[67]:


# generate a report for the given gene
# save it with the given filepath prefix
# use the virids color way
# color cell lines by metadata colors
# restrict data shown to just the hffc6 cell line and replicates 3 and 2
# order transcripts based on genomic location of TSS
sg.gen_report('NIPAL3',
              prefix='figures/nipal3',
              cmap='viridis',
              metadata_cols=['cell_line', 'replicate'],
              datasets={'cell_line': 'hffc6', 'replicate': ['3', '2']},
              order='tss')


# As you can see here, the data displayed is limited to those belonging to hffc6 replicates 2 and 3, and we display replicates 3 and 2 in a specific order.

# In[68]:


# display(Image(filename='figures/nipal3_ENSG00000001461.16_report.png'))


# And here I'll show all the hffc6 replicates ordered by tes.

# In[69]:


# generate a report for the given gene
# save it with the given filepath prefix
# use the virids color way
# color cell lines and replicates by metadata colors
# restrict data shown to just the hffc6 datasets
# order transcripts based on genomic location of TES
# use browser-style representation
sg.gen_report('NIPAL3',
              prefix='figures/nipal3',
              cmap='viridis',
              metadata_cols=['cell_line', 'replicate'],
              datasets={'cell_line': 'hffc6'},
              order='tes',
              browser=True)


# In[70]:


# display(Image(filename='figures/nipal3_browser_hffc6_ENSG00000001461.16_report.png'))
