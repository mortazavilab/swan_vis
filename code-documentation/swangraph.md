# swan\_vis.SwanGraph\(\)

## swan.Swangraph\(\)


`SwanGraph(sc=False)`
:   A graph class to represent a transcriptome and perform
    plotting and analysis from it

```text    
    Attributes:

            datasets (list of str):
                    Names of datasets in the Graph
            annotation (bool):
                    Whether an annotation transcriptome has been added.
            abundance (bool):
                    Whether abundance data has been added to the SwanGraph
            loc_df (pandas DataFrame):
                    DataFrame of all unique observed genomic
                    coordinates in the transcriptome
            edge_df (pandas DataFrame):
                    DataFrame of all unique observed exonic or intronic
                    combinations of splice sites in the transcriptome
            t_df (pandas DataFrame):
                    DataFrame of all unique transcripts found
                    in the transcriptome
            pg (swan PlottedGraph):
                    The PlottedGraph holds the information from the most
                    recently made plot
            adata (anndata AnnData):
                    Annotated data object to hold transcript expression values
                    and metadata
            edge_adata (anndata AnnData):
                    Annotated data object to hold edge expression values and metadata
            tss_adata (anndata AnnData):
                    Annotated data object to hold TSS expression values and metadata
            tes_adata (anndata AnnData):
                    Annotated data object to hold TES expression values and metadata
```

  ### Methods

  `add_abundance(self, counts_file)`
  :   Adds abundance from a counts matrix to the SwanGraph. Transcripts in the
      SwanGraph but not in the counts matrix will be assigned 0 counts.
      Transcripts in the abundance matrix but not in the SwanGraph will not
      have expression added.

      Parameters:
              counts_file (str): Path to TSV expression file where first column is
                      the transcript ID and following columns name the added datasets and
                      their counts in each dataset, OR to a TALON abundance matrix.

  `add_annotation(self, fname, verbose=False)`
  :   Adds an annotation from input fname to the SwanGraph.

      Parameters:
              fname (str): Path to annotation GTF
              verbose (bool): Display progress
                      Default: False


  `add_metadata(self, fname, overwrite=False)`
  :   Adds metadata to the SwanGraph from a tsv.

      Parameters:
              fname (str): Path / filename of tab-separated input file to add as
                      metadata for the datasets in the SwanGraph. Must contain column
                      'dataset' which contains dataset names that match those already
                      in the SwanGraph.
              overwrite (bool): Whether or not to overwrite duplicated columns
                      already present in the SwanGraph.
                      Default: False

  `add_multi_groupby(self, groupby)`
  :   Adds a groupby column that is comprised of multiple other columns. For
      instance, if 'sex' and 'age' are already in the obs table, add an
      additional column that's comprised of sex and age.

      Parameters:
              groupby (list of str): List of column names to turn into a multi
                      groupby column

  `add_transcriptome(self, fname, pass_list=None, include_isms=False, verbose=False)`
  :   Adds a whole transcriptome from a set of samples. 

      Parameters:
              fname (str): Path to GTF or TALON db
              pass_list (str): Path to pass list file (if passing a TALON DB)
              include_isms (bool): Include ISMs from input dataset
                      Default: False
              verbose (bool): Display progress
                      Default: False

  `de_gene_test(self, obs_col, obs_conditions=None)`
  :   Runs a differential expression test on the gene level.

      Parameters:
              obs_col (str): Metadata column from self.adata.obs to perform
                      the test on
              obs_conditions (list of str, len 2): Categories from 'obs_col' in
                      self.adata.obs to perform the test on.

      Returns:
              test (pandas DataFrame): A summary table of the differential
                      expression test, including p and q-values, as well
                      as log fold change.

  `de_transcript_test(self, obs_col, obs_conditions=None)`
  :   Runs a differential expression test on the transcript level.

      Parameters:
              obs_col (str): Metadata column from self.adata.obs to perform
                      the test on
              obs_conditions (list of str, len 2): Categories from 'obs_col' in
                      self.adata.obs to perform the test on.

      Returns:
              test (pandas DataFrame): A summary table of the differential
                      expression test, including p and q-values, as well
                      as log fold change.

  `die_gene_test(self, kind='iso', obs_col='dataset', obs_conditions=None, rc_thresh=10, verbose=False)`
  :   Finds genes with differential isoform expression between two conditions
      that are in the obs table. If there are more than 2 unique values in
      `obs_col`, the specific categories must be specified in `obs_conditions`

      Parameters:
              obs_col (str): Column name from self.adata.obs table to group on.
                      Default: 'dataset'
              obs_conditions (list of str, len 2): Which conditions from obs_col
                      to compare? Required if obs_col has more than 2 unqiue values.
              rc_thresh (int): Number of reads required for each conditions
                      in order to test the gene.
                      Default: 10
              verbose (bool): Display progress

      Returns:
              test (pandas DataFrame): A summary table of the differential
                      isoform expression test, including p-values and adjusted
                      p-values, as well as change in percent isoform usage (dpi) for
                      all tested genes.

  `find_es_genes(self, verbose=False)`
  :   Finds all unique genes containing novel exon skipping events.
      Requires that an annotation has been added to the SwanGraph.

      Parameters:
              verbose (bool): Display output

      Returns:
              es_df (pandas DataFrame): DataFrame detailing discovered novel
                      exon-skipping edges and the transcripts and genes they
                      come from

  `find_ir_genes(self, verbose=False)`
  :   Finds all unique genes containing novel intron retention events.
      Requires that an annotation has been added to the SwanGraph.

      Parameters:
              verbose (bool): Display output

      Returns:
              ir_df (pandas DataFrame): DataFrame detailing discovered novel
                      intron retention edges and the transcripts and genes they
                      come from

  `gen_report(self, gid, prefix, datasets=None, groupby=None, metadata_cols=None, novelty=False, layer='tpm', cmap='Spectral_r', include_qvals=False, q=0.05, log2fc=1, qval_obs_col=None, qval_obs_conditions=None, include_unexpressed=False, indicate_novel=False, display_numbers=False, transcript_name=False, browser=False, order='expression')`
  :   Generates a PDF report for a given gene or list of genes according
      to the user's input.

      Parameters:
              gid (str): Gene id or name to generate
                      reports for
              prefix (str): Path and/or filename prefix to save PDF and
                      images used to generate the PDF
              datasets (dict of lists): Dictionary of {'metadata_col':
                      ['metadata_category_1', 'metadata_category_2'...]} to represent
                      datasets and their order to include in the report.
                      Default: Include columns for all datasets / groupby category
              groupby (str): Column in self.adata.obs to group expression
                      values by
                      Default: None
              metadata_cols (list of str): Columns from metadata tables to include
                      as colored bars. Requires that colors have been set using
                      set_metadata_colors
              novelty (bool): Include a column to dipslay novelty type of
                      each transcript. Requires that a TALON GTF or DB has
                      been used to load data in
                      Default: False
              layer (str): Layer to plot expression from. Choose 'tpm' or 'pi'
              cmap (str): Matplotlib color map to display heatmap values
                      in.
                      Default: 'Spectral_r'
              include_qvals (bool): Display q-val of differential expression
                      for each transcript and bold entries found to be
                      differentially expressed. Requires that de_transcript_test
                      has been run, and that abundance information has been
                      added to the SwanGraph
                      Default: False
              q (float): Q-value significance threshold to use when
                      bolding transcripts if include_qvals = True.
                      Default: 0.05
              log2fc (float): Log2fc significance threshold to use when
                      bolding transcripts if include_qvals = True
              qval_obs_col (str): Metadata column from self.adata
              include_unexpressed (bool): Add transcript entries to report
                      that are not expressed in any input dataset.
                      Default: False
              indicate_novel (bool): Emphasize novel nodes and edges by
                      outlining them and dashing them respectively
                      Incompatible with indicate_dataset
                      Default: False
              browser (bool): Plot transcript models in genome browser-
                      style format. Incompatible with indicate_dataset and
                      indicate_novel
              display_numbers (bool): Display TPM or pi values atop each cell
                      Default: False
              order (str): Order to display transcripts in the report.
                      Options are
                              'tid': alphabetically by transcript ID
                              'expression': cumulative expression from high to low
                                      Requires that abundance information has been
                                      added to the SwanGraph
                              'tss': genomic coordinate of transcription start site
                              'tes': genomic coordinate of transcription end site
                      Default: 'expression' if abundance information is present,
                                       'tid' if not

  `get_de_genes(self, obs_col, obs_conditions=None, q=0.05, log2fc=1)`
  :   Subsets the differential gene expression test summary table based
      on a q-value and log2fc cutoff. Requires that de_gene_test has
      already been run.

      Parameters:
              obs_col (str): Column name from self.adata.obs table to group on.
                      Default: 'dataset'
              obs_conditions (list of str, len 2): Which conditions from obs_col
                      to compare? Required if obs_col has more than 2 unique values.
              q (float): q-value threshold to declare a gene as significant
                      Default: 0.05
              log2fc (float): log2fc threshold to declare a gene significant
                      Default: 1

      Returns:
              test (pandas DataFrame): Summary table of genes that pass the
                      significance threshold

  `get_de_transcripts(self, obs_col, obs_conditions=None, q=0.05, log2fc=1)`
  :   Subsets the differential gene expression test summary table based
      on a q-value and log2fc cutoff. Requires that de_gene_test has
      already been run.

      Parameters:
              obs_col (str): Column name from self.adata.obs table to group on.
                      Default: 'dataset'
              obs_conditions (list of str, len 2): Which conditions from obs_col
                      to compare? Required if obs_col has more than 2 unique values.
              q (float): q-value threshold to declare a gene as significant
                      Default: 0.05
              log2fc (float): log2fc threshold to declare a gene significant
                      Default: 1

      Returns:
              test (pandas DataFrame): Summary table of genes that pass the
                      significance threshold

  `get_die_genes(self, kind='iso', obs_col='dataset', obs_conditions=None, p=0.05, dpi=10)`
  :   Filters differential isoform expression test results based on adj.
      p-value and change in percent isoform usage (dpi).

      Parameters:
              kind (str): Choose from 'iso', 'tss' or 'tes'
                      Default: 'iso'
              obs_col (str): Column name from self.adata.obs table to group on.
                      Default: 'dataset'
              obs_conditions (list of str, len 2): Which conditions from obs_col
                      to compare? Required if obs_col has more than 2 unique values.
              p (float): Adj. p-value threshold to declare a gene as isoform
                      switching / having DIE.
                      Default: 0.05
              dpi (float): DPI (in percent) value to threshold genes with DIE
                      Default: 10

      Returns:
              test (pandas DataFrame): Summary table of genes that pass
                      the significance threshold

  `get_edge_abundance(self, prefix=None, kind='counts')`
  :   Gets edge expression from the current SwanGraph in a DataFrame
      complete information about where edge is.

      Parameters:
              prefix (str): Path and filename prefix. Resulting file will
                      be saved as prefix_edge_abundance.tsv
                      Default: None (will not save)
              kind (str): Choose "tpm" or "counts"

      Returns:
              df (pandas DataFrame): Abundance and metadata information about
                      each edge.

  `get_tes_abundance(self, prefix=None, kind='counts')`
  :   Gets TES expression from the current SwanGraph in a DataFrame
      complete information about where TES is.

      Parameters:
              prefix (str): Path and filename prefix. Resulting file will
                      be saved as prefix_tes_abundance.tsv
                      Default: None (will not save)
              kind (str): Choose "tpm" or "counts"

      Returns:
              df (pandas DataFrame): Abundance and metadata information about
                      each TSS.

  `get_tss_abundance(self, prefix=None, kind='counts')`
  :   Gets TSS expression from the current SwanGraph in a DataFrame with
      complete information about where TSS is.

      Parameters:
              prefix (str): Path and filename prefix. Resulting file will
                      be saved as prefix_tss_abundance.tsv
              kind (str): Choose "tpm" or "counts"

      Returns:
              df (pandas DataFrame): Abundance and metadata information about
                      each TSS.


  `plot_each_transcript(self, tids, prefix, indicate_dataset=False, indicate_novel=False, browser=False)`
  :   Plot each input transcript and automatically save figures

      Parameters:
              tids (list of str): List of transcript ids to plot
              prefix (str): Path and file prefix to automatically save
                      the plotted figures
              indicate_dataset (str): Dataset name from SwanGraph to
                      highlight with outlined nodes and dashed edges
                      Incompatible with indicate_novel
                      Default: False (no highlighting)
              indicate_novel (bool): Highlight novel nodes and edges by
                      outlining them and dashing them respectively
                      Incompatible with indicate_dataset
                      Default: False
              browser (bool): Plot transcript models in genome browser-
                      style format. Incompatible with indicate_dataset and
                      indicate_novel

  `plot_each_transcript_in_gene(self, gid, prefix, indicate_dataset=False, indicate_novel=False, browser=False)`
  :   Plot each transcript in a given gene and automatically save figures

      Parameters:
              gid (str): Gene id or gene name to plot transcripts from
              prefix (str): Path and file prefix to automatically save
                      the plotted figures
              indicate_dataset (str): Dataset name from SwanGraph to
                      highlight with outlined nodes and dashed edges
                      Incompatible with indicate_novel
                      Default: False (no highlighting)
              indicate_novel (bool): Highlight novel nodes and edges by
                      outlining them and dashing them respectively
                      Incompatible with indicate_dataset
                      Default: False
              browser (bool): Plot transcript models in genome browser-
                      style format. Incompatible with indicate_dataset and
                      indicate_novel

  `plot_graph(self, gid, indicate_dataset=False, indicate_novel=False, prefix=None)`
  :   Plots a gene summary SwanGraph for an input gene.
      Does not automatically save the figure by default!

      Parameters:
              gid (str): Gene ID to plot for (can also be gene name but
                      we've seen non-unique gene names so use at your own risk!)
              indicate_dataset (str): Dataset name from SwanGraph to
                      highlight with outlined nodes and dashed edges
                      Incompatible with indicate_novel
                      Default: False (no highlighting)
              indicate_novel (bool): Highlight novel nodes and edges by
                      outlining them and dashing them respectively
                      Incompatible with indicate_dataset
                      Default: False
              prefix (str): Path and file prefix to automatically save
                      the plotted figure
                      Default: None, won't automatically save
              display (bool): Display the plot during runtime
                      Default: True

  `plot_transcript_path(self, tid, indicate_dataset=False, indicate_novel=False, browser=False, prefix=None, display=True)`
  :   Plots a path of a single transcript isoform through a gene summary
      SwanGraph.

      Parameters:
              tid (str): Transcript id of transcript to plot
              indicate_dataset (str): Dataset name from SwanGraph to
                      highlight with outlined nodes and dashed edges
                      Incompatible with indicate_novel
                      Default: False (no highlighting)
              indicate_novel (bool): Highlight novel nodes and edges by
                      outlining them and dashing them respectively
                      Incompatible with indicate_dataset
                      Default: False
              browser (bool): Plot transcript models in genome browser-
                      style format. Incompatible with indicate_dataset and
                      indicate_novel
              prefix (str): Path and file prefix to automatically save
                      the plotted figure
                      Default: None, won't automatically save
              display (bool): Display the plot during runtime
                      Default: True

  `save_graph(self, prefix)`
  :   Saves the current SwanGraph in pickle format with the .p extension

      Parameters:
              prefix (str): Path and filename prefix. Resulting file will
                      be saved as prefix.p

  `set_metadata_colors(self, obs_col, cmap)`
  :   Set plotting colors for datasets based on a column in the metadata
      table.

      Parameters:
              obs_col (str): Name of metadata column to set colors for
              cmap (dict): Dictionary of metadata value : color (hex code with #
                      character or named matplotlib color)

  `set_plotting_colors(self, cmap=None, default=False)`
  :   Set plotting colors for SwanGraph and browser models.

      Parameters:
              cmap (dict): Dictionary of metadata value : color (hex code with #
                      character or named matplotlib color). Keys should be a subset
                      or all of ['tss', 'tes', 'internal', 'exon', 'intron', 'browser']
                      Default: None
              default (bool): Whether to revert to default colors
