# swan\_vis.SwanGraph\(\)

## swan.Swangraph\(\)

```python
SwanGraph(file=None)
```

A graph class to represent a transcriptome and perform plotting and analysis from it

```text
Attributes:

        datasets (list of str):
                Names of datasets in the Graph
        counts (list of str):
                Names of columns holding counts in the Graph
        tpm (list of str):
                Names of columns holding tpm values in the Graph
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
        deg_test (pandas DataFrame): 
                A summary table of the results of a differential gene
                expression test
        deg_test_groups (list of str, len 2):
                The configuration of groupings used to run the differential
                gene expression test
                det_test (pandas DataFrame): 
                A summary table of the results of a differential transcript
                expression test
        det_test_groups (list of str, len 2):
                The configuration of groupings used to run the differential
                transcript expression test
```

## Methods

### add\_abundance

```python
add_abundance(self, counts_file, count_cols, dataset_name, tid_col='annot_transcript_id')
```

Adds abundance information to an existing dataset in the SwanGraph.

```text
Parameters:

        counts_file (str): Path to tsv counts matrix
        count_cols (str or list of str): Column names in counts_file to use
        dataset_name (str): Name of SwanGraph dataset to associate abundance with
        tid_col (str): Column name in counts_file containing transcript id
                Default='annot_transcript_id'
```

### add\_annotation

```python
add_annotation(self, fname)
```

Adds an annotation from input fname to the SwanGraph.

```text
Parameters:
        fname (str): Path to annotation GTF
```

### add\_dataset

```python
add_dataset(self, col, fname, dataset_name=None, whitelist=None, annot=None, counts_file=None, count_cols=None, tid_col='annot_transcript_id', include_isms=False)
```

Add transcripts from a dataset from either a GTF or a TALON database.

```text
Parameters:

    col (str): Name of column to add data to in the SwanGraph
    fname (str): Path to GTF or TALON db

    # Only for loading from TALON
    dataset_name (str): Dataset name in TALON db to add transcripts from
            Default=None
    whitelist (str): TALON whitelist of transcripts to add.
            Default: None
    annot (str): TALON annotation name in database to 
            add transcripts from
            Default: None

    # Only if also adding abundance
    counts_file (str): Path to tsv counts matrix
            Default=None
    count_cols (str or list of str): Column names in counts_file to use
            Default=None
    tid_col (str): Column name in counts_file containing transcript id
            Default='annot_transcript_id'

    include_isms (bool): Include ISMs from input dataset
            Default=False
```

### de\_gene\_test

```python
de_gene_test(self, dataset_groups)
```

Runs a differential expression test on the gene level.

```text
Parameters:

        dataset_groups (list of list of str, len 2): Grouping of datasets 
                from the SwanGraph to be used in the differential
                expression test
                Example: [['data1','data2'],['data3','data4']]

Returns: 

        test (pandas DataFrame): A summary table of the differential
                expression test, including p and q-values, as well 
                as log fold change.
```

### de\_transcript\_test

```python
de_transcript_test(self, dataset_groups)
```

Runs a differential expression test on the transcript level.

```text
Parameters:

        dataset_groups (list of list of str, len 2): Grouping of datasets 
                from the SwanGraph to be used in the differential
                expression test
                Example: [['data1','data2'],['data3','data4']]

Returns: 

        test (pandas DataFrame): A summary table of the differential
                expression test, including p and q-values, as well 
                as log fold change.
```

### find\_es\_genes

```python
find_es_genes(self)
```

Finds all unique genes containing novel exon skipping events. Requires that an annotation has been added to the SwanGraph.

```text
Returns:

        es_genes (list of str): A list of gene ids from the SwanGraph with 
                at least one novel exon skipping event
        es_transcripts (list of str): A list of transcript ids from the 
                SwanGraph with at least one novel exon skipping event
```

### find\_ir\_genes

```python
find_ir_genes(self)
```

Finds all unique genes containing novel intron retention events. Requires that an annotation has been added to the SwanGraph.

```text
Returns:

        ir_genes (list of str): A list of gene ids from the SwanGraph with 
                at least one novel intron retention event
        ir_transcripts (list of str): A list of transcript ids from the 
                SwanGraph with at least one novel intron retention event
```

### find\_isoform\_switching\_genes

```python
find_isoform_switching_genes(self, q=0.05, n_genes=None)
```

Finds isoform switching genes; genes that are not differentially expressed but contain at least one transcript that is. Requires that de\_gene\_test and de\_transcript\_test have been run.

```text
Parameters:

        q (float): q-value threshold to declare a gene/transcript 
                as significant
                Default: 0.05
        n_genes (int): Number of results to return. 
                Default: None (returns all found significant)

Returns:

        genes (list of str): List of gene names that are categorized as 
                isoform switching
        switches (pandas DataFrame): Summary table of genes that are 
                categorized as isoform switching
```

### gen\_report

```python
gen_report(self, gids, prefix, datasets='all', dataset_groups=False, dataset_group_names=False, novelty=False, heatmap=False, tpm=False, include_qvals=False, q=0.05, include_unexpressed=False, indicate_dataset=False, indicate_novel=False, browser=False, order='expression', threads=1)
```

Generates a PDF report for a given gene or list of genes according to the user's input.

```text
Parameters: 

        gids (str or list of str): Gene ids or names to generate
                reports for
        prefix (str): Path and/or filename prefix to save PDF and
                images used to generate the PDF

        datasets (list of str): Datasets to include in the report
                Default: Include columns for all datasets
        dataset_groups (list of list of str): Datasets to average
                together in the report and display as one column
                Example: [['group1_1','group1_2'],['group2_1','group2_2']]
        dataset_group_names (list of str): Names to give each group 
                given by dataset_groups. Must be the same length as 
                dataset_groups
                Example: ['group1', 'group2']
                Default: Will assign numbers 1 through length(dataset_group)

        novelty (bool): Include a column to dipslay novelty type of
                each transcript. Requires that a TALON GTF or DB has 
                been used to load data in
                Default: False

        heatmap (bool): Display expression values in a heatmap
                format. Requires that abundance information has been 
                added to the SwanGraph
                Default: False
        tpm (bool): Display TPM value of each transcript/dataset 
                combination, instead of presence/absence of each 
                transcript. Requires that abundance information has
                been added to the SwanGraph
                Default:False

        include_qvals (bool): Display q-val of differential expression 
                for each transcript and bold entries found to be
                differentially expressed. Requires that de_transcript_test
                has been run, and that abundance information has been
                added to the SwanGraph
                Default: False
        q (float): Q-value significance threshold to use when 
                bolding transcripts if include_qvals = True.
                Default: 0.05

        include_unexpressed (bool): Add transcript entries to report
                that are not expressed in any input dataset.
                Default: False

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
                                 
        threads (int): Number of threads to use. Multithreading is 
                recommended when making a large number of gene reports.
                Default: 1. 
```

### get\_de\_genes

```python
get_de_genes(self, q=0.05, n_genes=None)
```

Subsets the differential gene expression test summary table based on a q-value cutoff. Requires that de\_gene\_test has already been run.

```text
Parameters:

        q (float): q-value threshold to declare a gene as significant
                Default: 0.05
        n_genes (int): Number of results to return. 
                Default: None (returns all found significant)

Returns:

        genes (list of str): List of gene names that pass the 
                significance threshold
        test (pandas DataFrame): Summary table of genes that pass the
                significance threshold
```

### get\_de\_transcripts

```python
get_de_transcripts(self, q=0.05, n_transcripts=None)
```

Subsets the differential transcript expression test summary table based on a q-value cutoff. Requires that de\_transcript\_test has already been run.

```text
Parameters:

        q (float): q-value threshold to declare a transcript as significant
                Default: 0.05
        n_transcripts (int): Number of results to return. 
                Default: None (returns all found significant)

Returns:

        tids (list of str): List of transcript ids that pass the 
                significance threshold
        test (pandas DataFrame): Summary table of transcripts that pass
                the significance threshold
```

### plot\_each\_transcript

```python
plot_each_transcript(self, tids, prefix, indicate_dataset=False, indicate_novel=False, browser=False)
```

Plot each input transcript and automatically save figures.

```text
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
```

### plot\_each\_transcript\_in\_gene

```python
plot_each_transcript_in_gene(self, gid, prefix, indicate_dataset=False, indicate_novel=False, browser=False)
```

Plot each transcript in a given gene and automatically save figures.

```text
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
```

### plot\_graph

```python
plot_graph(self, gid, indicate_dataset=False, indicate_novel=False, prefix=None)
```

Plots a gene summary SwanGraph for an input gene. Does not automatically save the figure by default!

```text
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
```

### plot\_transcript\_path

```python
plot_transcript_path(self, tid, indicate_dataset=False, indicate_novel=False, browser=False, prefix=None)
```

Plots a path of a single transcript isoform through a gene summary SwanGraph.

```text
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
```

### save\_graph

```python
save_graph(self, prefix)
```

Saves the current SwanGraph in pickle format with the .p extension.

```text
Parameters: 

        prefix (str): Path and filename prefix. Resulting file will 
                be saved as prefix.p
```

