# Getting started

First, if you haven't already, make sure to [install Swan](../#installation). After installing, you'll be able to run Swan from Python.

Then, download the data and the reference transcriptome annotation from [here](https://hpc.oit.uci.edu/~freese/swan_files/). The bash commands to do so are given below.

Swan offers two main ways for loading transcriptomes. You can either load models from [GTFs](getting_started.md#adding-transcript-models-gtf-and-abundance-information-at-the-same-time), or from a [TALON db](getting_started.md#adding-transcript-models-talon-db-and-abundance-information)

Table of contents

* [Example data download](getting_started.md#download-example-data)
* [Starting and initializing your SwanGraph](getting_started.md#starting-up-swan-and-initializing-your-swangraph)
* [Add transcript models \(GTF\) and abundance info](getting_started.md#adding-transcript-models-gtf-and-abundance-information-at-the-same-time)
* [Saving and loading your SwanGraph](getting_started.md#saving-and-loading-your-swangraph)
* [Adding transcript models \(GTF\) and abundance info separately](getting_started.md#adding-transcript-models-gtf-and-abundance-information-separately)
* [Adding transcript models \(TALON db\) and abundance info](getting_started.md#adding-transcript-models-talon-db-and-abundance-information)

## Download example data

```bash
mkdir data
mkdir figures
cd data/

# gencode v29 human annotation
wget https://hpc.oit.uci.edu/~freese/swan_files/gencode.v29.annotation.gtf

# hepg2 data
wget https://hpc.oit.uci.edu/~freese/swan_files/hepg2_1_talon.gtf
wget https://hpc.oit.uci.edu/~freese/swan_files/hepg2_2_talon.gtf

# hffc6 data
wget https://hpc.oit.uci.edu/~freese/swan_files/hffc6_1_talon.gtf
wget https://hpc.oit.uci.edu/~freese/swan_files/hffc6_2_talon.gtf
wget https://hpc.oit.uci.edu/~freese/swan_files/hffc6_3_talon.gtf

# abundance file
wget https://hpc.oit.uci.edu/~freese/swan_files/all_talon_abundance_filtered.tsv

# example saved SwanGraph
wget https://hpc.oit.uci.edu/~freese/swan_files/swan.p

# talon database
wget https://hpc.oit.uci.edu/~freese/swan_files/talon.db

# talon whitelists
wget https://hpc.oit.uci.edu/~freese/swan_files/hepg2_whitelist.csv
wget https://hpc.oit.uci.edu/~freese/swan_files/hffc6_whitelist.csv

cd ../
```

## Starting up Swan and initializing your SwanGraph

```python
import swan_vis as swan

annot_gtf = 'data/gencode.v29.annotation.gtf'
hep_1_gtf = 'data/hepg2_1_talon.gtf'
hep_2_gtf = 'data/hepg2_2_talon.gtf'
hff_1_gtf = 'data/hffc6_1_talon.gtf'
hff_2_gtf = 'data/hffc6_2_talon.gtf'
hff_3_gtf = 'data/hffc6_3_talon.gtf'
ab_file = 'data/all_talon_abundance_filtered.tsv'
talon_db = 'data/talon.db'
```

Initialize an empty SwanGraph and add the transcriptome annotation to the SwanGraph.

```python
sg = swan.SwanGraph()
sg.add_annotation(annot_gtf)
```

## Adding transcript models \(GTF\) and abundance information at the same time

Add each dataset to the SwanGraph, along with the corresponding abundance information from the abundance matrix. The `count_cols` variable refers to the column name in the abundance file that corresponds to the counts for the input dataset.

```python
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
```

## Saving and loading your SwanGraph

Following this, you can save your SwanGraph so you can easily work with it again without re-adding all the data.

```python
sg.save_graph('swan')
```

And you can reload the graph again.

```python
sg = swan.SwanGraph('swan.p')
```

## Adding transcript models \(GTF\) and abundance information separately

Swan can also run without abundance information, although many of Swan's analysis functions depend on abundance information. To load just the transcript models, simply just leave out the `counts_file` and `count_cols` arguments to the `add_dataset()` function as shown below.

```python
sg = swan.SwanGraph()
sg.add_annotation(annot_gtf)
sg.add_dataset('HepG2_1', hep_1_gtf)
sg.add_dataset('HepG2_2', hep_2_gtf)
sg.add_dataset('HFFc6_1', hff_1_gtf)
sg.add_dataset('HFFc6_2', hff_2_gtf)
sg.add_dataset('HFFc6_3', hff_3_gtf)
```

If you have just added transcript models to the graph via `add_dataset()` and wish to add abundance information, this can be done using the `add_abundance()` function as seen below. Here, the string passed to `count_cols` is the column in the abundance file that corresponds to the dataset, and the argument passed to `dataset_name` is the name of the dataset that has already been added to the SwanGraph in the previous code block.

```python
sg.add_abundance(ab_file, count_cols='hepg2_1', dataset_name='HepG2_1')
sg.add_abundance(ab_file, count_cols='hepg2_2', dataset_name='HepG2_2')
sg.add_abundance(ab_file, count_cols='hffc6_1', dataset_name='HFFc6_1')
sg.add_abundance(ab_file, count_cols='hffc6_2', dataset_name='HFFc6_2')
sg.add_abundance(ab_file, count_cols='hffc6_3', dataset_name='HFFc6_3')
```

## Adding transcript models \(TALON db\) and abundance information

Swan is also directly compatible with TALON databases and can pull transcript models directly from them.

```python
sg = swan.SwanGraph()
sg.add_annotation(annot_gtf)

hepg2_whitelist='data/hepg2_whitelist.csv'
hffc6_whitelist='data/hffc6_whitelist.csv'
```

```python
sg.add_dataset('HepG2_1', talon_db,
    dataset_name='hepg2_1',
    whitelist=hepg2_whitelist,
    counts_file=ab_file,
    count_cols='hepg2_1')
sg.add_dataset('HepG2_2', talon_db,
    dataset_name='hepg2_2',
    whitelist=hepg2_whitelist,
    counts_file=ab_file,
    count_cols='hepg2_2')

sg.add_dataset('Hffc6_1', talon_db,
    dataset_name='hffc6_1',
    whitelist=hffc6_whitelist,
    counts_file=ab_file,
    count_cols='hffc6_1')
sg.add_dataset('Hffc6_2', talon_db,
    dataset_name='hffc6_2',
    whitelist=hffc6_whitelist,
    counts_file=ab_file,
    count_cols='hffc6_2')
sg.add_dataset('Hffc6_3', talon_db,
    dataset_name='hffc6_3',
    whitelist=hffc6_whitelist,
    counts_file=ab_file,
    count_cols='hffc6_3')
```

