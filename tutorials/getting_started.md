#  Getting started: initializing, adding data, and saving your SwanGraph 

First, if you haven't already, make sure to [install Swan](../#installation).
After installing, you'll be able to run Swan from Python.

Then, download the data and the reference transcriptome annotation from [here](https://hpc.oit.uci.edu/~freese/swan_files/). The bash commands to do so are given below.

Swan offers two main ways for loading transcriptomes. You can either load models from [a properly-formatted GTF](getting_started.md#adding-transcript-models-gtf-and-abundance-information-at-the-same-time), or from a [TALON db](getting_started.md#adding-transcript-models-talon-db-and-abundance-information).
Please see the [input file format documentation](../faqs/file_formats.md) for specifics on how these files should be formatted.

We've provided four examples on different ways you can add data to your SwanGraph in the following tutorial. You only need to run one!
* [Using a GTF and abundance table together](getting_started.md#adding-transcript-models-gtf-and-abundance-information-at-the-same-time)
* [Using a GTF and abundance table separately](getting_started.md#adding-transcript-models-gtf-and-abundance-information-separately)
* [Using a TALON database and abundance table together](getting_started.md#adding-transcript-models-talon-db-and-abundance-information)
* [Batch adding datasets](getting_started.md#batch-adding-datasets)


Other sections: 
* [Example data download](getting_started.md#download-example-data)
* [Starting and initializing your SwanGraph](getting_started.md#starting-up-swan-and-initializing-your-swangraph)
* [Saving and loading your SwanGraph](getting_started.md#saving-and-loading-your-swangraph)

This page can also be read from top to bottom, just know that you may be running things more than once!

Running this tutorial (with only one of the dataset addition options) on my laptop took around 7 minutes and 5 GB of RAM. 

View and download the Jupyter Notebook for this tutorial [here](https://drive.google.com/file/d/1xLEQU_nWAPb5es4lOR4-99kftW6cu8SI/view?usp=sharing)

## Download example data

```bash
# # run this block in your bash terminal
# mkdir data
# mkdir figures
# cd data/

# # download files
# wget http://crick.bio.uci.edu/freese/swan_files.tgz
    
# # expand files 
# tar xzf swan_files.tgz
# mv swan_files/* .
# rm -r swan_files/

# # download reference annotation
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
# gunzip gencode.v29.annotation.gtf.gz

# cd ../
```

Alternatively, if you wish to use a smaller example dataset (just chr20), run the following. The rest of the commands in the tutorial will work with the smaller data too!

```bash
# # run this block in your bash terminal
# mkdir data
# mkdir figures
# cd data/

# # download files
# wget http://crick.bio.uci.edu/freese/swan_files_mini.tgz
    
# # expand files 
# tar xzf swan_files_mini.tgz
# mv swan_files_mini/* .
# rm -r swan_files_mini/

# cd ../
```

## Starting up Swan and initializing your SwanGraph

The rest of the code should be done in the Python shell, or run from a `.py` file. 


```python
import swan_vis as swan
```


```python
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
# initialize a new SwanGraph
sg = swan.SwanGraph() 

# add an annotation transcriptome 
sg.add_annotation(annot_gtf)
```

## Adding transcript models \(GTF\) and abundance information at the same time

Add each dataset to the SwanGraph, along with the corresponding abundance information from the abundance matrix. The `count_cols` variable refers to the column name in the abundance file that corresponds to the counts for the input dataset.

```python
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
```

## Saving and loading your SwanGraph

Following this, you can save your SwanGraph so you can easily work with it again without re-adding all the data.

```python
# save the SwanGraph as a Python pickle file
sg.save_graph('swan')
```

And you can reload the graph again.

```python
# load up a saved SwanGraph from a pickle file
sg = swan.SwanGraph('swan.p')
```

## Adding transcript models \(GTF\) and abundance information separately

Swan can also run without abundance information, although many of Swan's analysis functions depend on abundance information. To load just the transcript models, simply just leave out the `counts_file` and `count_cols` arguments to the `add_dataset()` function as shown below.

```python
# for this new example, create a new empty SwanGraph
sg = swan.SwanGraph()
# and add the annotation transcriptome to it
sg.add_annotation(annot_gtf)
```

```python
# add transcriptome datasets from GTF files without
# corresponding abundance information
sg.add_dataset('HepG2_1', hep_1_gtf)
sg.add_dataset('HepG2_2', hep_2_gtf)
sg.add_dataset('HFFc6_1', hff_1_gtf)
sg.add_dataset('HFFc6_2', hff_2_gtf)
sg.add_dataset('HFFc6_3', hff_3_gtf)
```

If you have just added transcript models to the graph via `add_dataset()` and wish to add abundance information, this can be done using the `add_abundance()` function as seen below. Here, the string passed to `count_cols` is the column in the abundance file that corresponds to the dataset, and the argument passed to `dataset_name` is the name of the dataset that has already been added to the SwanGraph in the previous code block.

```python
# add abundance information corresponding to each of the datasets
# we've already added to the SwanGraph
# dataset_name must be a dataset that is already present in the SwanGraph
sg.add_abundance(ab_file, count_cols='hepg2_1', dataset_name='HepG2_1')
sg.add_abundance(ab_file, count_cols='hepg2_2', dataset_name='HepG2_2')
sg.add_abundance(ab_file, count_cols='hffc6_1', dataset_name='HFFc6_1')
sg.add_abundance(ab_file, count_cols='hffc6_2', dataset_name='HFFc6_2')
sg.add_abundance(ab_file, count_cols='hffc6_3', dataset_name='HFFc6_3')
```

## Adding transcript models \(TALON db\) and abundance information

Swan is also directly compatible with TALON databases and can pull transcript models directly from them.

```python
# for this new example, create a new empty SwanGraph
sg = swan.SwanGraph()
# and add the annotation transcriptome to it
sg.add_annotation(annot_gtf)
```

```python
hepg2_whitelist='data/hepg2_whitelist.csv'
hffc6_whitelist='data/hffc6_whitelist.csv'
```

```python
# add datasets directly from a TALON database and abundance
# information from an abundance table
# whitelist option is output from the talon_filter_transcripts
# step, which filters novel isoforms based on their reproducibility
# and for those that exhibit internal priming
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

## Batch adding datasets

If you wish to add multiple datasets to the SwanGraph with a single command, you can use the `add_datasets()` function with a config file. The format of the config file is detailed [here](../faqs/file_formats.md#batch-config-file). You can provide datasets from both a TALON db or a GTF in the config file, as well as the annotation dataset. Below is an example of a config file that contains an annotation to be added to the SwanGraph, as well datasets from GTF files and from a TALON db.

```python
import pandas as pd
config_df = pd.read_csv('config.csv', sep='\t')
config_df
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>fname</th>
      <th>col</th>
      <th>counts_file</th>
      <th>count_cols</th>
      <th>tid_col</th>
      <th>dataset_name</th>
      <th>whitelist</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>data/gencode.v29.annotation.gtf</td>
      <td>annotation</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>data/talon.db</td>
      <td>HepG2_1</td>
      <td>data/all_talon_abundance_filtered.tsv</td>
      <td>hepg2_1</td>
      <td>annot_transcript_id</td>
      <td>hepg2_1</td>
      <td>data/hepg2_whitelist.csv</td>
    </tr>
    <tr>
      <th>2</th>
      <td>data/talon.db</td>
      <td>HepG2_2</td>
      <td>data/all_talon_abundance_filtered.tsv</td>
      <td>hepg2_2</td>
      <td>annot_transcript_id</td>
      <td>hepg2_2</td>
      <td>data/hepg2_whitelist.csv</td>
    </tr>
    <tr>
      <th>3</th>
      <td>data/hffc6_1_talon.gtf</td>
      <td>HFFc6_1</td>
      <td>data/all_talon_abundance_filtered.tsv</td>
      <td>hffc6_1</td>
      <td>annot_transcript_id</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>4</th>
      <td>data/hffc6_2_talon.gtf</td>
      <td>HFFc6_2</td>
      <td>data/all_talon_abundance_filtered.tsv</td>
      <td>hffc6_2</td>
      <td>annot_transcript_id</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>5</th>
      <td>data/hffc6_3_talon.gtf</td>
      <td>HFFc6_3</td>
      <td>data/all_talon_abundance_filtered.tsv</td>
      <td>hffc6_3</td>
      <td>annot_transcript_id</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>

As you can see, the header columns in the config file are equivalent to arguments you would pass into the `add_dataset()` function, and unneccessary arguments based on the input file type (TALON db vs. GTF) can be left blank.


```python
# for this new example, create a new empty SwanGraph
sg = swan.SwanGraph()
```


```python
# add each dataset from the config file with the corresponding input
# settings to the SwanGraph
sg.add_datasets('config.csv')
```

