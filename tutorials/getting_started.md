#  Getting started: initializing, adding data, and saving your SwanGraph

First, if you haven't already, make sure to [install Swan](https://github.com/fairliereese/swan_vis/wiki#installation).
After installing, you'll be able to run Swan from Python.

Then, download the data and the reference transcriptome annotation from [here](http://crick.bio.uci.edu/freese/swan_files_example/). The bash commands to do so are given below.

The main workflow to get started with Swan consists of:

0. [Starting up Swan and initializing your SwanGraph](getting_started.md#starting-up-swan-and-initializing-your-swangraph)
1. [Adding a reference transcriptome (optional)](getting_started.md#adding-a-reference-transcriptome)
2. Adding a transcriptome for your samples
    * [From a GTF](getting_started.md#adding-transcript-models-from-a-gtf)
    * [From a TALON db](getting_started.md#adding-transcript-models-from-a-talon-db)
3. [Adding datasets and their expression values](getting_started.md#adding-datasets-and-their-abundance)
4. [Adding metadata to your datasets](getting_started.md#adding-metadata)

Other sections:
* [Example data download](getting_started.md#download-example-data)
* [Saving and loading your SwanGraph](getting_started.md#saving-and-loading-your-swangraph)

This page can also be read from top to bottom, just know that you may be running things more than once!

For information on the file formats needed to use Swan, please read the [file format specifications FAQ](https://freese.gitbook.io/swan/faqs/file_formats).

<!-- Running this tutorial (with only one of the dataset addition options) on my laptop took around 7 minutes and 5 GB of RAM.  -->

## <a name="data_download"></a> Download example data

This data is the data used in the [Swan publication](https://academic.oup.com/bioinformatics/article/37/9/1322/5912931)

Run this block in your bash terminal
```bash
mkdir data
mkdir figures
cd data/

# download files
wget http://crick.bio.uci.edu/freese/swan_files.tgz

# expand files
tar -xzf swan_files.tgz
mv swan_files/* .
rm -r swan_files/

# download reference annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
gunzip gencode.v29.annotation.gtf.gz

cd ../
```

Alternatively, just run on a smaller example, chr20.

Run this block in your bash terminal

```bash
mkdir data
mkdir figures
cd data/

# download files
wget http://crick.bio.uci.edu/freese/swan_files_example.tgz

# expand files
tar -xzf swan_files_example.tar.gz
mv swan_files_example/* .
rm -r swan_files_example/

cd ../
```

## <a name="init"></a>Starting up Swan and initializing your SwanGraph

The rest of the code in this tutorial should be run in using Python

Initialize an empty SwanGraph and add the transcriptome annotation to the SwanGraph.


```python
import swan_vis as swan

# initialize a new SwanGraph
sg = swan.SwanGraph()
```

**Note:** to initialize a SwanGraph in single-cell mode (which will avoid calculating percent isoform use \[pi\] numbers for each cell), use the following code:

```python
sg = swan.SwanGraph(sc=True)
```


```python
annot_gtf = 'data/gencode.v29.annotation.gtf'
data_gtf = 'data/all_talon_observedOnly.gtf'
ab_file = 'data/all_talon_abundance_filtered.tsv'
talon_db = 'data/talon.db'
pass_list = 'data/all_pass_list.csv'
meta = 'data/metadata.tsv'
```

## <a name="add_trans"></a>Adding a reference transcriptome


```python
# add an annotation transcriptome
sg.add_annotation(annot_gtf)
```


    Adding annotation to the SwanGraph


## <a name="add_gtf"></a>Adding transcript models from a GTF

Add all filtered transcript models to the SwanGraph.


```python
# add a dataset's transcriptome and abundance information to
# the SwanGraph
sg.add_transcriptome(data_gtf)
```


    Adding transcriptome to the SwanGraph


## <a name="add_ab"></a>Adding datasets and their abundance

Use an abundance matrix with columns for each desired dataset to add datasets to the SwanGraph.


```python
# add each dataset's abundance information to the SwanGraph
sg.add_abundance(ab_file)
```


    Adding abundance for datasets hepg2_1, hepg2_2, hffc6_1, hffc6_2, hffc6_3 to SwanGraph.

##  <a name="save_load"></a>Saving and loading your SwanGraph

Following this, you can save your SwanGraph so you can easily work with it again without re-adding all the data.


```python
# save the SwanGraph as a Python pickle file
sg.save_graph('swan')
```

    Saving graph as swan.p


And you can reload the graph again.


```python
# load up a saved SwanGraph from a pickle file
sg = swan.read('swan.p')
```

    Read in graph from swan.p


##  <a name="add_db"></a>Adding transcript models from a TALON DB

Swan is also directly compatible with TALON databases and can pull transcript models directly from them. You can also optionally pass in a list of isoforms from [`talon_filter_transcripts`](https://github.com/mortazavilab/TALON#talon_filter) to filter your input transcript models.


```python
# for this new example, create a new empty SwanGraph
sg = swan.SwanGraph()

# and add the annotation transcriptome to it
sg.add_annotation(annot_gtf)

# add transcriptome from TALON db
sg.add_transcriptome(talon_db, pass_list=pass_list)

# add each dataset's abundance information to the SwanGraph
sg.add_abundance(ab_file)
```


    Adding annotation to the SwanGraph

    Adding transcriptome to the SwanGraph

    Adding abundance for datasets hepg2_1, hepg2_2, hffc6_1, hffc6_2, hffc6_3 to SwanGraph.

##  <a name="add_meta"></a>Adding metadata

Swan provides functionality to perform tests and plotting on the basis of metadata categories.


```python
sg.add_metadata(meta)
```

```python
sg.adata.obs
```

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>dataset</th>
      <th>cell_line</th>
      <th>replicate</th>
    </tr>
    <tr>
      <th>index</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>hepg2_1</th>
      <td>hepg2_1</td>
      <td>hepg2</td>
      <td>1</td>
    </tr>
    <tr>
      <th>hepg2_2</th>
      <td>hepg2_2</td>
      <td>hepg2</td>
      <td>2</td>
    </tr>
    <tr>
      <th>hffc6_1</th>
      <td>hffc6_1</td>
      <td>hffc6</td>
      <td>1</td>
    </tr>
    <tr>
      <th>hffc6_2</th>
      <td>hffc6_2</td>
      <td>hffc6</td>
      <td>2</td>
    </tr>
    <tr>
      <th>hffc6_3</th>
      <td>hffc6_3</td>
      <td>hffc6</td>
      <td>3</td>
    </tr>
  </tbody>
</table>
</div>
