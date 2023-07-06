#  Getting started: initializing, adding data, and saving your SwanGraph

First, if you haven't already, make sure to [install Swan](https://github.com/fairliereese/swan_vis/wiki#installation).
After installing, you'll be able to run Swan from Python.

Then, download the data and the reference transcriptome annotation from [here](http://crick.bio.uci.edu/freese/swan_files_example/). The bash commands to do so are given below.

The main workflow to get started with Swan consists of:
1. [Adding a reference transcriptome](#add_trans)
2. Adding a transcriptome for your samples
    * [From a GTF](#add_gtf)
    * [From a TALON db](#add_db)
3. [Adding abundance information](#add_ab)
    * [From a TSV](#add_ab_tsv)
    * [From an AnnData](#add_ab_ad)
    * [On the gene level](#add_ab_gene)
4. [Adding metadata to your datasets](#add_meta)

Other sections:
* [Example data download](#data_download)
* [Starting and initializing your SwanGraph](#init)
* [Saving and loading your SwanGraph](#save_load)
* [Behavior with Cerberus](#cerberus)

This page can also be read from top to bottom, just know that you may be running things more than once!

For information on the file formats needed to use Swan, please read the [file format specifications FAQ](https://freese.gitbook.io/swan/faqs/file_formats).

<!-- Running this tutorial (with only one of the dataset addition options) on my laptop took around 7 minutes and 5 GB of RAM.  -->

## <a name="data_download"></a> Download example data

This data is the data used in the [Swan publication](https://academic.oup.com/bioinformatics/article/37/9/1322/5912931)

Run this block in your bash terminal
```bash
mkdir figures

# download files
wget https://zenodo.org/record/8118614/files/data.tgz

# expand files
tar -xzf data.tgz
```

<!-- Alternatively, just run on a smaller example, chr20.

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
``` -->

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
adata_file = 'data/swan_anndata.h5ad'
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
# add a dataset's transcriptome to the SwanGraph
sg.add_transcriptome(data_gtf)
```


    Adding transcriptome to the SwanGraph


## <a name="add_ab"></a>Adding abundance information

### <a name="add_ab_tsv"></a>Adding abundance from a TSV

You can use an abundance matrix with columns for each desired dataset to add datasets to the SwanGraph. The file format is specified [here](https://freese.gitbook.io/swan/faqs/file_formats#abundance-matrix).


```python
# add each dataset's abundance information to the SwanGraph
sg.add_abundance(ab_file)
```


    Adding abundance for datasets hepg2_1, hepg2_2, hffc6_1, hffc6_2, hffc6_3 to SwanGraph.


    /Users/fairliereese/miniconda3/lib/python3.7/site-packages/anndata/_core/anndata.py:120: ImplicitModificationWarning: Transforming to str index.
      warnings.warn("Transforming to str index.", ImplicitModificationWarning)


### <a name="add_ab_ad"></a>Adding abundance from an AnnData

If you have abundance information and metadata information in AnnData format, you can use this as direct input into Swan. This will help circumvent the dense matrix representation of the TSV in the case of very large datasets or single-cell data.


```python
# add abundance for each dataset from the AnnData into the SwanGraph
sg = swan.SwanGraph()
sg.add_annotation(annot_gtf)
sg.add_transcriptome(data_gtf)
sg.add_adata(adata_file)
```


    Adding annotation to the SwanGraph

    Adding transcriptome to the SwanGraph

    Adding abundance for datasets hepg2_1, hepg2_2, hffc6_1, hffc6_2, hffc6_3 to SwanGraph.
    Calculating TPM...
    Calculating PI...
    Calculating edge usage...


    /Users/fairliereese/Documents/programming/mortazavi_lab/bin/swan_vis/swan_vis/swangraph.py:828: FutureWarning: X.dtype being converted to np.float32 from float64. In the next version of anndata (0.9) conversion will not be automatic. Pass dtype explicitly to avoid this warning. Pass `AnnData(X, dtype=X.dtype, ...)` to get the future behavour.
      adata = anndata.AnnData(var=var, obs=obs, X=X)
    /Users/fairliereese/miniconda3/envs/scanpy_2/lib/python3.7/site-packages/anndata/_core/anndata.py:121: ImplicitModificationWarning: Transforming to str index.
      warnings.warn("Transforming to str index.", ImplicitModificationWarning)


    Calculating TSS usage...


    /Users/fairliereese/Documents/programming/mortazavi_lab/bin/swan_vis/swan_vis/swangraph.py:759: FutureWarning: X.dtype being converted to np.float32 from float64. In the next version of anndata (0.9) conversion will not be automatic. Pass dtype explicitly to avoid this warning. Pass `AnnData(X, dtype=X.dtype, ...)` to get the future behavour.
      adata = anndata.AnnData(var=var, obs=obs, X=X)


    Calculating TES usage...


By adding abundance information from either an AnnData or TSV file, Swan will also automatically calculate the counts and TPM for each TSS, TES, and intron or exon. If you had previously used `add_transcriptome()` to add a GTF that was generated by [Cerberus](https://github.com/mortazavilab/cerberus/tree/master) or uses Cerberus-style transcript IDs (ie. \<gene_id\>\[1,1,1\]), Swan will also calculate intron chain counts and TPM automatically.

### <a name="add_ab_gene"></a>Adding gene-level abundance

You can also store gene expression in the SwanGraph. This can either be done from a TALON abundance TSV that contains transcript-level counts where the counts for each transcript will be summed up across the gene. Alternatively, supply this function a gene-level counts matrix where the first column is the gene ID rather than the transcript ID, but otherwise follows the [input abundance TSV format](https://freese.gitbook.io/swan/faqs/file_formats#abundance-matrix).


```python
# add gene-level abundance to the SwanGraph
sg.add_abundance(ab_file, how='gene')
```

    /Users/fairliereese/Documents/programming/mortazavi_lab/bin/swan_vis/swan_vis/swangraph.py:363: FutureWarning: X.dtype being converted to np.float32 from int64. In the next version of anndata (0.9) conversion will not be automatic. Pass dtype explicitly to avoid this warning. Pass `AnnData(X, dtype=X.dtype, ...)` to get the future behavour.
      adata = anndata.AnnData(var=var, obs=obs, X=X)



    Adding abundance for datasets hepg2_1, hepg2_2, hffc6_1, hffc6_2, hffc6_3 to SwanGraph.
    Calculating TPM...


##  <a name="save_load"></a>Saving and loading your SwanGraph

Following this, you can save your SwanGraph so you can easily work with it again without re-adding all the data.


```python
# save the SwanGraph as a Python pickle file
sg.save_graph('data/swan')
```

    Saving graph as data/swan.p


And you can reload the graph again.


```python
# load up a saved SwanGraph from a pickle file
sg = swan.read('data/swan.p')
```

    Read in graph from data/swan.p


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


    /Users/fairliereese/Documents/programming/mortazavi_lab/bin/swan_vis/swan_vis/swangraph.py:346: FutureWarning: X.dtype being converted to np.float32 from int64. In the next version of anndata (0.9) conversion will not be automatic. Pass dtype explicitly to avoid this warning. Pass `AnnData(X, dtype=X.dtype, ...)` to get the future behavour.
      adata = anndata.AnnData(var=var, obs=obs, X=X)



    Adding abundance for datasets hepg2_1, hepg2_2, hffc6_1, hffc6_2, hffc6_3 to SwanGraph.
    Calculating TPM...
    Calculating PI...
    Calculating edge usage...


    /Users/fairliereese/Documents/programming/mortazavi_lab/bin/swan_vis/swan_vis/swangraph.py:810: FutureWarning: X.dtype being converted to np.float32 from float64. In the next version of anndata (0.9) conversion will not be automatic. Pass dtype explicitly to avoid this warning. Pass `AnnData(X, dtype=X.dtype, ...)` to get the future behavour.

    /Users/fairliereese/miniconda3/envs/scanpy_2/lib/python3.7/site-packages/anndata/_core/anndata.py:121: ImplicitModificationWarning: Transforming to str index.
      warnings.warn("Transforming to str index.", ImplicitModificationWarning)


    Calculating TSS usage...


    /Users/fairliereese/Documents/programming/mortazavi_lab/bin/swan_vis/swan_vis/swangraph.py:741: FutureWarning: X.dtype being converted to np.float32 from float64. In the next version of anndata (0.9) conversion will not be automatic. Pass dtype explicitly to avoid this warning. Pass `AnnData(X, dtype=X.dtype, ...)` to get the future behavour.



    Calculating TES usage...


##  <a name="add_meta"></a>Adding metadata

Swan provides functionality to perform tests and plotting on the basis of metadata categories. Add metadata by calling the `SwanGraph.add_metadata()` function, or use the `SwanGraph.add_adata()` function to add both expression information and metadata at the same time.


```python
sg = swan.read('data/swan.p')
sg.add_metadata(meta)
```

    Read in graph from data/swan.p


    /Users/fairliereese/miniconda3/envs/scanpy_2/lib/python3.7/site-packages/anndata/_core/anndata.py:798: UserWarning:
    AnnData expects .obs.index to contain strings, but got values like:
        [0, 1, 2, 3, 4]

        Inferred to be: integer

      value_idx = self._prep_dim_index(value.index, attr)



```python
sg.adata.obs
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
      <th>cell_line</th>
      <th>replicate</th>
      <th>dataset</th>
      <th>total_counts</th>
      <th>description</th>
    </tr>
    <tr>
      <th>index</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>hepg2_1</th>
      <td>hepg2</td>
      <td>1</td>
      <td>hepg2_1</td>
      <td>499647.0</td>
      <td>liver</td>
    </tr>
    <tr>
      <th>hepg2_2</th>
      <td>hepg2</td>
      <td>2</td>
      <td>hepg2_2</td>
      <td>848447.0</td>
      <td>liver</td>
    </tr>
    <tr>
      <th>hffc6_1</th>
      <td>hffc6</td>
      <td>1</td>
      <td>hffc6_1</td>
      <td>761493.0</td>
      <td>fibroblast</td>
    </tr>
    <tr>
      <th>hffc6_2</th>
      <td>hffc6</td>
      <td>2</td>
      <td>hffc6_2</td>
      <td>787967.0</td>
      <td>fibroblast</td>
    </tr>
    <tr>
      <th>hffc6_3</th>
      <td>hffc6</td>
      <td>3</td>
      <td>hffc6_3</td>
      <td>614921.0</td>
      <td>fibroblast</td>
    </tr>
  </tbody>
</table>
</div>



## <a name="cerberus"></a> Behavior with Cerberus

When you use a [Cerberus](https://github.com/mortazavilab/cerberus/tree/master) GTF in `SwanGraph.add_annotation()` or `SwanGraph.add_transcriptome()`, keep in mind the following:

* Swan will use the TSS / TES assignments as dictated by Cerberus to define unique entries in `SwanGraph.tss_adata` and `SwanGraph.tes_adata`. For instance, if the same vertex is used in more than one gene, they will still be treated as separate vertices in the TSS / TES AnnDatas.
* Swan will automatically pull intron chain information from the transcript triplet in Cerberus and use it to generate an AnnData tracking the expression of intron chains separately from the transcripts they come from in `SwanGraph.ic_adata`. This can also be used to perform isoform switching tests.
* Currently, Swan does not parse Cerberus novelty categories. We are hoping to support this in a future release.


```python
sg = swan.read('data/swan_modelad.p')
sg.ic_adata.var.tail()
```

    Read in graph from data/swan_modelad.p





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
      <th>gid</th>
      <th>gname</th>
      <th>ic_name</th>
      <th>n_cells</th>
    </tr>
    <tr>
      <th>ic_id</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENSMUSG00000118369_2</th>
      <td>ENSMUSG00000118369</td>
      <td>Gm30541</td>
      <td>Gm30541_2</td>
      <td>14</td>
    </tr>
    <tr>
      <th>ENSMUSG00000118380_3</th>
      <td>ENSMUSG00000118380</td>
      <td>Gm36037</td>
      <td>Gm36037_3</td>
      <td>1</td>
    </tr>
    <tr>
      <th>ENSMUSG00000118382_1</th>
      <td>ENSMUSG00000118382</td>
      <td>Gm8373</td>
      <td>Gm8373_1</td>
      <td>2</td>
    </tr>
    <tr>
      <th>ENSMUSG00000118383_1</th>
      <td>ENSMUSG00000118383</td>
      <td>Gm50321</td>
      <td>Gm50321_1</td>
      <td>14</td>
    </tr>
    <tr>
      <th>ENSMUSG00000118390_1</th>
      <td>ENSMUSG00000118390</td>
      <td>Gm50102</td>
      <td>Gm50102_1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>
