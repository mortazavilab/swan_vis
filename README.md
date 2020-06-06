# Swan
<img align="center" width="450" src="figures/swan_logo.png">

Swan is a Python library designed for the analysis and visualization of transcriptomes, especially with long-read transcriptomes in mind. Users can merge transcriptomes from different datasets and explore transcript models distict splicing and expression patterns across datasets.
<!-- 
## Table of contents
* [Installation](#installation)
* [Swan]
* [Vignette](#vignette)
* [Input File Formats](#fileformats) -->

## What can Swan do?
Swan can make informative plots, find differentially expressed genes and transcripts, find isoform-switching genes, and discover novel exon skipping and intron retention events.

## <a name="installation"></a>Installation
Swan can be installed directly from PyPi. Run `pip install swan_vis` to install Swan's most recent release. Alternatively, the most recent commits can be installed by git cloning this repo and running `pip install .` in the downloaded directory.

After installation, to enable visualizations using dashed edges, run `swan_patch_networkx` from anywhere in the terminal.

## Tutorials
* [Getting started](https://github.com/fairliereese/swan_vis/blob/master/tutorials/Getting%20started.ipynb): how to load data into Swan
* [Visualization tools](https://github.com/fairliereese/swan_vis/blob/master/tutorials/Visualization.ipynb): make gene and transcript-level plots to visualize the complexity of alternative splicing
* [Analysis tools](https://github.com/fairliereese/swan_vis/blob/master/tutorials/Analysis%20tools.ipynb): find differentially expressed genes and transcripts; find isoform-switching genes, discover novel intron retention and exon skipping events




<!-- # <a name="vignette"></a>Vignette
Here is an example of how to run Swan. 

```python
# import Swan library
import swan_vis as swan

# annotation and dataset transcriptome gtfs
annot_gtf = 'example/annot_adrm1.gtf'
hepg2_gtf = 'example/hepg2_adrm1.gtf'
hffc6_gtf = 'example/hffc6_adrm1.gtf'

# abundance table with counts for each transcript corresponding to datasets
ab_file = 'example/abundance.tsv'

# initialize an empty SwanGraph object
sg = swan.SwanGraph()

# add the annotation to the SwanGraph
sg.add_annotation(annot_gtf)

# add HepG2 and HFFc6 datasets to the SwanGraph along with the corresponding
# abundance information
sg.add_dataset('HepG2', hepg2_gtf,
	counts_file=ab_file,
	counts_cols='hepg2')
sg.add_dataset('HFFc6', hffc6_gtf,
	counts_file=ab_file,
	counts_cols='hffc6')
```

# <a name="fileformats"></a>Input File Formats
Currently, Swan works with several file types:
* GTF, used to add transcript models
	* Must have transcript entries. Future releases will only require exon entries.
	* Last column must have gene_id, gene_name, and transcript_id entries. These entries must correspond across datasets to be properly compared.
	* Swan includes a GTF validator that can be used to see if all the correct entires and fields are present in an input GTF it can be run from python as follows
```python
import swan_vis as swan
swan.validate_gtf('example/annot.gtf')
```
* TALON database, used to add transcript models
	* Currently, Swan has been tested with TALON databases from v5.0+.
* Abundance matrix
	* Each row in an abundance matrix corresponds to one transcript's expression. Each row must contain an entry corresponding to a transcript id from a dataset or annotation that has already been added to the SwanGraph.
Please also see the files in the [example](https://github.com/fairliereese/swan_vis/tree/master/example) folder.
 -->

