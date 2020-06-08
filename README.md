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
Swan can be installed directly from PyPi. To install Swan's most recent release, run 

`pip install swan_vis` 

Alternatively, the most recent commits can be installed by git cloning this repo, moving to the swan_vis directory, and running

`pip install .`


After installation with pip, to enable visualizations using dashed edges, run the following command from anywhere in the terminal

`swan_patch_networkx`


## Tutorials
* [Getting started](https://github.com/fairliereese/swan_vis/blob/master/tutorials/Getting%20started.ipynb): how to load data into Swan
* [Visualization tools](https://github.com/fairliereese/swan_vis/blob/master/tutorials/Visualization.ipynb): make gene and transcript-level plots to visualize the complexity of alternative splicing
* [Analysis tools](https://github.com/fairliereese/swan_vis/blob/master/tutorials/Analysis%20tools.ipynb): find differentially expressed genes and transcripts; find isoform-switching genes, discover novel intron retention and exon skipping events


## Wiki links
* [Understanding Swan visualizations](https://github.com/fairliereese/swan_vis/wiki/Understanding-Swan-visualizations)
* [Input file format specifications](https://github.com/fairliereese/swan_vis/wiki/File-format-specifications)


Logo by the wonderful [Eamonn Casey](https://www.instagram.com/designsbyeamonn/)
