# Swan

![](.gitbook/assets/swan_logo.png)

Swan is a Python library designed for the analysis and visualization of transcriptomes, especially with long-read transcriptomes in mind. Users can add transcriptomes from different datasets and explore distinct splicing and expression patterns across datasets.

Please visit the [Swan repository](https://github.com/mortazavilab/swan_vis) to download and view the source code

Also see the [Swan manuscript repository](https://github.com/fairliereese/swan_paper) for the exact commands used to do the analysis in our [publication](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa836/5912931).

Also see our [website](https://freese.gitbook.io/swan/) for in-depth tutorials and documentation

## What can Swan do?

Swan can make informative plots, find differentially expressed genes and transcripts, find isoform-switching genes, and discover novel exon skipping and intron retention events.

## Installation

Swan can be installed directly from PyPi. To install Swan's most recent release, run

`pip install swan_vis`

Alternatively, the most recent commits can be installed by git cloning [the Swan repository](https://github.com/fairliereese/swan_vis), moving to the swan\_vis directory, and running

`pip install .`

<!-- After installation with pip, to enable visualizations using dashed edges, run the following command from anywhere in the terminal

`swan_patch_networkx` -->

## Tutorials

* [Data preprocessing with TALON](tutorials/data_processing.md)
* [Getting started](tutorials/getting_started.md): how to load data into Swan
* [Visualization tools](tutorials/visualization.md): make gene and transcript-level plots to visualize the complexity of alternative splicing
* [Analysis tools](tutorials/analysis_tools.md): find differentially expressed genes and transcripts; find isoform-switching genes, discover novel intron retention and exon skipping events
* [Scanpy compatibility](tutorials/scanpy_compatibility.md): Some brief examples of how to use external Scanpy plotting functions on Swan objects

## FAQs

* [Understanding Swan visualizations](faqs/understanding_swan_vis.md)
* [Additional utilities](faqs/utilities.md)
* [SwanGraph data structure](faqs/data_structure.md)
* [Input file format specifications](faqs/file_formats.md)

## Comprehensive Documentation

For full documentaion, please visit [our website](https://freese.gitbook.io/swan/)

Logo by the wonderful [Eamonn Casey](https://www.instagram.com/designsbyeamonn/)
