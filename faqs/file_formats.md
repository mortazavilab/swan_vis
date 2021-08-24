# File format specifications

File formats in bioinformatics are notoriously hard to standardize. We hope that this documentation provides the user with a clear idea of what is need as input into Swan.

## Table of contents

* [GTF](file_formats.md#gtf)
* [Abundance matrix](file_formats.md#abundance-matrix)
* [TALON db](file_formats.md#talon-db)
* [Batch add datasets config file](file_formats.md#batch-config-file)

## GTF

In Swan, transcript models are loaded from GTFs. To work with Swan, GTFs must adhere to the following specifications:

* transcript and exon entries in column 3 - this is a dependency we would like to remove in the future but for now this is the way it works
* gene\_id and transcript\_id attributes \(for transcripts and exons\) in column 9.
* recommended: including the transcript\_name and gene\_name field will enable you to plot genes and transcript with their human-readable names as well
* any non-data header lines must begin with \#

Here is an example of what the first few lines of a GTF should look like:

```text
##description: evidence-based annotation of the human genome (GRCh38), version 29 (Ensembl 94)
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2018-08-30
chr1    HAVANA    gene    11869    14409    .    +    .    gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";
chr1    HAVANA    transcript    11869    14409    .    +    .    gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; level 2; transcript_support_level "1"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1    HAVANA    exon    11869    12227    .    +    .    gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 1; exon_id "ENSE00002234944.1"; level 2; transcript_support_level "1"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
```

If you are having trouble with your GTF, Swan includes a quick GTF validator which can tell you if your file seems to have an unconventional header or lacks entries needed to run Swan. It cannot tell you if your gene/transcript names/ids match across datasets, or if your exon entries are in the correct order after the corresponding transcript entry. The validator can be run as follows:

```python
import swan_vis as swan
swan.validate_gtf('test.gtf')
```

## Abundance matrix

Swan can load abundance information for more meaningful analysis and visualizations. To work with Swan, abundance matrices must:

* Be tab-separated
* First column are transcript ids that are the same as those loaded via GTF or TALON db
* Columns labelled by their dataset names containing raw counts for each transcript
* Alternatively, a [TALON abundance file](https://github.com/mortazavilab/TALON#talon_abundance) can be used in its unaltered form

Sample abundance file:

| transcript_id | dataset1 | dataset2 |
| ---------- | ----------- | ----------- |
| ENST00000416931.1 | 0 | 1 |
| ENST00000414273.1 | 0 | 2 |
| ENST00000621981.1 | 0 | 0 |
| ENST00000514057.1 | 0 | 1 |
| ENST00000411249.1 | 0 | 0 |
| ENST00000445118.6 | 1 | 0 |
| ENST00000441765.5 | 0 | 0 |


## TALON db

Swan currently works with TALON databases created with TALON v5.0+

## Metadata file

Metadata files must:

* Contain a column labeled `dataset` whose entries correspond to the datasets from an already-added abundance file
* Be tab-separated

Sample metadata file (corresponds to above abundance file):

| dataset | sex | tissue |
| ------- | --- | ------ |
| dataset1 | M | heart |
| dataset2 | F | liver |
