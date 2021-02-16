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
* recommended: including the gene\_name field will enable you to plot genes from their human-readable names as well
* any non-data header lines must begin with \#
* gene\_ids, gene\_names, and transcript\_ids must be the same across datasets for proper dataset merging 
* exons must be in order under the transcript entry to which they belong

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
* Have a column containing transcript ids that are the same as those loaded via GTF or TALON db
* Have a column containing counts of each transcript for a given dataset column name

Luckily, the names of the column names to obtain transcript ids and counts from are flexible. If you were to add abundance to your SwanGraph with the following line, for instance

```python
sg = swan.SwanGraph('swan.p')
sg.add_abundance('counts_file.tsv', \
    count_cols='counts_dataset', \
    dataset_name='sg_dataset', \
    tid_col='transcript_id')
```

The corresponding abundance file should look something like this:

```text
transcript_id    counts_dataset
ENST00000623083.4    1
ENST00000416931.1    0
ENST00000457540.1    0
ENST00000414273.1    0
ENST00000621981.1    0
ENST00000514057.1    0
ENST00000411249.1    0
ENST00000445118.6    1
ENST00000441765.5    0
```

## TALON db

Swan currently works with TALON databases created with TALON v5.0+

## Batch config file

If you wish to add your datasets at the same time you can use tab-separated configuration file. 

Each column should correspond to a different argument you can pass to the `add_dataset()` function, with the argument name detailed in the header. Columns can be in any order. 

You can provide datasets from both a TALON db or a GTF in the config file, as well as the annotation dataset. 

Arguments that you wish to use the default value for or arguments that are not necessary can be left blank. For instance, if adding from a GTF, there is no need for the `dataset_name` or `whitelist` arguments and these columns can be left blank. 

Below is an example of a config file. As you can see, it adds the annotation transcriptome as well as data from a TALON database and different GTF files.

<div>
<!-- <style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style> -->
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
      <td>hffc6_1</td>
      <td>data/hffc6_whitelist.csv</td>
    </tr>
    <tr>
      <th>4</th>
      <td>data/hffc6_2_talon.gtf</td>
      <td>HFFc6_2</td>
      <td>data/all_talon_abundance_filtered.tsv</td>
      <td>hffc6_2</td>
      <td>annot_transcript_id</td>
      <td>hffc6_2</td>
      <td>data/hffc6_whitelist.csv</td>
    </tr>
    <tr>
      <th>5</th>
      <td>data/hffc6_3_talon.gtf</td>
      <td>HFFc6_3</td>
      <td>data/all_talon_abundance_filtered.tsv</td>
      <td>hffc6_3</td>
      <td>annot_transcript_id</td>
      <td>hffc6_3</td>
      <td>data/hffc6_whitelist.csv</td>
    </tr>
  </tbody>
</table>
</div>

