# Analysis

Swan has several analysis options to use.

Table of contents

* [Differential gene expression](analysis_tools.md#differential-gene-expression-tests)
* [Differential transcript expression](analysis_tools.md#differential-transcript-expression-tests)
* [Isoform switching](analysis_tools.md#isoform-switching)
* [Exon skipping and intron retention](analysis_tools.md#exon-skipping-and-intron-retention)

```python
import swan_vis as swan

sg = swan.SwanGraph('data/swan.p')
```

```text
Graph from data/swan.p loaded
```

## Differential gene expression tests

Differential gene expression testing in Swan is implemented via [diffxpy](https://github.com/theislab/diffxpy). To run the test, first partition the datasets that you have added to your SwanGraph into biological replicates. Then, use this grouping to run the differential expression test.

```python
dataset_groups = [['HepG2_1','HepG2_2'],['HFFc6_1','HFFc6_2','HFFc6_3']]
sg.de_gene_test(dataset_groups)
```

|  | gid | pval | qval | log2fc | mean | zero\_mean | grad | coef\_mle | coef\_sd | ll | gname |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 93003 | ENSG00000158874.11 | 0.0 | 0.0 | -9.798466 | 7203.040430 | False | 8.236582e-07 | -9.798466 | 0.597998 | -21.964623 | APOA2 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 99184 | ENSG00000163631.16 | 0.0 | 0.0 | -9.657541 | 6256.314453 | False | 2.962665e-07 | -9.657541 | 0.583054 | -20.338973 | ALB |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 79865 | ENSG00000145192.12 | 0.0 | 0.0 | -9.320415 | 4466.045703 | False | 4.583478e-08 | -9.320415 | 0.591606 | -20.613694 | AHSG |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 4810 | ENSG00000026025.15 | 0.0 | 0.0 | 9.226718 | 9419.225325 | False | 1.069652e-06 | 9.226718 | 0.569277 | -24.014993 | VIM |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 139314 | ENSG00000197249.13 | 0.0 | 0.0 | -9.106863 | 5051.564227 | False | 1.033368e-06 | -9.106863 | 0.496072 | -20.785889 | SERPINA1 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 176368 | ENSG00000242968.1 | 1.0 | 1.0 | 0.000000 | 1.000000 | False | 0.000000e+00 | 0.000000 | 0.912871 | 0.000000 | AC096992.1 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 176367 | ENSG00000242963.1 | 1.0 | 1.0 | 0.000000 | 1.000000 | False | 0.000000e+00 | 0.000000 | 0.912871 | 0.000000 | AC026336.1 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 176366 | ENSG00000242960.1 | 1.0 | 1.0 | 0.000000 | 1.000000 | False | 0.000000e+00 | 0.000000 | 0.912871 | 0.000000 | FTH1P23 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 176365 | ENSG00000242958.1 | 1.0 | 1.0 | 0.000000 | 1.000000 | False | 0.000000e+00 | 0.000000 | 0.912871 | 0.000000 | AC040975.1 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 163669 | ENSG00000230383.1 | 1.0 | 1.0 | 0.000000 | 1.000000 | False | 0.000000e+00 | 0.000000 | 0.912871 | 0.000000 | AC009245.1 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


The results of this test are stored in `sg.deg_test` so they can be accessed later as follows. Test results will also be stored if `save_graph()` is run again so the user can easily load the results up.

```python
sg.deg_test.head()
```

|  | gid | pval | qval | log2fc | mean | zero\_mean | grad | coef\_mle | coef\_sd | ll | gname |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 93003 | ENSG00000158874.11 | 0.0 | 0.0 | -9.798466 | 7203.040430 | False | 8.236582e-07 | -9.798466 | 0.597998 | -21.964623 | APOA2 |
| 99184 | ENSG00000163631.16 | 0.0 | 0.0 | -9.657541 | 6256.314453 | False | 2.962665e-07 | -9.657541 | 0.583054 | -20.338973 | ALB |
| 79865 | ENSG00000145192.12 | 0.0 | 0.0 | -9.320415 | 4466.045703 | False | 4.583478e-08 | -9.320415 | 0.591606 | -20.613694 | AHSG |
| 4810 | ENSG00000026025.15 | 0.0 | 0.0 | 9.226718 | 9419.225325 | False | 1.069652e-06 | 9.226718 | 0.569277 | -24.014993 | VIM |
| 139314 | ENSG00000197249.13 | 0.0 | 0.0 | -9.106863 | 5051.564227 | False | 1.033368e-06 | -9.106863 | 0.496072 | -20.785889 | SERPINA1 |

Swan can also automatically subset the test summary table to pull out genes that pass a certain significance threshold. These genes can be directly passed into Swan's gene plotting functions, `gen_report()` or `plot_graph()`

```python
gene_ids, gene_summary = sg.get_de_genes(q=0.05)
```

```python
print(gene_ids[:5])
gene_summary.head()
```

```text
['APOA2', 'ALB', 'AHSG', 'VIM', 'SERPINA1']
```

|  | gid | pval | qval | log2fc | mean | zero\_mean | grad | coef\_mle | coef\_sd | ll | gname |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 93003 | ENSG00000158874.11 | 0.0 | 0.0 | -9.798466 | 7203.040430 | False | 8.236582e-07 | -9.798466 | 0.597998 | -21.964623 | APOA2 |
| 99184 | ENSG00000163631.16 | 0.0 | 0.0 | -9.657541 | 6256.314453 | False | 2.962665e-07 | -9.657541 | 0.583054 | -20.338973 | ALB |
| 79865 | ENSG00000145192.12 | 0.0 | 0.0 | -9.320415 | 4466.045703 | False | 4.583478e-08 | -9.320415 | 0.591606 | -20.613694 | AHSG |
| 4810 | ENSG00000026025.15 | 0.0 | 0.0 | 9.226718 | 9419.225325 | False | 1.069652e-06 | 9.226718 | 0.569277 | -24.014993 | VIM |
| 139314 | ENSG00000197249.13 | 0.0 | 0.0 | -9.106863 | 5051.564227 | False | 1.033368e-06 | -9.106863 | 0.496072 | -20.785889 | SERPINA1 |

## Differential transcript expression tests

Similarly, Swan can run tests to find differentially expressed transcript isoforms. The input and output to these functions are identical to that of the differential gene tests.

```python
dataset_groups = [['HepG2_1','HepG2_2'],['HFFc6_1','HFFc6_2','HFFc6_3']]
sg.de_transcript_test(dataset_groups);
sg.det_test.head()
```

|  | tid | pval | qval | log2fc | mean | zero\_mean | grad | coef\_mle | coef\_sd | ll | gid | gname |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 21203 | ENST00000367990.7 | 0.0 | 0.0 | -9.641046 | 6153.972461 | False | 2.680981e-07 | -9.641046 | 0.599458 | -21.723931 | ENSG00000158874.11 | APOA2 |
| 6792 | ENST00000295897.8 | 0.0 | 0.0 | -9.599193 | 5901.748828 | False | 7.406762e-07 | -9.599193 | 0.583689 | -20.330176 | ENSG00000163631.16 | ALB |
| 31598 | ENST00000393087.8 | 0.0 | 0.0 | -9.435849 | 5012.441016 | False | 6.445624e-07 | -9.435849 | 0.584451 | -20.120285 | ENSG00000197249.13 | SERPINA1 |
| 39289 | ENST00000411641.6 | 0.0 | 0.0 | -9.236873 | 4108.149609 | False | 4.956935e-08 | -9.236873 | 0.592742 | -20.527872 | ENSG00000145192.12 | AHSG |
| 1044 | ENST00000224237.9 | 0.0 | 0.0 | 9.203737 | 9205.245637 | False | 1.456279e-06 | 9.203737 | 0.569203 | -23.507290 | ENSG00000026025.15 | VIM |

And Swan can subset the results for you based on a q-value significance threshold. The resultant transcript ids can then be passed into Swan's transcript plotting function, `plot_transcript_path()`.

```python
transcript_ids, transcript_summary = sg.get_de_transcripts(q=0.05)
```

```python
print(transcript_ids[:5])
transcript_summary.head()
```

```text
['ENST00000367990.7', 'ENST00000295897.8', 'ENST00000393087.8', 'ENST00000411641.6', 'ENST00000224237.9']
```

|  | tid | pval | qval | log2fc | mean | zero\_mean | grad | coef\_mle | coef\_sd | ll | gid | gname |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 21203 | ENST00000367990.7 | 0.0 | 0.0 | -9.641046 | 6153.972461 | False | 2.680981e-07 | -9.641046 | 0.599458 | -21.723931 | ENSG00000158874.11 | APOA2 |
| 6792 | ENST00000295897.8 | 0.0 | 0.0 | -9.599193 | 5901.748828 | False | 7.406762e-07 | -9.599193 | 0.583689 | -20.330176 | ENSG00000163631.16 | ALB |
| 31598 | ENST00000393087.8 | 0.0 | 0.0 | -9.435849 | 5012.441016 | False | 6.445624e-07 | -9.435849 | 0.584451 | -20.120285 | ENSG00000197249.13 | SERPINA1 |
| 39289 | ENST00000411641.6 | 0.0 | 0.0 | -9.236873 | 4108.149609 | False | 4.956935e-08 | -9.236873 | 0.592742 | -20.527872 | ENSG00000145192.12 | AHSG |
| 1044 | ENST00000224237.9 | 0.0 | 0.0 | 9.203737 | 9205.245637 | False | 1.456279e-06 | 9.203737 | 0.569203 | -23.507290 | ENSG00000026025.15 | VIM |

## Isoform switching

We wanted to include a module to conduct rudimentary isoform switching analysis as well. We define a gene that exhibits isoform switching for our purposes as a gene that is not differentially expressed that has transcript isoforms that are differentially expressed. We have provided code to detect such instances. To run it, `de_gene_test()` and `de_transcript_test()` must first be run.

```python
is_genes, is_table = sg.find_isoform_switching_genes(q=0.05)
```

```python
print(is_genes[:5])
is_table.head()
```

```text
['ENSG00000197746.13', 'ENSG00000067225.17', 'ENSG00000117450.13', 'ENSG00000105254.11', 'ENSG00000177600.8']
```

|  | tid | pval | qval | log2fc | mean | zero\_mean | grad | coef\_mle | coef\_sd | ll | gid | gname |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 0 | TALONT000283514 | 5.099119e-10 | 8.195811e-08 | 5.538936 | 153.044266 | False | 2.116293e-07 | 5.538936 | 0.891074 | -21.442096 | ENSG00000197746.13 | PSAP |
| 1 | ENST00000394936.7 | 9.324034e-05 | 5.520899e-03 | -0.970022 | 278.640723 | False | 5.421814e-08 | -0.970022 | 0.248244 | -28.092922 | ENSG00000197746.13 | PSAP |
| 2 | TALONT000316712 | 2.490252e-09 | 3.622447e-07 | 5.038395 | 92.933382 | False | 2.948318e-08 | 5.038395 | 0.845071 | -19.486628 | ENSG00000067225.17 | PKM |
| 3 | TALONT000375121 | 2.220446e-16 | 7.353454e-14 | -4.855414 | 51.973529 | False | 1.056665e-08 | -4.855414 | 0.592949 | -11.909066 | ENSG00000117450.13 | PRDX1 |
| 4 | ENST00000585910.5 | 1.398881e-13 | 3.448465e-11 | -3.687907 | 23.070891 | False | 4.299709e-01 | -3.687907 | 0.498608 | -7.304680 | ENSG00000105254.11 | TBCB |

## Exon skipping and intron retention

Swan can detect novel \(unannotated\) exon skipping and intron retention events.

To obtain a list of genes containing novel exon skipping events, run the following code:

```python
es_genes = sg.find_es_genes()
print(es_genes[:5])
```

```text
Analyzing 893 intronic edges for ES
Found 1021 novel es events from 285 genes.
['ENSG00000130706.12', 'ENSG00000111237.18', 'ENSG00000101363.12', 'ENSG00000163069.12', 'ENSG00000132677.12']
```

As usual, we can feed `es_genes` into `gen_report()` or individual gene ids from `es_genes` into `plot_graph()` to generate gene reports or gene summary graphs respectively.

To obtain a list of genes containing novel intron retention events, run the following code:

```python
ir_genes = sg.find_ir_genes()
print(ir_genes[:5])
```

```text
Analyzing 2185 exonic edges for IR
Found 73 novel ir events from 47 genes.
['ENSG00000213719.8', 'ENSG00000135480.15', 'ENSG00000196421.8', 'ENSG00000141505.11', 'ENSG00000030582.17']
```

As usual, we can feed `ir_genes` into `gen_report()` or individual gene ids from `ir_genes` into `plot_graph()` to generate gene reports or gene summary graphs respectively.

