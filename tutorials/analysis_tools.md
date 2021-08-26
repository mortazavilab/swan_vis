# Analysis tools

Swan has several analysis options to use.
* [Differential gene expression](analysis_tools.md#differential-gene-expression-tests)
* [Differential transcript expression](analysis_tools.md#differential-transcript-expression-tests)
* [Isoform switching / differential isoform expression](analysis_tools.md#isoform-switching-differential-isoform-expression-testing)
* [Combining metadata columns](analysis_tools.md#combining-metadata-columns)
* [Exon skipping and intron retention](analysis_tools.md#exon-skipping-and-intron-retention)
<!-- * [More differential expression](#more_de) -->

<!-- Running this tutorial on my laptop took around 30 minutes and 3 GB of RAM. The longest steps by far are running the differential gene and transcript expression tools. The diffxpy tools are multithreaded, and my laptop has 8 cores. -->


```python
import swan_vis as swan

sg = swan.read('data/swan.p')
```

    Read in graph from data/swan.p


## <a name="deg"></a>Differential gene expression tests

Differential gene expression testing in Swan is implemented via [diffxpy](https://github.com/theislab/diffxpy). To run the test, choose a metadata column from `sg.adata.obs` to test on using the `obs_col` argument. If there are more than 2 unique values in this column, further specify the conditions to test using the `obs_conditions` arguments.

The differential expression test that is run is [diffxpy's Wald test](https://diffxpy.readthedocs.io/en/latest/api/diffxpy.api.test.wald.html#diffxpy.api.test.wald), which checks if a "a certain coefficient introduces a significant difference in the expression of a gene". This test is performed on the normalized TPM for each gene.

<!-- For individuals wanting to run a different diffxpy differential test, see [this section](#more_de). -->


```python
obs_col = 'cell_line'
obs_conditions = ['hepg2', 'hffc6']

# perform a differential gene expression
# Wald test on the provided metadata column and conditions
test = sg.de_gene_test(obs_col, obs_conditions=obs_conditions)
```

The output in `test` is a summary table for the differential expression test.


```python
test.head(2)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gid</th>
      <th>pval</th>
      <th>qval</th>
      <th>log2fc</th>
      <th>mean</th>
      <th>zero_mean</th>
      <th>grad</th>
      <th>coef_mle</th>
      <th>coef_sd</th>
      <th>ll</th>
      <th>gname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>68444</th>
      <td>ENSG00000137204.14</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-297.776029</td>
      <td>2.873145</td>
      <td>False</td>
      <td>1.328783e-08</td>
      <td>-297.776029</td>
      <td>2.222759e-162</td>
      <td>-5.835233</td>
      <td>SLC22A7</td>
    </tr>
    <tr>
      <th>132616</th>
      <td>ENSG00000186204.14</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-297.776029</td>
      <td>9.517887</td>
      <td>False</td>
      <td>2.018855e-07</td>
      <td>-297.776029</td>
      <td>2.222759e-162</td>
      <td>-5.016340</td>
      <td>CYP4F12</td>
    </tr>
  </tbody>
</table>
</div>



The results of this test are also stored in an automatically-generated key in `sg.adata.uns`, and will be saved to the SwanGraph if you save it. You can regenerate this key and access the summary table by running the following code:


```python
# deg - differential gene expression
uns_key = swan.make_uns_key('deg',
                            obs_col=obs_col,
                            obs_conditions=obs_conditions)
test = sg.adata.uns[uns_key]
test.head(2)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gid</th>
      <th>pval</th>
      <th>qval</th>
      <th>log2fc</th>
      <th>mean</th>
      <th>zero_mean</th>
      <th>grad</th>
      <th>coef_mle</th>
      <th>coef_sd</th>
      <th>ll</th>
      <th>gname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>68444</th>
      <td>ENSG00000137204.14</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-297.776029</td>
      <td>2.873145</td>
      <td>False</td>
      <td>1.328783e-08</td>
      <td>-297.776029</td>
      <td>2.222759e-162</td>
      <td>-5.835233</td>
      <td>SLC22A7</td>
    </tr>
    <tr>
      <th>132616</th>
      <td>ENSG00000186204.14</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-297.776029</td>
      <td>9.517887</td>
      <td>False</td>
      <td>2.018855e-07</td>
      <td>-297.776029</td>
      <td>2.222759e-162</td>
      <td>-5.016340</td>
      <td>CYP4F12</td>
    </tr>
  </tbody>
</table>
</div>



Swan can also automatically subset the test summary table to pull out genes that pass a certain significance threshold.


```python
# return a table of significantly differentially-expressed genes
# for a given q val + log2fc threshold
de_genes = sg.get_de_genes(obs_col, obs_conditions=obs_conditions,
                           q=0.05, log2fc=1)
```


```python
de_genes.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gid</th>
      <th>pval</th>
      <th>qval</th>
      <th>log2fc</th>
      <th>mean</th>
      <th>zero_mean</th>
      <th>grad</th>
      <th>coef_mle</th>
      <th>coef_sd</th>
      <th>ll</th>
      <th>gname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>57195</th>
      <td>ENSG00000129245.11</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>283.913085</td>
      <td>7.659833</td>
      <td>False</td>
      <td>1.200001</td>
      <td>283.913085</td>
      <td>2.222759e-162</td>
      <td>-59.551343</td>
      <td>FXR2</td>
    </tr>
    <tr>
      <th>56566</th>
      <td>ENSG00000128656.13</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>283.913085</td>
      <td>81.457545</td>
      <td>False</td>
      <td>1.199999</td>
      <td>283.913085</td>
      <td>2.222759e-162</td>
      <td>-74.030688</td>
      <td>CHN1</td>
    </tr>
    <tr>
      <th>101233</th>
      <td>ENSG00000164318.17</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>283.913085</td>
      <td>4.777928</td>
      <td>False</td>
      <td>0.800000</td>
      <td>283.913085</td>
      <td>2.222759e-162</td>
      <td>-39.453406</td>
      <td>EGFLAM</td>
    </tr>
    <tr>
      <th>57196</th>
      <td>ENSG00000129245.11</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>283.913085</td>
      <td>7.659833</td>
      <td>False</td>
      <td>1.200001</td>
      <td>283.913085</td>
      <td>2.222759e-162</td>
      <td>-59.551343</td>
      <td>FXR2</td>
    </tr>
    <tr>
      <th>101230</th>
      <td>ENSG00000164318.17</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>283.913085</td>
      <td>4.777928</td>
      <td>False</td>
      <td>0.800000</td>
      <td>283.913085</td>
      <td>2.222759e-162</td>
      <td>-39.453406</td>
      <td>EGFLAM</td>
    </tr>
  </tbody>
</table>
</div>



## <a name="det"></a>Differential transcript expression tests

Similarly, Swan can run tests to find differentially expressed transcript isoforms. The input and output to these functions are identical to that of the differential gene tests.

The differential expression test that is run is [diffxpy's Wald test](https://diffxpy.readthedocs.io/en/latest/api/diffxpy.api.test.wald.html#diffxpy.api.test.wald), which checks if a "a certain coefficient introduces a significant difference in the expression of a transcript". This test is performed on the normalized TPM for each transcript.

<!-- For individuals wanting to run a different diffxpy differential test, see [this section](#more_de) -->


```python
obs_col = 'cell_line'
obs_conditions = ['hepg2', 'hffc6']

# perform a differential transcript expression
# Wald test on the provided metadata column and conditions
test = sg.de_transcript_test(obs_col, obs_conditions=obs_conditions)
```



The output in `test` is a summary table for the differential expression test.


```python
test.head(2)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>tid</th>
      <th>pval</th>
      <th>qval</th>
      <th>log2fc</th>
      <th>mean</th>
      <th>zero_mean</th>
      <th>grad</th>
      <th>coef_mle</th>
      <th>coef_sd</th>
      <th>ll</th>
      <th>gid</th>
      <th>gname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>21101</th>
      <td>ENST00000367818.3</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-297.776029</td>
      <td>0.400283</td>
      <td>False</td>
      <td>0.621579</td>
      <td>-297.776029</td>
      <td>2.222759e-162</td>
      <td>0.0</td>
      <td>ENSG00000143184.4</td>
      <td>XCL1</td>
    </tr>
    <tr>
      <th>136458</th>
      <td>ENST00000544590.1</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-297.776029</td>
      <td>0.235725</td>
      <td>False</td>
      <td>0.389695</td>
      <td>-297.776029</td>
      <td>2.222759e-162</td>
      <td>0.0</td>
      <td>ENSG00000109920.12</td>
      <td>FNBP4</td>
    </tr>
  </tbody>
</table>
</div>



The results of this test are similarly stored in an automatically-generated key in `sg.adata.uns`, and will be saved to the SwanGraph if you save it. You can regenerate this key and access the summary table by running the following code:


```python
# det - differential transcript expression
uns_key = swan.make_uns_key('det',
                            obs_col=obs_col,
                            obs_conditions=obs_conditions)
test = sg.adata.uns[uns_key]
test.head(2)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>tid</th>
      <th>pval</th>
      <th>qval</th>
      <th>log2fc</th>
      <th>mean</th>
      <th>zero_mean</th>
      <th>grad</th>
      <th>coef_mle</th>
      <th>coef_sd</th>
      <th>ll</th>
      <th>gid</th>
      <th>gname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>21101</th>
      <td>ENST00000367818.3</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-297.776029</td>
      <td>0.400283</td>
      <td>False</td>
      <td>0.621579</td>
      <td>-297.776029</td>
      <td>2.222759e-162</td>
      <td>0.0</td>
      <td>ENSG00000143184.4</td>
      <td>XCL1</td>
    </tr>
    <tr>
      <th>136458</th>
      <td>ENST00000544590.1</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-297.776029</td>
      <td>0.235725</td>
      <td>False</td>
      <td>0.389695</td>
      <td>-297.776029</td>
      <td>2.222759e-162</td>
      <td>0.0</td>
      <td>ENSG00000109920.12</td>
      <td>FNBP4</td>
    </tr>
  </tbody>
</table>
</div>



Again, Swan can also automatically subset the test summary table to pull out genes that pass a certain significance threshold.


```python
# return a table of significantly differentially-expressed genes
# for a given q val + log2fc threshold
de_transcripts = sg.get_de_transcripts(obs_col, obs_conditions=obs_conditions,
                           q=0.05, log2fc=1)
```


```python
de_transcripts.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>tid</th>
      <th>pval</th>
      <th>qval</th>
      <th>log2fc</th>
      <th>mean</th>
      <th>zero_mean</th>
      <th>grad</th>
      <th>coef_mle</th>
      <th>coef_sd</th>
      <th>ll</th>
      <th>gid</th>
      <th>gname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>91026</th>
      <td>ENST00000486541.1</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>283.913085</td>
      <td>7.222791</td>
      <td>False</td>
      <td>1.20000</td>
      <td>283.913085</td>
      <td>2.222759e-162</td>
      <td>-59.116265</td>
      <td>ENSG00000117318.8</td>
      <td>ID3</td>
    </tr>
    <tr>
      <th>81775</th>
      <td>ENST00000475122.1</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>283.913085</td>
      <td>0.904308</td>
      <td>False</td>
      <td>0.80000</td>
      <td>283.913085</td>
      <td>2.222759e-162</td>
      <td>-32.247542</td>
      <td>ENSG00000119812.18</td>
      <td>FAM98A</td>
    </tr>
    <tr>
      <th>91246</th>
      <td>ENST00000486828.6</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>283.913085</td>
      <td>0.253818</td>
      <td>False</td>
      <td>0.83105</td>
      <td>283.913085</td>
      <td>2.222759e-162</td>
      <td>0.000000</td>
      <td>ENSG00000196923.13</td>
      <td>PDLIM7</td>
    </tr>
    <tr>
      <th>190208</th>
      <td>ENST00000623250.1</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>283.913085</td>
      <td>0.841705</td>
      <td>False</td>
      <td>1.20000</td>
      <td>283.913085</td>
      <td>2.222759e-162</td>
      <td>-45.109253</td>
      <td>ENSG00000279348.1</td>
      <td>AC012513.3</td>
    </tr>
    <tr>
      <th>91032</th>
      <td>ENST00000486554.1</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>283.913085</td>
      <td>0.507635</td>
      <td>False</td>
      <td>0.40000</td>
      <td>283.913085</td>
      <td>2.222759e-162</td>
      <td>-16.496230</td>
      <td>ENSG00000157514.16</td>
      <td>TSC22D3</td>
    </tr>
  </tbody>
</table>
</div>



## <a name="is"></a>Isoform switching / Differential isoform expression testing

Isoform switching / differential isoform expression (DIE) testing is implemented according to the strategy in [Joglekar et. al., 2021](https://www.nature.com/articles/s41467-020-20343-5). DIE can roughly be described as finding statistically significant changes in isoform expression between two conditions along with a change in percent isoform usage per gene.

Pairwise comparisons can be set up using different columns in the metadata that was added to the SwanGraph with the `obs_col` and `obs_conditions` arguments.


```python
# look at valid metadata options
sg.adata.obs
```




<div>
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




```python
# find genes that exhibit DIE between HFFc6 and HepG2
obs_col = 'cell_line'
obs_conditions = ['hepg2', 'hffc6']
die_table = sg.die_gene_test(obs_col=obs_col,
                             obs_conditions=obs_conditions,
                             verbose=True)
```

    Testing for DIE for each gene: 100%|█████████▉| 58905/58906 [12:38<00:00, 87.74it/s]

The resultant table contains an entry for each gene with the p value (`p_val`), adjusted p value (`adj_p_val`), and change in percent isoform usage for the top two isoforms (`dpi`). Exact details on these calculations can be found in [Joglekar et. al., 2021](https://www.nature.com/articles/s41467-020-20343-5).


```python
die_table.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gid</th>
      <th>p_val</th>
      <th>dpi</th>
      <th>adj_p_val</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000004059.10</td>
      <td>0.370956</td>
      <td>0.659836</td>
      <td>0.545069</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000003509.15</td>
      <td>0.110749</td>
      <td>14.178416</td>
      <td>0.236859</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ENSG00000001630.16</td>
      <td>0.199904</td>
      <td>1.194254</td>
      <td>0.360028</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ENSG00000001461.16</td>
      <td>0.717673</td>
      <td>9.237480</td>
      <td>0.823473</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ENSG00000005801.17</td>
      <td>0.004018</td>
      <td>33.102501</td>
      <td>0.017566</td>
    </tr>
  </tbody>
</table>
</div>



As with differential expression testing, differential isoform expression testing results are stored automatically in `sg.adata.uns`, and will be saved to the SwanGraph if you save it. You can regenerate this key and access the summary table by running the following code:


```python
# die_iso - isoform level differential isoform expression test results
uns_key = swan.make_uns_key('die',
                            obs_col=obs_col,
                            obs_conditions=obs_conditions)
test = sg.adata.uns[uns_key]
test.head(2)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gid</th>
      <th>p_val</th>
      <th>dpi</th>
      <th>adj_p_val</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000004059.10</td>
      <td>0.370956</td>
      <td>0.659836</td>
      <td>0.545069</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000003509.15</td>
      <td>0.110749</td>
      <td>14.178416</td>
      <td>0.236859</td>
    </tr>
  </tbody>
</table>
</div>



Swan comes with an easy way to filter your DIE test results based on adjusted p value and dpi thresholds.


```python
test = sg.get_die_genes(obs_col=obs_col, obs_conditions=obs_conditions,
                       p=0.05, dpi=10)
test.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gid</th>
      <th>p_val</th>
      <th>dpi</th>
      <th>adj_p_val</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>4</th>
      <td>ENSG00000005801.17</td>
      <td>4.017539e-03</td>
      <td>33.102501</td>
      <td>0.017566</td>
    </tr>
    <tr>
      <th>6</th>
      <td>ENSG00000075790.10</td>
      <td>3.733015e-06</td>
      <td>34.215599</td>
      <td>0.000038</td>
    </tr>
    <tr>
      <th>8</th>
      <td>ENSG00000005175.9</td>
      <td>6.216589e-03</td>
      <td>14.736797</td>
      <td>0.025327</td>
    </tr>
    <tr>
      <th>11</th>
      <td>ENSG00000006282.20</td>
      <td>1.472361e-04</td>
      <td>19.920383</td>
      <td>0.001061</td>
    </tr>
    <tr>
      <th>13</th>
      <td>ENSG00000007376.7</td>
      <td>1.989310e-07</td>
      <td>29.276819</td>
      <td>0.000003</td>
    </tr>
  </tbody>
</table>
</div>



Swan also now automatically tracks transcription start site (TSS) and transcription end site (TES) usage, and find genes that exhibit DIE on the basis of their starts or ends. To do this, use the `kind` argument to `die_gene_test`.


```python
# find genes that exhibit DIE for TSSs between HFFc6 and HepG2
die_table = sg.die_gene_test(kind='tss',
                             obs_col=obs_col,
                             obs_conditions=obs_conditions,
                             verbose=True)
die_table.head()
```

    Testing for DIE for each gene: 100%|██████████| 58906/58906 [12:38<00:00, 77.66it/s]
    Testing for DIE for each gene: 100%|█████████▉| 58896/58906 [09:30<00:00, 108.98it/s]




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gid</th>
      <th>p_val</th>
      <th>dpi</th>
      <th>adj_p_val</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000001630.16</td>
      <td>0.286866</td>
      <td>2.580168</td>
      <td>0.710920</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000002330.13</td>
      <td>0.089441</td>
      <td>12.817429</td>
      <td>0.348247</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ENSG00000002586.19</td>
      <td>0.885980</td>
      <td>0.102093</td>
      <td>0.982149</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ENSG00000002822.15</td>
      <td>0.678511</td>
      <td>1.724138</td>
      <td>0.962778</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ENSG00000002919.14</td>
      <td>0.669412</td>
      <td>2.439026</td>
      <td>0.962778</td>
    </tr>
  </tbody>
</table>
</div>



To access the die results on the tss level, use `die_kind='tss'` as input to `make_uns_key()`.


```python
# die_iso - TSS level differential isoform expression test results
uns_key = swan.make_uns_key('die',
                            obs_col=obs_col,
                            obs_conditions=obs_conditions,
                            die_kind='tss')
test = sg.adata.uns[uns_key]
test.head(2)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gid</th>
      <th>p_val</th>
      <th>dpi</th>
      <th>adj_p_val</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000001630.16</td>
      <td>0.286866</td>
      <td>2.580168</td>
      <td>0.710920</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000002330.13</td>
      <td>0.089441</td>
      <td>12.817429</td>
      <td>0.348247</td>
    </tr>
  </tbody>
</table>
</div>



And provide the `kind='tss'` option to `get_die_genes()` when trying to filter your test results.


```python
test = sg.get_die_genes(kind='tss', obs_col=obs_col,
                        obs_conditions=obs_conditions,
                        p=0.05, dpi=10)
test.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gid</th>
      <th>p_val</th>
      <th>dpi</th>
      <th>adj_p_val</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>5</th>
      <td>ENSG00000003402.19</td>
      <td>4.039839e-06</td>
      <td>39.797981</td>
      <td>6.987697e-05</td>
    </tr>
    <tr>
      <th>6</th>
      <td>ENSG00000003436.15</td>
      <td>1.183656e-05</td>
      <td>28.073771</td>
      <td>1.856128e-04</td>
    </tr>
    <tr>
      <th>7</th>
      <td>ENSG00000004487.16</td>
      <td>1.135407e-03</td>
      <td>75.714287</td>
      <td>1.145036e-02</td>
    </tr>
    <tr>
      <th>36</th>
      <td>ENSG00000008952.16</td>
      <td>4.075852e-03</td>
      <td>26.086960</td>
      <td>3.342667e-02</td>
    </tr>
    <tr>
      <th>40</th>
      <td>ENSG00000010278.13</td>
      <td>3.077503e-22</td>
      <td>39.024387</td>
      <td>2.195799e-20</td>
    </tr>
  </tbody>
</table>
</div>



For TESs, use `kind='tes'` as input to `die_genes_test()`, `die_kind='tes'` to `make_uns_key()`, and `kind='tes'` to `get_die_genes()`.


```python
# find genes that exhibit DIE for TESs between HFFc6 and HepG2
die_table = sg.die_gene_test(kind='tes', obs_col='cell_line', obs_conditions=['hepg2', 'hffc6'])
die_table.head()
```

    Testing for DIE for each gene: 100%|██████████| 58906/58906 [09:31<00:00, 103.12it/s]





<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gid</th>
      <th>p_val</th>
      <th>dpi</th>
      <th>adj_p_val</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000000419.12</td>
      <td>0.749266</td>
      <td>0.096133</td>
      <td>0.904272</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000001461.16</td>
      <td>0.852914</td>
      <td>9.093739</td>
      <td>0.943997</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ENSG00000001630.16</td>
      <td>0.286866</td>
      <td>2.580168</td>
      <td>0.613824</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ENSG00000002330.13</td>
      <td>0.184635</td>
      <td>12.817430</td>
      <td>0.479316</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ENSG00000002549.12</td>
      <td>0.679148</td>
      <td>0.694543</td>
      <td>0.879387</td>
    </tr>
  </tbody>
</table>
</div>




```python
# die_iso - TSS level differential isoform expression test results
uns_key = swan.make_uns_key('die',
                            obs_col=obs_col,
                            obs_conditions=obs_conditions,
                            die_kind='tes')
test = sg.adata.uns[uns_key]
test.head(2)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gid</th>
      <th>p_val</th>
      <th>dpi</th>
      <th>adj_p_val</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000000419.12</td>
      <td>0.749266</td>
      <td>0.096133</td>
      <td>0.904272</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000001461.16</td>
      <td>0.852914</td>
      <td>9.093739</td>
      <td>0.943997</td>
    </tr>
  </tbody>
</table>
</div>




```python
test = sg.get_die_genes(kind='tes', obs_col=obs_col,
                        obs_conditions=obs_conditions,
                        p=0.05, dpi=10)
test.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gid</th>
      <th>p_val</th>
      <th>dpi</th>
      <th>adj_p_val</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>9</th>
      <td>ENSG00000003402.19</td>
      <td>7.338098e-15</td>
      <td>84.848488</td>
      <td>5.296535e-13</td>
    </tr>
    <tr>
      <th>10</th>
      <td>ENSG00000003436.15</td>
      <td>2.772096e-07</td>
      <td>32.637852</td>
      <td>7.003007e-06</td>
    </tr>
    <tr>
      <th>14</th>
      <td>ENSG00000004487.16</td>
      <td>1.135407e-03</td>
      <td>75.714287</td>
      <td>1.141621e-02</td>
    </tr>
    <tr>
      <th>32</th>
      <td>ENSG00000006282.20</td>
      <td>6.858641e-03</td>
      <td>19.819595</td>
      <td>4.915014e-02</td>
    </tr>
    <tr>
      <th>60</th>
      <td>ENSG00000010278.13</td>
      <td>3.191221e-22</td>
      <td>39.024387</td>
      <td>3.486194e-20</td>
    </tr>
  </tbody>
</table>
</div>



## <a name="multi_gb"></a>Combining metadata columns

What if none of the metadata columns you have summarize the comparisons you want to make? What if I want to find differentially expressed genes or transcripts, or find isoform-switching genes between hffc6 replicate 3 and hepg2 replicate 1?


```python
sg.adata.obs.head()
```




<div>
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



Let's ignore for a moment the fact that the `dataset` column does effectively capture both replicate as well as cell line metadata. This may not always be the case with more complex datasets! Swan has a function to concatenate columns together and add them as an additional column to the metadata tables. Use the following code to generate a new column that concatenates as many preexisting metadata columns as you wish:


```python
col_name = sg.add_multi_groupby(['cell_line', 'replicate'])

print(col_name)
print(sg.adata.obs.head())
```

    cell_line_replicate
             dataset cell_line replicate cell_line_replicate
    index                                                   
    hepg2_1  hepg2_1     hepg2         1             hepg2_1
    hepg2_2  hepg2_2     hepg2         2             hepg2_2
    hffc6_1  hffc6_1     hffc6         1             hffc6_1
    hffc6_2  hffc6_2     hffc6         2             hffc6_2
    hffc6_3  hffc6_3     hffc6         3             hffc6_3


The added column in `col_name` can then be used as the `obs_col` input to `de_gene_test(), de_transcript_test(), and die_gene_test()`, as in the following calls:


```python
obs_col = col_name
obs_conditions = ['hffc6_3', 'hepg2_1']

deg_summary = sg.de_gene_test(obs_col=obs_col,
                              obs_conditions=obs_conditions)
det_summary = sg.de_transcript_test(obs_col=obs_col,
                                    obs_conditions=obs_conditions)
die_summary = sg.die_gene_test(obs_col=obs_col,
                               obs_conditions=obs_conditions)
```

## <a name="es_ir"></a>Exon skipping and intron retention

Swan can detect novel (unannotated) exon skipping and intron retention events.

To obtain a dataframe of novel exon skipping events, run the following code:


```python
# returns a DataFrame of genes, transcripts, and specific edges in
# the SwanGraph with novel exon skipping events
es_df = sg.find_es_genes(verbose=True)
```

    Testing each novel edge for exon skipping: 100%|██████████| 855/855 [1:26:06<00:00,  6.13s/it]

    Found 529 novel es events in 149 transcripts.



```python
es_df.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gid</th>
      <th>tid</th>
      <th>egde_id</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000157916.19</td>
      <td>TALONT000218256</td>
      <td>952616</td>
    </tr>
    <tr>
      <th>0</th>
      <td>ENSG00000122406.13</td>
      <td>TALONT000425229</td>
      <td>952716</td>
    </tr>
    <tr>
      <th>0</th>
      <td>ENSG00000224093.5</td>
      <td>TALONT000434035</td>
      <td>952720</td>
    </tr>
    <tr>
      <th>0</th>
      <td>ENSG00000224093.5</td>
      <td>TALONT000434035</td>
      <td>952720</td>
    </tr>
    <tr>
      <th>0</th>
      <td>ENSG00000224093.5</td>
      <td>TALONT000434035</td>
      <td>952720</td>
    </tr>
  </tbody>
</table>
</div>



You can pass gene IDs from `es_df` into `gen_report()` or `plot_graph()` to visualize where these exon-skipping events are in gene reports or gene summary graphs respectively.

To obtain a list of genes containing novel intron retention events, run the following code:


```python
# returns a DataFrame of genes, transcripts, and specific edges in
# the SwanGraph with novel intron retaining events
ir_df = sg.find_ir_genes(verbose=True)
```

    Testing each novel edge for intron retention: 100%|██████████| 1186/1186 [1:59:16<00:00,  5.97s/it][A

    Found 35 novel ir events in 27 transcripts.


```python
ir_df.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gid</th>
      <th>tid</th>
      <th>egde_id</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000143753.12</td>
      <td>TALONT000482711</td>
      <td>952811</td>
    </tr>
    <tr>
      <th>0</th>
      <td>ENSG00000285053.1</td>
      <td>TALONT000483978</td>
      <td>952821</td>
    </tr>
    <tr>
      <th>0</th>
      <td>ENSG00000177042.14</td>
      <td>TALONT000213980</td>
      <td>954058</td>
    </tr>
    <tr>
      <th>0</th>
      <td>ENSG00000177042.14</td>
      <td>TALONT000213980</td>
      <td>954058</td>
    </tr>
    <tr>
      <th>0</th>
      <td>ENSG00000148926.9</td>
      <td>TALONT000251937</td>
      <td>954085</td>
    </tr>
  </tbody>
</table>
</div>


You can pass gene IDs from `ir_df` into `gen_report()` or `plot_graph()` to visualize where these intron retention events are in gene reports or gene summary graphs respectively.

<!-- ## <a name="more_de"></a>More differential expression

For users that are interested in using different differential expression tests, or tweaking the input parameters, we encourage them to obtain an AnnData version of of their SwanGraph using `create_gene_anndata` or `create_transcript_anndata`, and exploring the numerous differential testing options that diffxpy supports.

[Diffxpy differential testing tutorials](https://diffxpy.readthedocs.io/en/latest/tutorials.html#differential-testing)

[More information on diffxpy differential expression tests](https://nbviewer.jupyter.org/github/theislab/diffxpy_tutorials/blob/master/diffxpy_tutorials/test/introduction_differential_testing.ipynb)


```python
dataset_groups = [['HepG2_1','HepG2_2'],
                  ['HFFc6_1','HFFc6_2','HFFc6_3']]

# create a gene-level AnnData object compatible with diffxpy
# that assigns different condition labels to the given dataset groups
gene_adata = sg.create_gene_anndata(dataset_groups)
```

    Transforming to str index.



```python
# create a transcript-level AnnData object compatible with diffxpy
# that assigns different condition labels to the given dataset groups
transcript_adata = sg.create_transcript_anndata(dataset_groups)
```

    Transforming to str index. -->
