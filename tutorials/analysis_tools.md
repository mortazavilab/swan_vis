# Analysis tools

Swan has several analysis options to use.
* [Differential expression tests](#de)
* [Isoform switching / differential isoform expression](#is)
* [Combining metadata columns](#multi_gb)
* [Exon skipping and intron retention](#es_ir)
* [More differential expression](#more_de)

<!-- Running this tutorial on my laptop took around 30 minutes and 3 GB of RAM. The longest steps by far are running the differential gene and transcript expression tools. The diffxpy tools are multithreaded, and my laptop has 8 cores. -->


```python
import swan_vis as swan

sg = swan.read('data/swan.p')
```

    Read in graph from data/swan.p


## <a name="de"></a>Differential expression tests

Swan's old differential gene and transcript expression tests using `diffxpy` have now been deprecated as it seems that the library is unsupported. I recommend that users interested in running differential gene or transcript expression test either use [PyDESeq2](https://github.com/owkin/PyDESeq2) or Scanpy's [`rank_genes_groups`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html) test, which both support the AnnData format for simple compatibility with Swan's AnnData expression representation.

<!-- ### <a name="de"></a>Using Scanpy's `rank_genes_groups` -->

<!-- `rank_genes_groups` expects logarithmized data, so make sure you transform your data before running the test on whichever [AnnData](https://freese.gitbook.io/swan/faqs/data_structure#anndata) you want that's in your SwanGraph. -->

<!-- ```python
import scanpy as sc

sg.adata.X = sg.adata.layers['tpm']
sc.pp.log1p(sg.adata)
sg.adata.layers['log_norm'] = sg.adata.X.copy()
sc.tl.rank_genes_groups(sg.adata,
                        groupby=<obs_col>,
                        groups=<obs_conditions>,
                        layer='log_norm',
                        method='wilcoxon')

results_df = sc.get.rank_genes_groups_df(sg.adata, <obs_condition>)
results_df.head()
``` -->

### <a name="de"></a>Using PyDESeq2

Please read the [PyDESeq2 documentation](https://pydeseq2.readthedocs.io/en/latest/) for details on how to use one of the SwanGraph AnnData objects to obtain differential expression results. Below is an example on how to find differentially expressed transcripts between cell lines.


```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import numpy as np

adata = sg.adata.copy()

# PyDESeq2 currently doesn't support column names with underscores, so change that
adata.obs.rename({'cell_line': 'cellline'}, axis=1, inplace=True)
obs_col = 'cellline'

threads = 8

# densify matrix
adata.X = np.array(adata.X.todense())

# run test
dds = DeseqDataSet(adata=adata,
               design_factors=obs_col,
               n_cpus=threads,
               refit_cooks=True)
dds.deseq2()
stat_res = DeseqStats(dds,
                  n_cpus=threads)
stat_res.summary()

df = stat_res.results_df
```

    Fitting size factors...
    ... done in 0.00 seconds.

    Fitting dispersions...
    ... done in 34.09 seconds.

    Fitting dispersion trend curve...
    ... done in 12.01 seconds.

    Fitting MAP dispersions...
    ... done in 13.97 seconds.

    Fitting LFCs...
    ... done in 4.86 seconds.

    Refitting 0 outliers.

    Running Wald tests...
    ... done in 2.36 seconds.

    Log2 fold change & Wald test p-value: cellline hffc6 vs hepg2



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
      <th>baseMean</th>
      <th>log2FoldChange</th>
      <th>lfcSE</th>
      <th>stat</th>
      <th>pvalue</th>
      <th>padj</th>
    </tr>
    <tr>
      <th>tid</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>TALONT000296400</th>
      <td>3.497574</td>
      <td>2.416436</td>
      <td>1.472616</td>
      <td>1.640914</td>
      <td>0.100815</td>
      <td>0.193840</td>
    </tr>
    <tr>
      <th>ENST00000591581.1</th>
      <td>0.162115</td>
      <td>0.624956</td>
      <td>5.002765</td>
      <td>0.124922</td>
      <td>0.900585</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>ENST00000546893.5</th>
      <td>8.173445</td>
      <td>-0.053334</td>
      <td>0.793179</td>
      <td>-0.067241</td>
      <td>0.946390</td>
      <td>0.968527</td>
    </tr>
    <tr>
      <th>ENST00000537289.1</th>
      <td>5.140216</td>
      <td>-1.248175</td>
      <td>0.978927</td>
      <td>-1.275044</td>
      <td>0.202294</td>
      <td>0.328286</td>
    </tr>
    <tr>
      <th>ENST00000258382.9</th>
      <td>2.012104</td>
      <td>-0.011819</td>
      <td>1.493866</td>
      <td>-0.007912</td>
      <td>0.993687</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>ENST00000506914.1</th>
      <td>1.022692</td>
      <td>0.947588</td>
      <td>2.272813</td>
      <td>0.416923</td>
      <td>0.676735</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>ENST00000571080.1</th>
      <td>0.473263</td>
      <td>-2.824409</td>
      <td>3.426473</td>
      <td>-0.824291</td>
      <td>0.409774</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>ENST00000378615.7</th>
      <td>1.334873</td>
      <td>-1.199810</td>
      <td>1.790781</td>
      <td>-0.669993</td>
      <td>0.502862</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>ENST00000409586.7</th>
      <td>0.338615</td>
      <td>1.464364</td>
      <td>4.023896</td>
      <td>0.363917</td>
      <td>0.715920</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>ENST00000370278.3</th>
      <td>1.328549</td>
      <td>-1.201246</td>
      <td>1.962358</td>
      <td>-0.612144</td>
      <td>0.540443</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
<p>34814 rows × 6 columns</p>
</div>



```python
df.head()
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
      <th>baseMean</th>
      <th>log2FoldChange</th>
      <th>lfcSE</th>
      <th>stat</th>
      <th>pvalue</th>
      <th>padj</th>
    </tr>
    <tr>
      <th>tid</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>TALONT000296400</th>
      <td>3.497574</td>
      <td>2.416436</td>
      <td>1.472616</td>
      <td>1.640914</td>
      <td>0.100815</td>
      <td>0.193840</td>
    </tr>
    <tr>
      <th>ENST00000591581.1</th>
      <td>0.162115</td>
      <td>0.624956</td>
      <td>5.002765</td>
      <td>0.124922</td>
      <td>0.900585</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>ENST00000546893.5</th>
      <td>8.173445</td>
      <td>-0.053334</td>
      <td>0.793179</td>
      <td>-0.067241</td>
      <td>0.946390</td>
      <td>0.968527</td>
    </tr>
    <tr>
      <th>ENST00000537289.1</th>
      <td>5.140216</td>
      <td>-1.248175</td>
      <td>0.978927</td>
      <td>-1.275044</td>
      <td>0.202294</td>
      <td>0.328286</td>
    </tr>
    <tr>
      <th>ENST00000258382.9</th>
      <td>2.012104</td>
      <td>-0.011819</td>
      <td>1.493866</td>
      <td>-0.007912</td>
      <td>0.993687</td>
      <td>NaN</td>
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
    </tr>
    <tr>
      <th>dataset</th>
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
    </tr>
    <tr>
      <th>hepg2_2</th>
      <td>hepg2</td>
      <td>2</td>
      <td>hepg2_2</td>
      <td>848447.0</td>
    </tr>
    <tr>
      <th>hffc6_1</th>
      <td>hffc6</td>
      <td>1</td>
      <td>hffc6_1</td>
      <td>761493.0</td>
    </tr>
    <tr>
      <th>hffc6_2</th>
      <td>hffc6</td>
      <td>2</td>
      <td>hffc6_2</td>
      <td>787967.0</td>
    </tr>
    <tr>
      <th>hffc6_3</th>
      <td>hffc6</td>
      <td>3</td>
      <td>hffc6_3</td>
      <td>614921.0</td>
    </tr>
  </tbody>
</table>
</div>




```python
# find genes that exhibit DIE between HFFc6 and HepG2
obs_col = 'cell_line'
obs_conditions = ['hepg2', 'hffc6']
die_table, die_results = sg.die_gene_test(obs_col=obs_col,
                                          obs_conditions=obs_conditions,
                                          verbose=True)
```

    Testing for DIE for each gene: 100%|██████████| 14684/14684 [03:50<00:00, 123.69it/s]

The resultant table contains an entry for each gene with the p value (`p_val`), adjusted p value (`adj_p_val`), and change in percent isoform usage for the top two isoforms (`dpi`), as well as the identities of the top isoforms involved in the switch. Exact details on these calculations can be found in [Joglekar et. al., 2021](https://www.nature.com/articles/s41467-020-20343-5).


```python
die_table.head()
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
      <th>gid</th>
      <th>p_val</th>
      <th>dpi</th>
      <th>pos_iso_1</th>
      <th>pos_iso_2</th>
      <th>pos_iso_1_dpi</th>
      <th>pos_iso_2_dpi</th>
      <th>neg_iso_1</th>
      <th>neg_iso_2</th>
      <th>neg_iso_1_dpi</th>
      <th>neg_iso_2_dpi</th>
      <th>adj_p_val</th>
      <th>gname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000130175</td>
      <td>0.206667</td>
      <td>39.264992</td>
      <td>TALONT000296399</td>
      <td>NaN</td>
      <td>39.264992</td>
      <td>NaN</td>
      <td>TALONT000296400</td>
      <td>ENST00000589838.5</td>
      <td>-20.116056</td>
      <td>-10.638298</td>
      <td>0.469420</td>
      <td>PRKCSH</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000130202</td>
      <td>0.000367</td>
      <td>11.893251</td>
      <td>ENST00000252485.8</td>
      <td>ENST00000252483.9</td>
      <td>9.684967</td>
      <td>2.208281</td>
      <td>TALONT000406668</td>
      <td>ENST00000591581.1</td>
      <td>-11.473083</td>
      <td>-0.420168</td>
      <td>0.003560</td>
      <td>NECTIN2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ENSG00000111371</td>
      <td>0.680435</td>
      <td>9.401713</td>
      <td>ENST00000398637.9</td>
      <td>ENST00000439706.5</td>
      <td>8.547012</td>
      <td>0.854701</td>
      <td>ENST00000546893.5</td>
      <td>NaN</td>
      <td>-9.401711</td>
      <td>NaN</td>
      <td>0.886243</td>
      <td>SLC38A1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ENSG00000181924</td>
      <td>0.028195</td>
      <td>9.452934</td>
      <td>ENST00000537289.1</td>
      <td>ENST00000545127.1</td>
      <td>7.298603</td>
      <td>2.154328</td>
      <td>ENST00000355693.4</td>
      <td>NaN</td>
      <td>-9.452934</td>
      <td>NaN</td>
      <td>0.122619</td>
      <td>COA4</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ENSG00000163468</td>
      <td>0.255788</td>
      <td>0.680048</td>
      <td>ENST00000295688.7</td>
      <td>TALONT000476055</td>
      <td>0.366966</td>
      <td>0.277693</td>
      <td>ENST00000368259.6</td>
      <td>ENST00000489870.1</td>
      <td>-0.568503</td>
      <td>-0.111545</td>
      <td>0.525066</td>
      <td>CCT3</td>
    </tr>
  </tbody>
</table>
</div>



Differential isoform expression testing results are stored automatically in `sg.adata.uns`, and will be saved to the SwanGraph if you save it. You can regenerate this key and access the summary table by running the following code:


```python
# die_iso - isoform level differential isoform expression test results
uns_key = swan.make_uns_key('die',
                            obs_col=obs_col,
                            obs_conditions=obs_conditions)
test = sg.adata.uns[uns_key]
test.head(2)
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
      <th>gid</th>
      <th>p_val</th>
      <th>dpi</th>
      <th>pos_iso_1</th>
      <th>pos_iso_2</th>
      <th>pos_iso_1_dpi</th>
      <th>pos_iso_2_dpi</th>
      <th>neg_iso_1</th>
      <th>neg_iso_2</th>
      <th>neg_iso_1_dpi</th>
      <th>neg_iso_2_dpi</th>
      <th>adj_p_val</th>
      <th>gname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000130175</td>
      <td>0.206667</td>
      <td>39.264992</td>
      <td>TALONT000296399</td>
      <td>NaN</td>
      <td>39.264992</td>
      <td>NaN</td>
      <td>TALONT000296400</td>
      <td>ENST00000589838.5</td>
      <td>-20.116056</td>
      <td>-10.638298</td>
      <td>0.46942</td>
      <td>PRKCSH</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000130202</td>
      <td>0.000367</td>
      <td>11.893251</td>
      <td>ENST00000252485.8</td>
      <td>ENST00000252483.9</td>
      <td>9.684967</td>
      <td>2.208281</td>
      <td>TALONT000406668</td>
      <td>ENST00000591581.1</td>
      <td>-11.473083</td>
      <td>-0.420168</td>
      <td>0.00356</td>
      <td>NECTIN2</td>
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
      <th>p_val</th>
      <th>dpi</th>
      <th>pos_iso_1</th>
      <th>pos_iso_2</th>
      <th>pos_iso_1_dpi</th>
      <th>pos_iso_2_dpi</th>
      <th>neg_iso_1</th>
      <th>neg_iso_2</th>
      <th>neg_iso_1_dpi</th>
      <th>neg_iso_2_dpi</th>
      <th>adj_p_val</th>
      <th>gname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>ENSG00000130202</td>
      <td>3.674234e-04</td>
      <td>11.893251</td>
      <td>ENST00000252485.8</td>
      <td>ENST00000252483.9</td>
      <td>9.684967</td>
      <td>2.208281</td>
      <td>TALONT000406668</td>
      <td>ENST00000591581.1</td>
      <td>-11.473083</td>
      <td>-0.420168</td>
      <td>3.560035e-03</td>
      <td>NECTIN2</td>
    </tr>
    <tr>
      <th>8</th>
      <td>ENSG00000025156</td>
      <td>1.857411e-03</td>
      <td>35.714287</td>
      <td>ENST00000452194.5</td>
      <td>ENST00000465214.2</td>
      <td>25.714287</td>
      <td>10.000000</td>
      <td>ENST00000368455.8</td>
      <td>NaN</td>
      <td>-35.714287</td>
      <td>NaN</td>
      <td>1.405048e-02</td>
      <td>HSF2</td>
    </tr>
    <tr>
      <th>20</th>
      <td>ENSG00000105254</td>
      <td>4.779408e-38</td>
      <td>25.419458</td>
      <td>ENST00000585910.5</td>
      <td>TALONT000366329</td>
      <td>20.394853</td>
      <td>4.702971</td>
      <td>ENST00000221855.7</td>
      <td>ENST00000589996.5</td>
      <td>-25.016304</td>
      <td>-0.403154</td>
      <td>7.343219e-36</td>
      <td>TBCB</td>
    </tr>
    <tr>
      <th>23</th>
      <td>ENSG00000105379</td>
      <td>1.802782e-307</td>
      <td>81.136333</td>
      <td>ENST00000354232.8</td>
      <td>NaN</td>
      <td>81.136333</td>
      <td>NaN</td>
      <td>ENST00000309244.8</td>
      <td>ENST00000596253.1</td>
      <td>-80.857935</td>
      <td>-0.278394</td>
      <td>1.551113e-304</td>
      <td>ETFB</td>
    </tr>
    <tr>
      <th>28</th>
      <td>ENSG00000148180</td>
      <td>0.000000e+00</td>
      <td>89.319039</td>
      <td>ENST00000373818.8</td>
      <td>TALONT000419680</td>
      <td>85.245850</td>
      <td>4.073189</td>
      <td>TALONT000418752</td>
      <td>ENST00000373808.8</td>
      <td>-68.603180</td>
      <td>-7.768958</td>
      <td>0.000000e+00</td>
      <td>GSN</td>
    </tr>
  </tbody>
</table>
</div>



Swan also now automatically tracks transcription start site (TSS) and transcription end site (TES) usage, and can find genes that exhibit DIE on the basis of their starts or ends. To do this, use the `kind` argument to `die_gene_test`.


```python
# find genes that exhibit DIE for TSSs between HFFc6 and HepG2
die_table, die_results = sg.die_gene_test(kind='tss',
                                          obs_col=obs_col,
                                          obs_conditions=obs_conditions,
                                          verbose=True)
die_table.head()
```

    Testing for DIE for each gene: 100%|█████████▉| 14674/14684 [02:18<00:00, 162.60it/s]




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
      <th>p_val</th>
      <th>dpi</th>
      <th>pos_iso_1</th>
      <th>pos_iso_2</th>
      <th>pos_iso_1_dpi</th>
      <th>pos_iso_2_dpi</th>
      <th>neg_iso_1</th>
      <th>neg_iso_2</th>
      <th>neg_iso_1_dpi</th>
      <th>neg_iso_2_dpi</th>
      <th>adj_p_val</th>
      <th>gname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000000419</td>
      <td>0.249561</td>
      <td>10.002187</td>
      <td>ENSG00000000419_2</td>
      <td>ENSG00000000419_1</td>
      <td>9.906052</td>
      <td>0.096133</td>
      <td>ENSG00000000419_3</td>
      <td>ENSG00000000419_4</td>
      <td>-7.218704</td>
      <td>-2.783483</td>
      <td>0.536906</td>
      <td>DPM1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000001461</td>
      <td>0.852914</td>
      <td>9.093739</td>
      <td>ENSG00000001461_5</td>
      <td>NaN</td>
      <td>9.093739</td>
      <td>NaN</td>
      <td>ENSG00000001461_1</td>
      <td>ENSG00000001461_3</td>
      <td>-5.543442</td>
      <td>-1.775148</td>
      <td>1.000000</td>
      <td>NIPAL3</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ENSG00000001497</td>
      <td>0.891496</td>
      <td>3.846153</td>
      <td>ENSG00000001497_1</td>
      <td>NaN</td>
      <td>3.846146</td>
      <td>NaN</td>
      <td>ENSG00000001497_2</td>
      <td>NaN</td>
      <td>-3.846153</td>
      <td>NaN</td>
      <td>1.000000</td>
      <td>LAS1L</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ENSG00000001630</td>
      <td>0.286866</td>
      <td>2.580168</td>
      <td>ENSG00000001630_3</td>
      <td>NaN</td>
      <td>2.580168</td>
      <td>NaN</td>
      <td>ENSG00000001630_2</td>
      <td>ENSG00000001630_1</td>
      <td>-1.549232</td>
      <td>-1.030928</td>
      <td>0.582636</td>
      <td>CYP51A1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ENSG00000002330</td>
      <td>0.184635</td>
      <td>12.817431</td>
      <td>ENSG00000002330_2</td>
      <td>ENSG00000002330_1</td>
      <td>11.868484</td>
      <td>0.948944</td>
      <td>ENSG00000002330_4</td>
      <td>ENSG00000002330_3</td>
      <td>-12.603584</td>
      <td>-0.213847</td>
      <td>0.454053</td>
      <td>BAD</td>
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
      <th>p_val</th>
      <th>dpi</th>
      <th>pos_iso_1</th>
      <th>pos_iso_2</th>
      <th>pos_iso_1_dpi</th>
      <th>pos_iso_2_dpi</th>
      <th>neg_iso_1</th>
      <th>neg_iso_2</th>
      <th>neg_iso_1_dpi</th>
      <th>neg_iso_2_dpi</th>
      <th>adj_p_val</th>
      <th>gname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000000419</td>
      <td>0.249561</td>
      <td>10.002187</td>
      <td>ENSG00000000419_2</td>
      <td>ENSG00000000419_1</td>
      <td>9.906052</td>
      <td>0.096133</td>
      <td>ENSG00000000419_3</td>
      <td>ENSG00000000419_4</td>
      <td>-7.218704</td>
      <td>-2.783483</td>
      <td>0.536906</td>
      <td>DPM1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000001461</td>
      <td>0.852914</td>
      <td>9.093739</td>
      <td>ENSG00000001461_5</td>
      <td>NaN</td>
      <td>9.093739</td>
      <td>NaN</td>
      <td>ENSG00000001461_1</td>
      <td>ENSG00000001461_3</td>
      <td>-5.543442</td>
      <td>-1.775148</td>
      <td>1.000000</td>
      <td>NIPAL3</td>
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
      <th>p_val</th>
      <th>dpi</th>
      <th>pos_iso_1</th>
      <th>pos_iso_2</th>
      <th>pos_iso_1_dpi</th>
      <th>pos_iso_2_dpi</th>
      <th>neg_iso_1</th>
      <th>neg_iso_2</th>
      <th>neg_iso_1_dpi</th>
      <th>neg_iso_2_dpi</th>
      <th>adj_p_val</th>
      <th>gname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>10</th>
      <td>ENSG00000003402</td>
      <td>7.338098e-15</td>
      <td>84.848486</td>
      <td>ENSG00000003402_3</td>
      <td>ENSG00000003402_7</td>
      <td>77.272728</td>
      <td>7.575758</td>
      <td>ENSG00000003402_1</td>
      <td>ENSG00000003402_5</td>
      <td>-31.717173</td>
      <td>-22.222223</td>
      <td>4.390279e-13</td>
      <td>CFLAR</td>
    </tr>
    <tr>
      <th>11</th>
      <td>ENSG00000003436</td>
      <td>3.812087e-06</td>
      <td>34.072281</td>
      <td>ENSG00000003436_1</td>
      <td>ENSG00000003436_3</td>
      <td>28.073771</td>
      <td>4.564083</td>
      <td>ENSG00000003436_4</td>
      <td>NaN</td>
      <td>-34.072281</td>
      <td>NaN</td>
      <td>7.064168e-05</td>
      <td>TFPI</td>
    </tr>
    <tr>
      <th>16</th>
      <td>ENSG00000004487</td>
      <td>1.135407e-03</td>
      <td>75.714287</td>
      <td>ENSG00000004487_1</td>
      <td>NaN</td>
      <td>75.714287</td>
      <td>NaN</td>
      <td>ENSG00000004487_2</td>
      <td>NaN</td>
      <td>-75.714285</td>
      <td>NaN</td>
      <td>1.020405e-02</td>
      <td>KDM1A</td>
    </tr>
    <tr>
      <th>33</th>
      <td>ENSG00000006282</td>
      <td>4.544048e-03</td>
      <td>19.819595</td>
      <td>ENSG00000006282_2</td>
      <td>ENSG00000006282_1</td>
      <td>13.679245</td>
      <td>6.140350</td>
      <td>ENSG00000006282_3</td>
      <td>NaN</td>
      <td>-19.819595</td>
      <td>NaN</td>
      <td>3.236475e-02</td>
      <td>SPATA20</td>
    </tr>
    <tr>
      <th>43</th>
      <td>ENSG00000007376</td>
      <td>1.256379e-04</td>
      <td>27.898551</td>
      <td>ENSG00000007376_3</td>
      <td>ENSG00000007376_1</td>
      <td>23.454107</td>
      <td>4.444445</td>
      <td>ENSG00000007376_5</td>
      <td>ENSG00000007376_7</td>
      <td>-14.130434</td>
      <td>-7.512077</td>
      <td>1.525134e-03</td>
      <td>RPUSD1</td>
    </tr>
  </tbody>
</table>
</div>



For TESs, use `kind='tes'` as input to `die_genes_test()`, `die_kind='tes'` to `make_uns_key()`, and `kind='tes'` to `get_die_genes()`.


```python
# find genes that exhibit DIE for TESs between HFFc6 and HepG2
die_table, die_results = sg.die_gene_test(kind='tes',
                                          obs_col='cell_line',
                                          obs_conditions=['hepg2', 'hffc6'])
die_table.head()
```

    Testing for DIE for each gene: 100%|██████████| 14684/14684 [02:19<00:00, 105.51it/s]





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
      <th>p_val</th>
      <th>dpi</th>
      <th>pos_iso_1</th>
      <th>pos_iso_2</th>
      <th>pos_iso_1_dpi</th>
      <th>pos_iso_2_dpi</th>
      <th>neg_iso_1</th>
      <th>neg_iso_2</th>
      <th>neg_iso_1_dpi</th>
      <th>neg_iso_2_dpi</th>
      <th>adj_p_val</th>
      <th>gname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000000419</td>
      <td>1.000000</td>
      <td>0.096133</td>
      <td>ENSG00000000419_2</td>
      <td>NaN</td>
      <td>0.096133</td>
      <td>NaN</td>
      <td>ENSG00000000419_1</td>
      <td>NaN</td>
      <td>-0.096130</td>
      <td>NaN</td>
      <td>1.000000</td>
      <td>DPM1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000001461</td>
      <td>0.852914</td>
      <td>9.093739</td>
      <td>ENSG00000001461_3</td>
      <td>NaN</td>
      <td>9.093739</td>
      <td>NaN</td>
      <td>ENSG00000001461_4</td>
      <td>ENSG00000001461_1</td>
      <td>-5.543442</td>
      <td>-1.775148</td>
      <td>1.000000</td>
      <td>NIPAL3</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ENSG00000001630</td>
      <td>0.286866</td>
      <td>2.580168</td>
      <td>ENSG00000001630_2</td>
      <td>NaN</td>
      <td>2.580168</td>
      <td>NaN</td>
      <td>ENSG00000001630_1</td>
      <td>ENSG00000001630_3</td>
      <td>-1.549232</td>
      <td>-1.030928</td>
      <td>0.618077</td>
      <td>CYP51A1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ENSG00000002330</td>
      <td>0.184635</td>
      <td>12.817431</td>
      <td>ENSG00000002330_1</td>
      <td>ENSG00000002330_4</td>
      <td>11.868484</td>
      <td>0.948944</td>
      <td>ENSG00000002330_2</td>
      <td>ENSG00000002330_3</td>
      <td>-12.603584</td>
      <td>-0.213847</td>
      <td>0.480860</td>
      <td>BAD</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ENSG00000002549</td>
      <td>0.679148</td>
      <td>0.694543</td>
      <td>ENSG00000002549_2</td>
      <td>NaN</td>
      <td>0.694542</td>
      <td>NaN</td>
      <td>ENSG00000002549_1</td>
      <td>NaN</td>
      <td>-0.694543</td>
      <td>NaN</td>
      <td>0.926889</td>
      <td>LAP3</td>
    </tr>
  </tbody>
</table>
</div>




```python
# die_iso - TES level differential isoform expression test results
uns_key = swan.make_uns_key('die',
                            obs_col=obs_col,
                            obs_conditions=obs_conditions,
                            die_kind='tes')
test = sg.adata.uns[uns_key]
test.head(2)
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
      <th>gid</th>
      <th>p_val</th>
      <th>dpi</th>
      <th>pos_iso_1</th>
      <th>pos_iso_2</th>
      <th>pos_iso_1_dpi</th>
      <th>pos_iso_2_dpi</th>
      <th>neg_iso_1</th>
      <th>neg_iso_2</th>
      <th>neg_iso_1_dpi</th>
      <th>neg_iso_2_dpi</th>
      <th>adj_p_val</th>
      <th>gname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000000419</td>
      <td>1.000000</td>
      <td>0.096133</td>
      <td>ENSG00000000419_2</td>
      <td>NaN</td>
      <td>0.096133</td>
      <td>NaN</td>
      <td>ENSG00000000419_1</td>
      <td>NaN</td>
      <td>-0.096130</td>
      <td>NaN</td>
      <td>1.0</td>
      <td>DPM1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000001461</td>
      <td>0.852914</td>
      <td>9.093739</td>
      <td>ENSG00000001461_3</td>
      <td>NaN</td>
      <td>9.093739</td>
      <td>NaN</td>
      <td>ENSG00000001461_4</td>
      <td>ENSG00000001461_1</td>
      <td>-5.543442</td>
      <td>-1.775148</td>
      <td>1.0</td>
      <td>NIPAL3</td>
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
      <th>p_val</th>
      <th>dpi</th>
      <th>pos_iso_1</th>
      <th>pos_iso_2</th>
      <th>pos_iso_1_dpi</th>
      <th>pos_iso_2_dpi</th>
      <th>neg_iso_1</th>
      <th>neg_iso_2</th>
      <th>neg_iso_1_dpi</th>
      <th>neg_iso_2_dpi</th>
      <th>adj_p_val</th>
      <th>gname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>9</th>
      <td>ENSG00000003402</td>
      <td>7.338098e-15</td>
      <td>84.848486</td>
      <td>ENSG00000003402_5</td>
      <td>ENSG00000003402_6</td>
      <td>77.272728</td>
      <td>7.575758</td>
      <td>ENSG00000003402_9</td>
      <td>ENSG00000003402_7</td>
      <td>-31.717173</td>
      <td>-22.222223</td>
      <td>5.296535e-13</td>
      <td>CFLAR</td>
    </tr>
    <tr>
      <th>10</th>
      <td>ENSG00000003436</td>
      <td>2.772096e-07</td>
      <td>32.637854</td>
      <td>ENSG00000003436_1</td>
      <td>ENSG00000003436_4</td>
      <td>22.082711</td>
      <td>10.555142</td>
      <td>ENSG00000003436_2</td>
      <td>ENSG00000003436_5</td>
      <td>-30.435921</td>
      <td>-1.818182</td>
      <td>7.003007e-06</td>
      <td>TFPI</td>
    </tr>
    <tr>
      <th>14</th>
      <td>ENSG00000004487</td>
      <td>1.135407e-03</td>
      <td>75.714287</td>
      <td>ENSG00000004487_2</td>
      <td>NaN</td>
      <td>75.714287</td>
      <td>NaN</td>
      <td>ENSG00000004487_1</td>
      <td>NaN</td>
      <td>-75.714285</td>
      <td>NaN</td>
      <td>1.141621e-02</td>
      <td>KDM1A</td>
    </tr>
    <tr>
      <th>32</th>
      <td>ENSG00000006282</td>
      <td>6.858641e-03</td>
      <td>19.819595</td>
      <td>ENSG00000006282_2</td>
      <td>NaN</td>
      <td>19.819595</td>
      <td>NaN</td>
      <td>ENSG00000006282_1</td>
      <td>NaN</td>
      <td>-19.819595</td>
      <td>NaN</td>
      <td>4.915014e-02</td>
      <td>SPATA20</td>
    </tr>
    <tr>
      <th>60</th>
      <td>ENSG00000010278</td>
      <td>3.191221e-22</td>
      <td>39.024387</td>
      <td>ENSG00000010278_2</td>
      <td>NaN</td>
      <td>39.024387</td>
      <td>NaN</td>
      <td>ENSG00000010278_3</td>
      <td>ENSG00000010278_1</td>
      <td>-38.605976</td>
      <td>-0.418410</td>
      <td>3.486194e-20</td>
      <td>CD9</td>
    </tr>
  </tbody>
</table>
</div>



Finally, for transcriptomes generated with [Cerberus](https://github.com/mortazavilab/cerberus), Swan automatically tracks intron chain usage, and you can perform intron chain switching analysis with `kind='ic'` as input to the `die_gene_test()` function.


```python
sg_brain = swan.read('swan_modelad.p')

# find genes that exhibit DIE for ICs between genotypes
die_table, die_results = sg_brain.die_gene_test(kind='ic',
                                                obs_col='genotype',
                                                obs_conditions=['b6n', '5xfad'])
die_table.head()
```

    Read in graph from swan_modelad.p





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
      <th>p_val</th>
      <th>dpi</th>
      <th>pos_iso_1</th>
      <th>pos_iso_2</th>
      <th>pos_iso_1_dpi</th>
      <th>pos_iso_2_dpi</th>
      <th>neg_iso_1</th>
      <th>neg_iso_2</th>
      <th>neg_iso_1_dpi</th>
      <th>neg_iso_2_dpi</th>
      <th>adj_p_val</th>
      <th>gname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENCODEMG000055991</td>
      <td>4.100695e-01</td>
      <td>6.071428</td>
      <td>ENCODEMG000055991_2</td>
      <td>ENCODEMG000055991_3</td>
      <td>3.571428</td>
      <td>2.500000</td>
      <td>ENCODEMG000055991_1</td>
      <td>NaN</td>
      <td>-6.071426</td>
      <td>NaN</td>
      <td>0.707992</td>
      <td>ENCODEMG000055991</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENCODEMG000055998</td>
      <td>6.190162e-01</td>
      <td>9.722223</td>
      <td>ENCODEMG000055998_2</td>
      <td>NaN</td>
      <td>9.722223</td>
      <td>NaN</td>
      <td>ENCODEMG000055998_1</td>
      <td>NaN</td>
      <td>-9.722218</td>
      <td>NaN</td>
      <td>0.836200</td>
      <td>ENCODEMG000055998</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ENCODEMG000056718</td>
      <td>8.081718e-07</td>
      <td>6.289159</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>2.329770</td>
      <td>1.953835</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-3.339438</td>
      <td>-2.949721</td>
      <td>0.000020</td>
      <td>ENCODEMG000056718</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ENCODEMG000056804</td>
      <td>2.107047e-03</td>
      <td>27.863049</td>
      <td>ENCODEMG000056804_1</td>
      <td>NaN</td>
      <td>27.863047</td>
      <td>NaN</td>
      <td>ENCODEMG000056804_2</td>
      <td>NaN</td>
      <td>-27.863049</td>
      <td>NaN</td>
      <td>0.021510</td>
      <td>ENCODEMG000056804</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ENCODEMG000063411</td>
      <td>7.699701e-01</td>
      <td>12.500000</td>
      <td>ENCODEMG000063411_1</td>
      <td>NaN</td>
      <td>12.500000</td>
      <td>NaN</td>
      <td>ENCODEMG000063411_2</td>
      <td>NaN</td>
      <td>-12.500000</td>
      <td>NaN</td>
      <td>0.919808</td>
      <td>ENCODEMG000063411</td>
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
    </tr>
    <tr>
      <th>dataset</th>
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
    </tr>
    <tr>
      <th>hepg2_2</th>
      <td>hepg2</td>
      <td>2</td>
      <td>hepg2_2</td>
      <td>848447.0</td>
    </tr>
    <tr>
      <th>hffc6_1</th>
      <td>hffc6</td>
      <td>1</td>
      <td>hffc6_1</td>
      <td>761493.0</td>
    </tr>
    <tr>
      <th>hffc6_2</th>
      <td>hffc6</td>
      <td>2</td>
      <td>hffc6_2</td>
      <td>787967.0</td>
    </tr>
    <tr>
      <th>hffc6_3</th>
      <td>hffc6</td>
      <td>3</td>
      <td>hffc6_3</td>
      <td>614921.0</td>
    </tr>
  </tbody>
</table>
</div>



Let's ignore for a moment the fact that the `dataset` column does effectively capture both replicate as well as cell line metadata. This may not always be the case with more complex datasets. Swan has a function to concatenate columns together and add them as an additional column to the metadata tables. Use the following code to generate a new column that concatenates as many preexisting metadata columns as you wish:


```python
col_name = sg.add_multi_groupby(['cell_line', 'replicate'])

print(col_name)
print(sg.adata.obs.head())
```

    cell_line_replicate
            cell_line replicate  dataset  total_counts cell_line_replicate
    dataset                                                               
    hepg2_1     hepg2         1  hepg2_1      499647.0             hepg2_1
    hepg2_2     hepg2         2  hepg2_2      848447.0             hepg2_2
    hffc6_1     hffc6         1  hffc6_1      761493.0             hffc6_1
    hffc6_2     hffc6         2  hffc6_2      787967.0             hffc6_2
    hffc6_3     hffc6         3  hffc6_3      614921.0             hffc6_3


The added column in `col_name` can then be used as the `obs_col` input to `die_gene_test()` as follows:


```python
obs_col = col_name
obs_conditions = ['hffc6_3', 'hepg2_1']


die_table, die_results = sg.die_gene_test(kind='iso',
                                          obs_col=obs_col,
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

    Testing each novel edge for exon skipping:   0%|          | 0/855 [00:00<?, ?it/s]

    Analyzing 855 intronic edges for ES


    Testing each novel edge for exon skipping: 100%|██████████| 855/855 [1:26:06<00:00,  6.13s/it]

    Found 529 novel es events in 149 transcripts.



```python
es_df.head()
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
      <th>gid</th>
      <th>tid</th>
      <th>edge_id</th>
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



To obtain a list of genes containing novel intron retention events, run the following code:


```python
# returns a DataFrame of genes, transcripts, and specific edges in
# the SwanGraph with novel intron retaining events
ir_df = sg.find_ir_genes(verbose=True)
```


      0%|          | 0/1186 [00:00<?, ?it/s][A
    Testing each novel edge for intron retention:   0%|          | 0/1186 [00:00<?, ?it/s][A

    Analyzing 1186 exonic edges for IR
    Testing each novel edge for intron retention: 100%|██████████| 1186/1186 [1:59:16<00:00,  5.97s/it][A

    Found 35 novel ir events in 27 transcripts.


You can pass gene IDs from `es_df` into `gen_report()` or `plot_graph()` to visualize where these exon-skipping events are in gene reports or gene summary graphs respectively.
