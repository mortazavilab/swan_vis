# Analysis tools

Swan has several analysis options to use. 
* [Differential gene expression](#deg)
* [Differential transcript expression](#det)
* [Isoform switching](#is)
* [Exon skipping and intron retention](#es_ir)


```python
import swan_vis as swan

sg = swan.SwanGraph('data/swan.p')
```

## <a name="deg"></a>Differential gene expression tests

Differential gene expression testing in Swan is implemented via [diffxpy](https://github.com/theislab/diffxpy). To run the test, first partition the datasets that you have added to your SwanGraph into biological replicates. Then, use this grouping to run the differential expression test.


```python
dataset_groups = [['HepG2_1','HepG2_2'],['HFFc6_1','HFFc6_2','HFFc6_3']]
sg.de_gene_test(dataset_groups)
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
      <th>93003</th>
      <td>ENSG00000158874.11</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.798466</td>
      <td>7203.040430</td>
      <td>False</td>
      <td>8.236582e-07</td>
      <td>-9.798466</td>
      <td>0.597998</td>
      <td>-21.964623</td>
      <td>APOA2</td>
    </tr>
    <tr>
      <th>99184</th>
      <td>ENSG00000163631.16</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.657541</td>
      <td>6256.314453</td>
      <td>False</td>
      <td>2.962665e-07</td>
      <td>-9.657541</td>
      <td>0.583054</td>
      <td>-20.338973</td>
      <td>ALB</td>
    </tr>
    <tr>
      <th>79865</th>
      <td>ENSG00000145192.12</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.320415</td>
      <td>4466.045703</td>
      <td>False</td>
      <td>4.583478e-08</td>
      <td>-9.320415</td>
      <td>0.591606</td>
      <td>-20.613694</td>
      <td>AHSG</td>
    </tr>
    <tr>
      <th>4810</th>
      <td>ENSG00000026025.15</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>9.226718</td>
      <td>9419.225325</td>
      <td>False</td>
      <td>1.069652e-06</td>
      <td>9.226718</td>
      <td>0.569277</td>
      <td>-24.014993</td>
      <td>VIM</td>
    </tr>
    <tr>
      <th>139314</th>
      <td>ENSG00000197249.13</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.106863</td>
      <td>5051.564227</td>
      <td>False</td>
      <td>1.033368e-06</td>
      <td>-9.106863</td>
      <td>0.496072</td>
      <td>-20.785889</td>
      <td>SERPINA1</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>176368</th>
      <td>ENSG00000242968.1</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.000000</td>
      <td>1.000000</td>
      <td>False</td>
      <td>0.000000e+00</td>
      <td>0.000000</td>
      <td>0.912871</td>
      <td>0.000000</td>
      <td>AC096992.1</td>
    </tr>
    <tr>
      <th>176367</th>
      <td>ENSG00000242963.1</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.000000</td>
      <td>1.000000</td>
      <td>False</td>
      <td>0.000000e+00</td>
      <td>0.000000</td>
      <td>0.912871</td>
      <td>0.000000</td>
      <td>AC026336.1</td>
    </tr>
    <tr>
      <th>176366</th>
      <td>ENSG00000242960.1</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.000000</td>
      <td>1.000000</td>
      <td>False</td>
      <td>0.000000e+00</td>
      <td>0.000000</td>
      <td>0.912871</td>
      <td>0.000000</td>
      <td>FTH1P23</td>
    </tr>
    <tr>
      <th>176365</th>
      <td>ENSG00000242958.1</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.000000</td>
      <td>1.000000</td>
      <td>False</td>
      <td>0.000000e+00</td>
      <td>0.000000</td>
      <td>0.912871</td>
      <td>0.000000</td>
      <td>AC040975.1</td>
    </tr>
    <tr>
      <th>163669</th>
      <td>ENSG00000230383.1</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.000000</td>
      <td>1.000000</td>
      <td>False</td>
      <td>0.000000e+00</td>
      <td>0.000000</td>
      <td>0.912871</td>
      <td>0.000000</td>
      <td>AC009245.1</td>
    </tr>
  </tbody>
</table>
<p>58906 rows × 11 columns</p>
</div>



The results of this test are stored in `sg.deg_test` so they can be accessed later as follows. Test results will also be stored if `save_graph()` is run again so the user can easily load the results up.


```python
sg.deg_test.head()
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
      <th>93003</th>
      <td>ENSG00000158874.11</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.798466</td>
      <td>7203.040430</td>
      <td>False</td>
      <td>8.236582e-07</td>
      <td>-9.798466</td>
      <td>0.597998</td>
      <td>-21.964623</td>
      <td>APOA2</td>
    </tr>
    <tr>
      <th>99184</th>
      <td>ENSG00000163631.16</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.657541</td>
      <td>6256.314453</td>
      <td>False</td>
      <td>2.962665e-07</td>
      <td>-9.657541</td>
      <td>0.583054</td>
      <td>-20.338973</td>
      <td>ALB</td>
    </tr>
    <tr>
      <th>79865</th>
      <td>ENSG00000145192.12</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.320415</td>
      <td>4466.045703</td>
      <td>False</td>
      <td>4.583478e-08</td>
      <td>-9.320415</td>
      <td>0.591606</td>
      <td>-20.613694</td>
      <td>AHSG</td>
    </tr>
    <tr>
      <th>4810</th>
      <td>ENSG00000026025.15</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>9.226718</td>
      <td>9419.225325</td>
      <td>False</td>
      <td>1.069652e-06</td>
      <td>9.226718</td>
      <td>0.569277</td>
      <td>-24.014993</td>
      <td>VIM</td>
    </tr>
    <tr>
      <th>139314</th>
      <td>ENSG00000197249.13</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.106863</td>
      <td>5051.564227</td>
      <td>False</td>
      <td>1.033368e-06</td>
      <td>-9.106863</td>
      <td>0.496072</td>
      <td>-20.785889</td>
      <td>SERPINA1</td>
    </tr>
  </tbody>
</table>
</div>



Swan can also automatically subset the test summary table to pull out genes that pass a certain significance threshold. These genes can be directly passed into Swan's gene plotting functions, `gen_report()` or `plot_graph()`


```python
gene_ids, gene_summary = sg.get_de_genes(q=0.05)
```


```python
print(gene_ids[:5])
gene_summary.head()
```

    ['APOA2', 'ALB', 'AHSG', 'VIM', 'SERPINA1']


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
      <th>93003</th>
      <td>ENSG00000158874.11</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.798466</td>
      <td>7203.040430</td>
      <td>False</td>
      <td>8.236582e-07</td>
      <td>-9.798466</td>
      <td>0.597998</td>
      <td>-21.964623</td>
      <td>APOA2</td>
    </tr>
    <tr>
      <th>99184</th>
      <td>ENSG00000163631.16</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.657541</td>
      <td>6256.314453</td>
      <td>False</td>
      <td>2.962665e-07</td>
      <td>-9.657541</td>
      <td>0.583054</td>
      <td>-20.338973</td>
      <td>ALB</td>
    </tr>
    <tr>
      <th>79865</th>
      <td>ENSG00000145192.12</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.320415</td>
      <td>4466.045703</td>
      <td>False</td>
      <td>4.583478e-08</td>
      <td>-9.320415</td>
      <td>0.591606</td>
      <td>-20.613694</td>
      <td>AHSG</td>
    </tr>
    <tr>
      <th>4810</th>
      <td>ENSG00000026025.15</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>9.226718</td>
      <td>9419.225325</td>
      <td>False</td>
      <td>1.069652e-06</td>
      <td>9.226718</td>
      <td>0.569277</td>
      <td>-24.014993</td>
      <td>VIM</td>
    </tr>
    <tr>
      <th>139314</th>
      <td>ENSG00000197249.13</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.106863</td>
      <td>5051.564227</td>
      <td>False</td>
      <td>1.033368e-06</td>
      <td>-9.106863</td>
      <td>0.496072</td>
      <td>-20.785889</td>
      <td>SERPINA1</td>
    </tr>
  </tbody>
</table>
</div>



## <a name="det"></a>Differential transcript expression tests

Similarly, Swan can run tests to find differentially expressed transcript isoforms. The input and output to these functions are identical to that of the differential gene tests.


```python
dataset_groups = [['HepG2_1','HepG2_2'],['HFFc6_1','HFFc6_2','HFFc6_3']]
sg.de_transcript_test(dataset_groups);
sg.det_test.head()
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
      <th>21203</th>
      <td>ENST00000367990.7</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.641046</td>
      <td>6153.972461</td>
      <td>False</td>
      <td>2.680981e-07</td>
      <td>-9.641046</td>
      <td>0.599458</td>
      <td>-21.723931</td>
      <td>ENSG00000158874.11</td>
      <td>APOA2</td>
    </tr>
    <tr>
      <th>6792</th>
      <td>ENST00000295897.8</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.599193</td>
      <td>5901.748828</td>
      <td>False</td>
      <td>7.406762e-07</td>
      <td>-9.599193</td>
      <td>0.583689</td>
      <td>-20.330176</td>
      <td>ENSG00000163631.16</td>
      <td>ALB</td>
    </tr>
    <tr>
      <th>31598</th>
      <td>ENST00000393087.8</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.435849</td>
      <td>5012.441016</td>
      <td>False</td>
      <td>6.445624e-07</td>
      <td>-9.435849</td>
      <td>0.584451</td>
      <td>-20.120285</td>
      <td>ENSG00000197249.13</td>
      <td>SERPINA1</td>
    </tr>
    <tr>
      <th>39289</th>
      <td>ENST00000411641.6</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.236873</td>
      <td>4108.149609</td>
      <td>False</td>
      <td>4.956935e-08</td>
      <td>-9.236873</td>
      <td>0.592742</td>
      <td>-20.527872</td>
      <td>ENSG00000145192.12</td>
      <td>AHSG</td>
    </tr>
    <tr>
      <th>1044</th>
      <td>ENST00000224237.9</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>9.203737</td>
      <td>9205.245637</td>
      <td>False</td>
      <td>1.456279e-06</td>
      <td>9.203737</td>
      <td>0.569203</td>
      <td>-23.507290</td>
      <td>ENSG00000026025.15</td>
      <td>VIM</td>
    </tr>
  </tbody>
</table>
</div>



And Swan can subset the results for you based on a q-value significance threshold. The resultant transcript ids can then be passed into Swan's transcript plotting function, `plot_transcript_path()`.


```python
transcript_ids, transcript_summary = sg.get_de_transcripts(q=0.05)
```


```python
print(transcript_ids[:5])
transcript_summary.head()
```

    ['ENST00000367990.7', 'ENST00000295897.8', 'ENST00000393087.8', 'ENST00000411641.6', 'ENST00000224237.9']


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
      <th>21203</th>
      <td>ENST00000367990.7</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.641046</td>
      <td>6153.972461</td>
      <td>False</td>
      <td>2.680981e-07</td>
      <td>-9.641046</td>
      <td>0.599458</td>
      <td>-21.723931</td>
      <td>ENSG00000158874.11</td>
      <td>APOA2</td>
    </tr>
    <tr>
      <th>6792</th>
      <td>ENST00000295897.8</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.599193</td>
      <td>5901.748828</td>
      <td>False</td>
      <td>7.406762e-07</td>
      <td>-9.599193</td>
      <td>0.583689</td>
      <td>-20.330176</td>
      <td>ENSG00000163631.16</td>
      <td>ALB</td>
    </tr>
    <tr>
      <th>31598</th>
      <td>ENST00000393087.8</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.435849</td>
      <td>5012.441016</td>
      <td>False</td>
      <td>6.445624e-07</td>
      <td>-9.435849</td>
      <td>0.584451</td>
      <td>-20.120285</td>
      <td>ENSG00000197249.13</td>
      <td>SERPINA1</td>
    </tr>
    <tr>
      <th>39289</th>
      <td>ENST00000411641.6</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>-9.236873</td>
      <td>4108.149609</td>
      <td>False</td>
      <td>4.956935e-08</td>
      <td>-9.236873</td>
      <td>0.592742</td>
      <td>-20.527872</td>
      <td>ENSG00000145192.12</td>
      <td>AHSG</td>
    </tr>
    <tr>
      <th>1044</th>
      <td>ENST00000224237.9</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>9.203737</td>
      <td>9205.245637</td>
      <td>False</td>
      <td>1.456279e-06</td>
      <td>9.203737</td>
      <td>0.569203</td>
      <td>-23.507290</td>
      <td>ENSG00000026025.15</td>
      <td>VIM</td>
    </tr>
  </tbody>
</table>
</div>



## <a name="is"></a>Isoform switching

We wanted to include a module to conduct rudimentary isoform switching analysis as well. We define a gene that exhibits isoform switching for our purposes as a gene that is not differentially expressed that has transcript isoforms that are differentially expressed. We have provided code to detect such instances. To run it, `de_gene_test()` and `de_transcript_test()` must first be run.


```python
is_genes, is_table = sg.find_isoform_switching_genes(q=0.05)
```


```python
print(is_genes[:5])
is_table.head()
```

    ['ENSG00000197746.13', 'ENSG00000067225.17', 'ENSG00000117450.13', 'ENSG00000105254.11', 'ENSG00000177600.8']



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
      <th>0</th>
      <td>TALONT000283514</td>
      <td>5.099119e-10</td>
      <td>8.195811e-08</td>
      <td>5.538936</td>
      <td>153.044266</td>
      <td>False</td>
      <td>2.116293e-07</td>
      <td>5.538936</td>
      <td>0.891074</td>
      <td>-21.442096</td>
      <td>ENSG00000197746.13</td>
      <td>PSAP</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENST00000394936.7</td>
      <td>9.324034e-05</td>
      <td>5.520899e-03</td>
      <td>-0.970022</td>
      <td>278.640723</td>
      <td>False</td>
      <td>5.421814e-08</td>
      <td>-0.970022</td>
      <td>0.248244</td>
      <td>-28.092922</td>
      <td>ENSG00000197746.13</td>
      <td>PSAP</td>
    </tr>
    <tr>
      <th>2</th>
      <td>TALONT000316712</td>
      <td>2.490252e-09</td>
      <td>3.622447e-07</td>
      <td>5.038395</td>
      <td>92.933382</td>
      <td>False</td>
      <td>2.948318e-08</td>
      <td>5.038395</td>
      <td>0.845071</td>
      <td>-19.486628</td>
      <td>ENSG00000067225.17</td>
      <td>PKM</td>
    </tr>
    <tr>
      <th>3</th>
      <td>TALONT000375121</td>
      <td>2.220446e-16</td>
      <td>7.353454e-14</td>
      <td>-4.855414</td>
      <td>51.973529</td>
      <td>False</td>
      <td>1.056665e-08</td>
      <td>-4.855414</td>
      <td>0.592949</td>
      <td>-11.909066</td>
      <td>ENSG00000117450.13</td>
      <td>PRDX1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ENST00000585910.5</td>
      <td>1.398881e-13</td>
      <td>3.448465e-11</td>
      <td>-3.687907</td>
      <td>23.070891</td>
      <td>False</td>
      <td>4.299709e-01</td>
      <td>-3.687907</td>
      <td>0.498608</td>
      <td>-7.304680</td>
      <td>ENSG00000105254.11</td>
      <td>TBCB</td>
    </tr>
  </tbody>
</table>
</div>



## <a name="es_ir"></a>Exon skipping and intron retention

Swan can detect novel (unannotated) exon skipping and intron retention events. 

To obtain a list of genes containing novel exon skipping events, run the following code:


```python
es_genes = sg.find_es_genes()
print(es_genes[:5])
```

    Analyzing 893 intronic edges for ES
    Found 1021 novel es events from 285 genes.
    ['ENSG00000130706.12', 'ENSG00000111237.18', 'ENSG00000101363.12', 'ENSG00000163069.12', 'ENSG00000132677.12']


As usual, we can feed `es_genes` into `gen_report()` or individual gene ids from `es_genes` into `plot_graph()` to generate gene reports or gene summary graphs respectively.

To obtain a list of genes containing novel intron retention events, run the following code:


```python
ir_genes = sg.find_ir_genes()
print(ir_genes[:5])
```

    Analyzing 2185 exonic edges for IR
    Found 73 novel ir events from 47 genes.
    ['ENSG00000213719.8', 'ENSG00000135480.15', 'ENSG00000196421.8', 'ENSG00000141505.11', 'ENSG00000030582.17']


As usual, we can feed `ir_genes` into `gen_report()` or individual gene ids from `ir_genes` into `plot_graph()` to generate gene reports or gene summary graphs respectively.
