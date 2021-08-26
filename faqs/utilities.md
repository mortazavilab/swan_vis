# Utilities

Swan now comes with several utilities that can be used fo compute and output various metrics using data in the SwanGraph.

## Table of contents

* [Calculating TPM values](utilities.md#calculating-tpm-values)
* [Calculating pi values](utilities.md#calculating-pi-values)
* [Obtaining edge abundance information](utilities.md#obtaining-edge-abundance-information)
* [Obtaining TSS/TES abundance information](utilities.md#obtaining-tss-tes-abundance-information)


We'll be using the same SwanGraph as the rest of the tutorial pages to demonstrate these utilities. Load it using the following code:


```python
import swan_vis as swan

# code to download this data is in the Getting started tutorial
sg = swan.read('../tutorials/data/swan.p')
```

    Read in graph from ../tutorials/data/swan.p


##  <a name="calc_tpm"></a>Calculating TPM values

Swan allows for users to calculate the TPM of their data using various groupby metrics using the `calc_tpm()` function. You can use this to calculate TPM of any of the AnnData SwanGraph objects (`SwanGraph.adata` for transcripts, `SwanGraph.tss_adata` for TSSs, `SwanGraph.tes_adata` for TESs, and `SwanGraph.edge_adata` for edges; see the [Data structure](data_structure.md) FAQ page for more information on these tables.

First, we'll calculate the TPM for each transcript in each dataset:


```python
df = swan.calc_tpm(sg.adata)
df.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>tid</th>
      <th>ENST00000000233.9</th>
      <th>ENST00000000412.7</th>
      <th>ENST00000000442.10</th>
      <th>ENST00000001008.5</th>
      <th>ENST00000001146.6</th>
      <th>ENST00000002125.8</th>
      <th>ENST00000002165.10</th>
      <th>ENST00000002501.10</th>
      <th>ENST00000002596.5</th>
      <th>ENST00000002829.7</th>
      <th>...</th>
      <th>TALONT000482711</th>
      <th>TALONT000482903</th>
      <th>TALONT000483195</th>
      <th>TALONT000483284</th>
      <th>TALONT000483315</th>
      <th>TALONT000483322</th>
      <th>TALONT000483327</th>
      <th>TALONT000483978</th>
      <th>TALONT000484004</th>
      <th>TALONT000484796</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>hepg2_1</th>
      <td>196.138474</td>
      <td>86.060760</td>
      <td>8.005652</td>
      <td>46.032497</td>
      <td>0.0</td>
      <td>16.011305</td>
      <td>258.182281</td>
      <td>60.042389</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.000000</td>
      <td>4.002826</td>
      <td>2.001413</td>
      <td>12.008478</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>4.002826</td>
      <td>14.009891</td>
      <td>8.005652</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>hepg2_2</th>
      <td>243.975174</td>
      <td>77.789185</td>
      <td>7.071744</td>
      <td>61.288448</td>
      <td>0.0</td>
      <td>12.964864</td>
      <td>380.695557</td>
      <td>64.824318</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>1.178624</td>
      <td>14.143488</td>
      <td>4.714496</td>
      <td>7.071744</td>
      <td>2.357248</td>
      <td>8.250368</td>
      <td>2.357248</td>
      <td>11.786240</td>
      <td>10.607616</td>
      <td>1.178624</td>
    </tr>
    <tr>
      <th>hffc6_1</th>
      <td>131.320969</td>
      <td>194.355042</td>
      <td>0.000000</td>
      <td>107.683197</td>
      <td>0.0</td>
      <td>6.566049</td>
      <td>278.400452</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>6.566049</td>
      <td>13.132097</td>
      <td>9.192468</td>
      <td>1.313210</td>
      <td>6.566049</td>
      <td>9.192468</td>
      <td>6.566049</td>
      <td>0.000000</td>
      <td>15.758516</td>
      <td>1.313210</td>
    </tr>
    <tr>
      <th>hffc6_2</th>
      <td>137.061584</td>
      <td>242.395935</td>
      <td>0.000000</td>
      <td>124.370689</td>
      <td>0.0</td>
      <td>8.883621</td>
      <td>219.552338</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>15.229064</td>
      <td>10.152709</td>
      <td>6.345443</td>
      <td>8.883621</td>
      <td>1.269089</td>
      <td>10.152709</td>
      <td>15.229064</td>
      <td>0.000000</td>
      <td>16.498154</td>
      <td>8.883621</td>
    </tr>
    <tr>
      <th>hffc6_3</th>
      <td>147.986496</td>
      <td>273.205841</td>
      <td>3.252450</td>
      <td>172.379868</td>
      <td>0.0</td>
      <td>9.757351</td>
      <td>200.025696</td>
      <td>1.626225</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>14.636026</td>
      <td>11.383576</td>
      <td>8.131125</td>
      <td>8.131125</td>
      <td>11.383576</td>
      <td>11.383576</td>
      <td>6.504900</td>
      <td>0.000000</td>
      <td>24.393377</td>
      <td>9.757351</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 208306 columns</p>
</div>



We can swap out the first argument with the different AnnData structures in the SwanGraph. For instance, say we want to calculate the TPM of each TSS:


```python
df = swan.calc_tpm(sg.tss_adata)
df.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>tss_id</th>
      <th>ENSG00000000003.14_1</th>
      <th>ENSG00000000003.14_2</th>
      <th>ENSG00000000003.14_3</th>
      <th>ENSG00000000003.14_4</th>
      <th>ENSG00000000005.5_1</th>
      <th>ENSG00000000005.5_2</th>
      <th>ENSG00000000419.12_1</th>
      <th>ENSG00000000419.12_2</th>
      <th>ENSG00000000457.13_1</th>
      <th>ENSG00000000457.13_2</th>
      <th>...</th>
      <th>TALONG000085596_1</th>
      <th>TALONG000085799_1</th>
      <th>TALONG000085978_1</th>
      <th>TALONG000086022_1</th>
      <th>TALONG000086057_1</th>
      <th>TALONG000086218_1</th>
      <th>TALONG000086443_1</th>
      <th>TALONG000086539_1</th>
      <th>TALONG000086553_1</th>
      <th>TALONG000086766_1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>hepg2_1</th>
      <td>0.0</td>
      <td>232.163910</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>54.038151</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>60.042389</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>6.004239</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>8.005652</td>
    </tr>
    <tr>
      <th>hepg2_2</th>
      <td>0.0</td>
      <td>276.976654</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>103.718910</td>
      <td>0.0</td>
      <td>2.357248</td>
      <td>...</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>95.468544</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>31.822847</td>
      <td>0.000000</td>
      <td>1.178624</td>
      <td>10.607616</td>
    </tr>
    <tr>
      <th>hffc6_1</th>
      <td>0.0</td>
      <td>45.962341</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>101.117149</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>...</td>
      <td>2.626419</td>
      <td>6.566049</td>
      <td>9.192468</td>
      <td>0.000000</td>
      <td>7.879258</td>
      <td>11.818888</td>
      <td>233.751328</td>
      <td>9.192468</td>
      <td>6.566049</td>
      <td>15.758516</td>
    </tr>
    <tr>
      <th>hffc6_2</th>
      <td>0.0</td>
      <td>53.301723</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>85.028938</td>
      <td>0.0</td>
      <td>1.269089</td>
      <td>...</td>
      <td>6.345443</td>
      <td>1.269089</td>
      <td>12.690886</td>
      <td>0.000000</td>
      <td>8.883621</td>
      <td>20.305418</td>
      <td>119.294334</td>
      <td>12.690886</td>
      <td>2.538177</td>
      <td>16.498154</td>
    </tr>
    <tr>
      <th>hffc6_3</th>
      <td>0.0</td>
      <td>68.301460</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>89.442383</td>
      <td>0.0</td>
      <td>1.626225</td>
      <td>...</td>
      <td>8.131125</td>
      <td>8.131125</td>
      <td>11.383576</td>
      <td>0.000000</td>
      <td>11.383576</td>
      <td>17.888477</td>
      <td>134.976685</td>
      <td>27.645828</td>
      <td>8.131125</td>
      <td>24.393377</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 130176 columns</p>
</div>



And finally, we can use an alternative metadata column to compute TPM on. For instance, we can use the `cell_line` column:


```python
df = swan.calc_tpm(sg.adata, obs_col='cell_line')
df.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>tid</th>
      <th>ENST00000000233.9</th>
      <th>ENST00000000412.7</th>
      <th>ENST00000000442.10</th>
      <th>ENST00000001008.5</th>
      <th>ENST00000001146.6</th>
      <th>ENST00000002125.8</th>
      <th>ENST00000002165.10</th>
      <th>ENST00000002501.10</th>
      <th>ENST00000002596.5</th>
      <th>ENST00000002829.7</th>
      <th>...</th>
      <th>TALONT000482711</th>
      <th>TALONT000482903</th>
      <th>TALONT000483195</th>
      <th>TALONT000483284</th>
      <th>TALONT000483315</th>
      <th>TALONT000483322</th>
      <th>TALONT000483327</th>
      <th>TALONT000483978</th>
      <th>TALONT000484004</th>
      <th>TALONT000484796</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>hepg2</th>
      <td>226.245346</td>
      <td>80.854897</td>
      <td>7.417881</td>
      <td>55.634102</td>
      <td>0.0</td>
      <td>14.093972</td>
      <td>335.288177</td>
      <td>63.051983</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.741788</td>
      <td>10.385033</td>
      <td>3.70894</td>
      <td>8.901457</td>
      <td>1.483576</td>
      <td>5.192516</td>
      <td>2.967152</td>
      <td>12.610396</td>
      <td>9.643245</td>
      <td>0.741788</td>
    </tr>
    <tr>
      <th>hffc6</th>
      <td>138.145737</td>
      <td>234.247116</td>
      <td>0.924052</td>
      <td>132.139404</td>
      <td>0.0</td>
      <td>8.316465</td>
      <td>234.709137</td>
      <td>0.462026</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>12.012672</td>
      <td>11.550647</td>
      <td>7.85444</td>
      <td>6.006336</td>
      <td>6.006336</td>
      <td>10.164569</td>
      <td>9.702543</td>
      <td>0.000000</td>
      <td>18.481035</td>
      <td>6.468362</td>
    </tr>
  </tbody>
</table>
<p>2 rows × 208306 columns</p>
</div>



##  <a name="calc_pi"></a>Calculating pi values

You can use the `calc_pi()` function to calculate percent isoform use (pi) per gene in nearly the exact same way that you can use `calc_tpm()`: you can run it on either the transcript, edge, TSS, or TES level, and you can choose the metadata variable to groupby. The only difference is that for `calc_pi()` you must also provide an additional DataFrame object as the second argument that tells Swan what gene each entry comes from. Below the corresponding DataFrame that must be provided is listed for each AnnData:

| AnnData | DataFrame |
| ------- | --------- |
| `SwanGraph.adata` | `SwanGraph.t_df` |
| `SwanGraph.tss_adata` | `SwanGraph.tss_adata.var` |
| `SwanGraph.tes_adata` | `SwanGraph.tes_adata.var` |

First, we'll calculate the pi value for each transcript in each dataset:


```python
df, sums = swan.calc_pi(sg.adata, sg.t_df)
df.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>tid</th>
      <th>ENST00000000233.9</th>
      <th>ENST00000000412.7</th>
      <th>ENST00000000442.10</th>
      <th>ENST00000001008.5</th>
      <th>ENST00000001146.6</th>
      <th>ENST00000002125.8</th>
      <th>ENST00000002165.10</th>
      <th>ENST00000002501.10</th>
      <th>ENST00000002596.5</th>
      <th>ENST00000002829.7</th>
      <th>...</th>
      <th>TALONT000482711</th>
      <th>TALONT000482903</th>
      <th>TALONT000483195</th>
      <th>TALONT000483284</th>
      <th>TALONT000483315</th>
      <th>TALONT000483322</th>
      <th>TALONT000483327</th>
      <th>TALONT000483978</th>
      <th>TALONT000484004</th>
      <th>TALONT000484796</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>hepg2_1</th>
      <td>100.000000</td>
      <td>100.0</td>
      <td>100.000000</td>
      <td>100.0</td>
      <td>0.0</td>
      <td>100.000000</td>
      <td>100.0</td>
      <td>93.750000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.000000</td>
      <td>1.904762</td>
      <td>6.666667</td>
      <td>13.043478</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>1.333333</td>
      <td>100.0</td>
      <td>100.0</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>hepg2_2</th>
      <td>99.519226</td>
      <td>100.0</td>
      <td>60.000004</td>
      <td>100.0</td>
      <td>0.0</td>
      <td>100.000000</td>
      <td>100.0</td>
      <td>80.882355</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>5.263158</td>
      <td>3.225806</td>
      <td>13.793103</td>
      <td>8.695652</td>
      <td>0.884956</td>
      <td>3.097345</td>
      <td>0.884956</td>
      <td>100.0</td>
      <td>100.0</td>
      <td>2.380952</td>
    </tr>
    <tr>
      <th>hffc6_1</th>
      <td>98.039215</td>
      <td>100.0</td>
      <td>0.000000</td>
      <td>100.0</td>
      <td>0.0</td>
      <td>100.000000</td>
      <td>100.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>2.604167</td>
      <td>2.092050</td>
      <td>16.279070</td>
      <td>1.428571</td>
      <td>0.854701</td>
      <td>1.196581</td>
      <td>0.854701</td>
      <td>0.0</td>
      <td>100.0</td>
      <td>1.886792</td>
    </tr>
    <tr>
      <th>hffc6_2</th>
      <td>99.082573</td>
      <td>100.0</td>
      <td>0.000000</td>
      <td>100.0</td>
      <td>0.0</td>
      <td>77.777779</td>
      <td>100.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>4.285715</td>
      <td>2.144772</td>
      <td>11.627908</td>
      <td>14.893617</td>
      <td>0.166667</td>
      <td>1.333333</td>
      <td>2.000000</td>
      <td>0.0</td>
      <td>100.0</td>
      <td>9.859155</td>
    </tr>
    <tr>
      <th>hffc6_3</th>
      <td>100.000000</td>
      <td>100.0</td>
      <td>100.000000</td>
      <td>100.0</td>
      <td>0.0</td>
      <td>85.714287</td>
      <td>100.0</td>
      <td>100.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>4.326923</td>
      <td>2.536232</td>
      <td>15.151516</td>
      <td>10.638298</td>
      <td>1.711491</td>
      <td>1.711491</td>
      <td>0.977995</td>
      <td>0.0</td>
      <td>100.0</td>
      <td>13.636364</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 208306 columns</p>
</div>



As a note, the `calc_pi()` function outputs not only a table of pi values but of counts per isoform per condition, which is used as an intermediate during DIE testing. To avoid recalculation, it is output here.


```python
sums.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ENST00000000233.9</th>
      <th>ENST00000000412.7</th>
      <th>ENST00000000442.10</th>
      <th>ENST00000001008.5</th>
      <th>ENST00000001146.6</th>
      <th>ENST00000002125.8</th>
      <th>ENST00000002165.10</th>
      <th>ENST00000002501.10</th>
      <th>ENST00000002596.5</th>
      <th>ENST00000002829.7</th>
      <th>...</th>
      <th>TALONT000482711</th>
      <th>TALONT000482903</th>
      <th>TALONT000483195</th>
      <th>TALONT000483284</th>
      <th>TALONT000483315</th>
      <th>TALONT000483322</th>
      <th>TALONT000483327</th>
      <th>TALONT000483978</th>
      <th>TALONT000484004</th>
      <th>TALONT000484796</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>hepg2_1</th>
      <td>98.0</td>
      <td>43.0</td>
      <td>4.0</td>
      <td>23.0</td>
      <td>0.0</td>
      <td>8.0</td>
      <td>129.0</td>
      <td>30.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>1.0</td>
      <td>6.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>7.0</td>
      <td>4.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>hepg2_2</th>
      <td>207.0</td>
      <td>66.0</td>
      <td>6.0</td>
      <td>52.0</td>
      <td>0.0</td>
      <td>11.0</td>
      <td>323.0</td>
      <td>55.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>1.0</td>
      <td>12.0</td>
      <td>4.0</td>
      <td>6.0</td>
      <td>2.0</td>
      <td>7.0</td>
      <td>2.0</td>
      <td>10.0</td>
      <td>9.0</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>hffc6_1</th>
      <td>100.0</td>
      <td>148.0</td>
      <td>0.0</td>
      <td>82.0</td>
      <td>0.0</td>
      <td>5.0</td>
      <td>212.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>5.0</td>
      <td>10.0</td>
      <td>7.0</td>
      <td>1.0</td>
      <td>5.0</td>
      <td>7.0</td>
      <td>5.0</td>
      <td>0.0</td>
      <td>12.0</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>hffc6_2</th>
      <td>108.0</td>
      <td>191.0</td>
      <td>0.0</td>
      <td>98.0</td>
      <td>0.0</td>
      <td>7.0</td>
      <td>173.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>12.0</td>
      <td>8.0</td>
      <td>5.0</td>
      <td>7.0</td>
      <td>1.0</td>
      <td>8.0</td>
      <td>12.0</td>
      <td>0.0</td>
      <td>13.0</td>
      <td>7.0</td>
    </tr>
    <tr>
      <th>hffc6_3</th>
      <td>91.0</td>
      <td>168.0</td>
      <td>2.0</td>
      <td>106.0</td>
      <td>0.0</td>
      <td>6.0</td>
      <td>123.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>9.0</td>
      <td>7.0</td>
      <td>5.0</td>
      <td>5.0</td>
      <td>7.0</td>
      <td>7.0</td>
      <td>4.0</td>
      <td>0.0</td>
      <td>15.0</td>
      <td>6.0</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 208306 columns</p>
</div>



We can also calculate the pi value for the TSSs and TESs in each dataset:


```python
df, sums = swan.calc_pi(sg.tss_adata, sg.tss_adata.var)
print(df.head())
print()

df, sums = swan.calc_pi(sg.tes_adata, sg.tes_adata.var)
print(df.head())
print()
```

    tss_id   ENSG00000000003.14_1  ENSG00000000003.14_2  ENSG00000000003.14_3  \
    hepg2_1                   0.0                 100.0                   0.0   
    hepg2_2                   0.0                 100.0                   0.0   
    hffc6_1                   0.0                 100.0                   0.0   
    hffc6_2                   0.0                 100.0                   0.0   
    hffc6_3                   0.0                 100.0                   0.0   

    tss_id   ENSG00000000003.14_4  ENSG00000000005.5_1  ENSG00000000005.5_2  \
    hepg2_1                   0.0                  0.0                  0.0   
    hepg2_2                   0.0                  0.0                  0.0   
    hffc6_1                   0.0                  0.0                  0.0   
    hffc6_2                   0.0                  0.0                  0.0   
    hffc6_3                   0.0                  0.0                  0.0   

    tss_id   ENSG00000000419.12_1  ENSG00000000419.12_2  ENSG00000000457.13_1  \
    hepg2_1                   0.0                 100.0                   0.0   
    hepg2_2                   0.0                 100.0                   0.0   
    hffc6_1                   0.0                 100.0                   0.0   
    hffc6_2                   0.0                 100.0                   0.0   
    hffc6_3                   0.0                 100.0                   0.0   

    tss_id   ENSG00000000457.13_2  ...  TALONG000085596_1  TALONG000085799_1  \
    hepg2_1                   0.0  ...                0.0                0.0   
    hepg2_2                 100.0  ...                0.0                0.0   
    hffc6_1                   0.0  ...              100.0              100.0   
    hffc6_2                 100.0  ...              100.0              100.0   
    hffc6_3                 100.0  ...              100.0              100.0   

    tss_id   TALONG000085978_1  TALONG000086022_1  TALONG000086057_1  \
    hepg2_1                0.0              100.0                0.0   
    hepg2_2                0.0              100.0                0.0   
    hffc6_1              100.0                0.0              100.0   
    hffc6_2              100.0                0.0              100.0   
    hffc6_3              100.0                0.0              100.0   

    tss_id   TALONG000086218_1  TALONG000086443_1  TALONG000086539_1  \
    hepg2_1                0.0              100.0                0.0   
    hepg2_2                0.0              100.0                0.0   
    hffc6_1              100.0              100.0              100.0   
    hffc6_2              100.0              100.0              100.0   
    hffc6_3              100.0              100.0              100.0   

    tss_id   TALONG000086553_1  TALONG000086766_1  
    hepg2_1                0.0              100.0  
    hepg2_2              100.0              100.0  
    hffc6_1              100.0              100.0  
    hffc6_2              100.0              100.0  
    hffc6_3              100.0              100.0  

    [5 rows x 130176 columns]

    tes_id   ENSG00000000003.14_1  ENSG00000000003.14_2  ENSG00000000003.14_3  \
    hepg2_1                   0.0                 100.0                   0.0   
    hepg2_2                   0.0                 100.0                   0.0   
    hffc6_1                   0.0                 100.0                   0.0   
    hffc6_2                   0.0                 100.0                   0.0   
    hffc6_3                   0.0                 100.0                   0.0   

    tes_id   ENSG00000000003.14_4  ENSG00000000003.14_5  ENSG00000000005.5_1  \
    hepg2_1                   0.0                   0.0                  0.0   
    hepg2_2                   0.0                   0.0                  0.0   
    hffc6_1                   0.0                   0.0                  0.0   
    hffc6_2                   0.0                   0.0                  0.0   
    hffc6_3                   0.0                   0.0                  0.0   

    tes_id   ENSG00000000005.5_2  ENSG00000000419.12_1  ENSG00000000419.12_2  \
    hepg2_1                  0.0             92.592590                   0.0   
    hepg2_2                  0.0             98.863640                   0.0   
    hffc6_1                  0.0             98.701302                   0.0   
    hffc6_2                  0.0             95.522385                   0.0   
    hffc6_3                  0.0             98.181824                   0.0   

    tes_id   ENSG00000000419.12_3  ...  TALONG000085596_1  TALONG000085799_1  \
    hepg2_1              7.407407  ...                0.0                0.0   
    hepg2_2              1.136364  ...                0.0                0.0   
    hffc6_1              1.298701  ...              100.0              100.0   
    hffc6_2              4.477612  ...              100.0              100.0   
    hffc6_3              1.818182  ...              100.0              100.0   

    tes_id   TALONG000085978_1  TALONG000086022_1  TALONG000086057_1  \
    hepg2_1                0.0              100.0                0.0   
    hepg2_2                0.0              100.0                0.0   
    hffc6_1              100.0                0.0              100.0   
    hffc6_2              100.0                0.0              100.0   
    hffc6_3              100.0                0.0              100.0   

    tes_id   TALONG000086218_1  TALONG000086443_1  TALONG000086539_1  \
    hepg2_1                0.0              100.0                0.0   
    hepg2_2                0.0              100.0                0.0   
    hffc6_1              100.0              100.0              100.0   
    hffc6_2              100.0              100.0              100.0   
    hffc6_3              100.0              100.0              100.0   

    tes_id   TALONG000086553_1  TALONG000086766_1  
    hepg2_1                0.0              100.0  
    hepg2_2              100.0              100.0  
    hffc6_1              100.0              100.0  
    hffc6_2              100.0              100.0  
    hffc6_3              100.0              100.0  

    [5 rows x 187454 columns]



And we can also choose to calculate pi values using a different metadata column, here shown on the `cell_line` column:


```python
df, sums = swan.calc_pi(sg.adata, sg.t_df, obs_col='cell_line')
df.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>tid</th>
      <th>ENST00000000233.9</th>
      <th>ENST00000000412.7</th>
      <th>ENST00000000442.10</th>
      <th>ENST00000001008.5</th>
      <th>ENST00000001146.6</th>
      <th>ENST00000002125.8</th>
      <th>ENST00000002165.10</th>
      <th>ENST00000002501.10</th>
      <th>ENST00000002596.5</th>
      <th>ENST00000002829.7</th>
      <th>...</th>
      <th>TALONT000482711</th>
      <th>TALONT000482903</th>
      <th>TALONT000483195</th>
      <th>TALONT000483284</th>
      <th>TALONT000483315</th>
      <th>TALONT000483322</th>
      <th>TALONT000483327</th>
      <th>TALONT000483978</th>
      <th>TALONT000484004</th>
      <th>TALONT000484796</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>hepg2</th>
      <td>99.673203</td>
      <td>100.0</td>
      <td>71.428574</td>
      <td>100.0</td>
      <td>0.0</td>
      <td>100.000000</td>
      <td>100.0</td>
      <td>85.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>4.545455</td>
      <td>2.935011</td>
      <td>11.363637</td>
      <td>10.434782</td>
      <td>0.531915</td>
      <td>1.861702</td>
      <td>1.06383</td>
      <td>100.0</td>
      <td>100.0</td>
      <td>1.785714</td>
    </tr>
    <tr>
      <th>hffc6</th>
      <td>99.006622</td>
      <td>100.0</td>
      <td>100.000000</td>
      <td>100.0</td>
      <td>0.0</td>
      <td>85.714287</td>
      <td>100.0</td>
      <td>100.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>3.823529</td>
      <td>2.218279</td>
      <td>14.285715</td>
      <td>7.926829</td>
      <td>0.815558</td>
      <td>1.380176</td>
      <td>1.31744</td>
      <td>0.0</td>
      <td>100.0</td>
      <td>8.333334</td>
    </tr>
  </tbody>
</table>
<p>2 rows × 208306 columns</p>
</div>



##  <a name="edge_ab"></a>Obtaining edge abundance information

In case you're interested in doing outside analyses on the level (For instance, using [intron counting](https://pubmed.ncbi.nlm.nih.gov/23172860/) to assess alternative splicing), Swan provides a tool to output a DataFrame with edge abundance on the dataset level.

If we just want to get access to the edge abundance DataFrame, just use the `get_edge_abundance()` function:


```python
df = sg.get_edge_abundance()
df.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>strand</th>
      <th>edge_type</th>
      <th>annotation</th>
      <th>chrom</th>
      <th>start</th>
      <th>stop</th>
      <th>hepg2_1</th>
      <th>hepg2_2</th>
      <th>hffc6_1</th>
      <th>hffc6_2</th>
      <th>hffc6_3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>+</td>
      <td>exon</td>
      <td>True</td>
      <td>chr1</td>
      <td>11869</td>
      <td>12227</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>+</td>
      <td>exon</td>
      <td>True</td>
      <td>chr1</td>
      <td>12010</td>
      <td>12057</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>+</td>
      <td>intron</td>
      <td>True</td>
      <td>chr1</td>
      <td>12057</td>
      <td>12179</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>+</td>
      <td>exon</td>
      <td>True</td>
      <td>chr1</td>
      <td>12179</td>
      <td>12227</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>+</td>
      <td>intron</td>
      <td>True</td>
      <td>chr1</td>
      <td>12227</td>
      <td>12613</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
</div>



You can also specify if you want the data to be output in raw counts (`kind='counts'`) or TPM (`kind='tpm`). By default, this function returns counts. Here's an example with TPM:


```python
df = sg.get_edge_abundance(kind='tpm')
df.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>strand</th>
      <th>edge_type</th>
      <th>annotation</th>
      <th>chrom</th>
      <th>start</th>
      <th>stop</th>
      <th>hepg2_1</th>
      <th>hepg2_2</th>
      <th>hffc6_1</th>
      <th>hffc6_2</th>
      <th>hffc6_3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>+</td>
      <td>exon</td>
      <td>True</td>
      <td>chr1</td>
      <td>11869</td>
      <td>12227</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>+</td>
      <td>exon</td>
      <td>True</td>
      <td>chr1</td>
      <td>12010</td>
      <td>12057</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>+</td>
      <td>intron</td>
      <td>True</td>
      <td>chr1</td>
      <td>12057</td>
      <td>12179</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>+</td>
      <td>exon</td>
      <td>True</td>
      <td>chr1</td>
      <td>12179</td>
      <td>12227</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>+</td>
      <td>intron</td>
      <td>True</td>
      <td>chr1</td>
      <td>12227</td>
      <td>12613</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
</div>



And finally, if you wish, you can provide the function with a `prefix` value which will indicate that you want the output DataFrame to be saved in TSV form.


```python
df = sg.get_edge_abundance(kind='tpm', prefix='test')
df.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>strand</th>
      <th>edge_type</th>
      <th>annotation</th>
      <th>chrom</th>
      <th>start</th>
      <th>stop</th>
      <th>hepg2_1</th>
      <th>hepg2_2</th>
      <th>hffc6_1</th>
      <th>hffc6_2</th>
      <th>hffc6_3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>+</td>
      <td>exon</td>
      <td>True</td>
      <td>chr1</td>
      <td>11869</td>
      <td>12227</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>+</td>
      <td>exon</td>
      <td>True</td>
      <td>chr1</td>
      <td>12010</td>
      <td>12057</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>+</td>
      <td>intron</td>
      <td>True</td>
      <td>chr1</td>
      <td>12057</td>
      <td>12179</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>+</td>
      <td>exon</td>
      <td>True</td>
      <td>chr1</td>
      <td>12179</td>
      <td>12227</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>+</td>
      <td>intron</td>
      <td>True</td>
      <td>chr1</td>
      <td>12227</td>
      <td>12613</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
</div>



The results will be saved in '{prefix}_edge_abundance.tsv'.

##  <a name="end_ab"></a>Obtaining TSS/TES abundance information

Similarly, if you wish to do analysis involving your TSS or TES data, you can also output these using the `get_tss_abundance()` and `get_tes_abundance()` functions respectively. These have identical options to `get_edge_abundance()` so they can either output counts or TPM and optionally save to an output file.

First, let's output the TSS TPM to a file:


```python
df = sg.get_tss_abundance(kind='tpm', prefix='test')
df.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>tss_id</th>
      <th>gid</th>
      <th>gname</th>
      <th>vertex_id</th>
      <th>tss_name</th>
      <th>chrom</th>
      <th>coord</th>
      <th>hepg2_1</th>
      <th>hepg2_2</th>
      <th>hffc6_1</th>
      <th>hffc6_2</th>
      <th>hffc6_3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000000003.14_1</td>
      <td>ENSG00000000003.14</td>
      <td>TSPAN6</td>
      <td>926111</td>
      <td>TSPAN6_1</td>
      <td>chrX</td>
      <td>100636191</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000000003.14_2</td>
      <td>ENSG00000000003.14</td>
      <td>TSPAN6</td>
      <td>926112</td>
      <td>TSPAN6_2</td>
      <td>chrX</td>
      <td>100636608</td>
      <td>232.16391</td>
      <td>276.976654</td>
      <td>45.962341</td>
      <td>53.301723</td>
      <td>68.30146</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ENSG00000000003.14_3</td>
      <td>ENSG00000000003.14</td>
      <td>TSPAN6</td>
      <td>926114</td>
      <td>TSPAN6_3</td>
      <td>chrX</td>
      <td>100636793</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ENSG00000000003.14_4</td>
      <td>ENSG00000000003.14</td>
      <td>TSPAN6</td>
      <td>926117</td>
      <td>TSPAN6_4</td>
      <td>chrX</td>
      <td>100639945</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ENSG00000000005.5_1</td>
      <td>ENSG00000000005.5</td>
      <td>TNMD</td>
      <td>926077</td>
      <td>TNMD_1</td>
      <td>chrX</td>
      <td>100585066</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
    </tr>
  </tbody>
</table>
</div>



Now we'll get the counts of each TES without saving to a file:


```python
df = sg.get_tes_abundance(kind='counts')
df.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>tes_id</th>
      <th>gid</th>
      <th>gname</th>
      <th>vertex_id</th>
      <th>tes_name</th>
      <th>chrom</th>
      <th>coord</th>
      <th>hepg2_1</th>
      <th>hepg2_2</th>
      <th>hffc6_1</th>
      <th>hffc6_2</th>
      <th>hffc6_3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000000003.14_1</td>
      <td>ENSG00000000003.14</td>
      <td>TSPAN6</td>
      <td>926092</td>
      <td>TSPAN6_1</td>
      <td>chrX</td>
      <td>100627109</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000000003.14_2</td>
      <td>ENSG00000000003.14</td>
      <td>TSPAN6</td>
      <td>926093</td>
      <td>TSPAN6_2</td>
      <td>chrX</td>
      <td>100628670</td>
      <td>116.0</td>
      <td>235.0</td>
      <td>35.0</td>
      <td>42.0</td>
      <td>42.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ENSG00000000003.14_3</td>
      <td>ENSG00000000003.14</td>
      <td>TSPAN6</td>
      <td>926097</td>
      <td>TSPAN6_3</td>
      <td>chrX</td>
      <td>100632063</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ENSG00000000003.14_4</td>
      <td>ENSG00000000003.14</td>
      <td>TSPAN6</td>
      <td>926100</td>
      <td>TSPAN6_4</td>
      <td>chrX</td>
      <td>100632541</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ENSG00000000003.14_5</td>
      <td>ENSG00000000003.14</td>
      <td>TSPAN6</td>
      <td>926103</td>
      <td>TSPAN6_5</td>
      <td>chrX</td>
      <td>100633442</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
</div>
