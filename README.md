# poola
> Python package for pooled screen analysis


## Install

Install from github for the latest development release:

`pip install git+git://github.com/gpp-rnd/poola.git#egg=poola`

Or install the most recent distribution from PyPi:

`pip install poola`

## How to use

Additional packages required for this tutorial can be install using `pip install -r requirements.txt`

```python
from poola import core as pool
import pandas as pd
import seaborn as sns
import gpplot
import matplotlib.pyplot as plt
```

To demonstrate the functionality of this module we'll use read counts from [Sanson et al. 2018](https://doi.org/10.1038/s41467-018-07901-8).

```python
supp_reads = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-07901-8/MediaObjects/41467_2018_7901_MOESM4_ESM.xlsx'
read_counts = pd.read_excel(supp_reads,
                            sheet_name = 'A375_orig_tracr raw reads', 
                            header = None,
                            skiprows = 3, 
                            names = ['sgRNA Sequence', 'pDNA', 'A375_RepA', 'A375_RepB'], 
                            engine='openpyxl')
guide_annotations = pd.read_excel(supp_reads,
                                  sheet_name='sgRNA annotations', 
                                  engine='openpyxl')
```

The input data has three columns with read counts and one column with sgRNA annotations

```python
read_counts.head()
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
      <th>sgRNA Sequence</th>
      <th>pDNA</th>
      <th>A375_RepA</th>
      <th>A375_RepB</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>AAAAAAAATCCGGACAATGG</td>
      <td>522</td>
      <td>729</td>
      <td>774</td>
    </tr>
    <tr>
      <th>1</th>
      <td>AAAAAAAGGATGGTGATCAA</td>
      <td>511</td>
      <td>1484</td>
      <td>1393</td>
    </tr>
    <tr>
      <th>2</th>
      <td>AAAAAAATGACATTACTGCA</td>
      <td>467</td>
      <td>375</td>
      <td>603</td>
    </tr>
    <tr>
      <th>3</th>
      <td>AAAAAAATGTCAGTCGAGTG</td>
      <td>200</td>
      <td>737</td>
      <td>506</td>
    </tr>
    <tr>
      <th>4</th>
      <td>AAAAAACACAAGCAAGACCG</td>
      <td>286</td>
      <td>672</td>
      <td>352</td>
    </tr>
  </tbody>
</table>
</div>



```python
lognorms = pool.lognorm_columns(reads_df=read_counts, columns=['pDNA', 'A375_RepA', 'A375_RepB'])
filtered_lognorms = pool.filter_pdna(lognorm_df=lognorms, pdna_cols=['pDNA'], z_low=-3)
print('Filtered ' + str(lognorms.shape[0] - filtered_lognorms.shape[0]) + ' columns due to low pDNA abundance')
```

    Filtered 576 columns due to low pDNA abundance


Note that the column names for the lognorms remain the same

```python
lognorms.head()
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
      <th>sgRNA Sequence</th>
      <th>pDNA</th>
      <th>A375_RepA</th>
      <th>A375_RepB</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>AAAAAAAATCCGGACAATGG</td>
      <td>4.192756</td>
      <td>3.373924</td>
      <td>3.521755</td>
    </tr>
    <tr>
      <th>1</th>
      <td>AAAAAAAGGATGGTGATCAA</td>
      <td>4.163726</td>
      <td>4.326828</td>
      <td>4.312620</td>
    </tr>
    <tr>
      <th>2</th>
      <td>AAAAAAATGACATTACTGCA</td>
      <td>4.041390</td>
      <td>2.540624</td>
      <td>3.196767</td>
    </tr>
    <tr>
      <th>3</th>
      <td>AAAAAAATGTCAGTCGAGTG</td>
      <td>2.930437</td>
      <td>3.388159</td>
      <td>2.973599</td>
    </tr>
    <tr>
      <th>4</th>
      <td>AAAAAACACAAGCAAGACCG</td>
      <td>3.388394</td>
      <td>3.268222</td>
      <td>2.528233</td>
    </tr>
  </tbody>
</table>
</div>



```python
lfc_df = pool.calculate_lfcs(lognorm_df=filtered_lognorms, ref_col='pDNA', target_cols=['A375_RepA', 'A375_RepB'])
```

We drop the pDNA column after calculating log-fold changes

```python
lfc_df.head()
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
      <th>sgRNA Sequence</th>
      <th>A375_RepA</th>
      <th>A375_RepB</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>AAAAAAAATCCGGACAATGG</td>
      <td>-0.818831</td>
      <td>-0.671000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>AAAAAAAGGATGGTGATCAA</td>
      <td>0.163102</td>
      <td>0.148894</td>
    </tr>
    <tr>
      <th>2</th>
      <td>AAAAAAATGACATTACTGCA</td>
      <td>-1.500766</td>
      <td>-0.844622</td>
    </tr>
    <tr>
      <th>3</th>
      <td>AAAAAAATGTCAGTCGAGTG</td>
      <td>0.457721</td>
      <td>0.043161</td>
    </tr>
    <tr>
      <th>4</th>
      <td>AAAAAACACAAGCAAGACCG</td>
      <td>-0.120172</td>
      <td>-0.860161</td>
    </tr>
  </tbody>
</table>
</div>



Since we only have two conditions it's easy to visualize replicates as a point densityplot using [gpplot](https://github.com/gpp-rnd/gpplot)

```python
plt.subplots(figsize=(4,4))
gpplot.point_densityplot(data=lfc_df, x='A375_RepA', y='A375_RepB')
gpplot.add_correlation(data=lfc_df, x='A375_RepA', y='A375_RepB')
sns.despine()
```


![png](docs/images/output_18_0.png)


Since we see a strong correlation, we'll average the log-fold change of each sgRNA across replicates

```python
avg_replicate_lfc_df = pool.average_replicate_lfcs(lfcs=lfc_df, guide_col='sgRNA Sequence', condition_indices=[0], sep='_')
```

After averaging log-fold changes our dataframe is melted, so the condition column specifies the experimental condition (A375 here) and the n_obs specifies the number of replicates

```python
avg_replicate_lfc_df.head()
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
      <th>sgRNA Sequence</th>
      <th>condition</th>
      <th>avg_lfc</th>
      <th>n_obs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>AAAAAAAATCCGGACAATGG</td>
      <td>A375</td>
      <td>-0.744916</td>
      <td>2</td>
    </tr>
    <tr>
      <th>1</th>
      <td>AAAAAAAGGATGGTGATCAA</td>
      <td>A375</td>
      <td>0.155998</td>
      <td>2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>AAAAAAATGACATTACTGCA</td>
      <td>A375</td>
      <td>-1.172694</td>
      <td>2</td>
    </tr>
    <tr>
      <th>3</th>
      <td>AAAAAAATGTCAGTCGAGTG</td>
      <td>A375</td>
      <td>0.250441</td>
      <td>2</td>
    </tr>
    <tr>
      <th>4</th>
      <td>AAAAAACACAAGCAAGACCG</td>
      <td>A375</td>
      <td>-0.490166</td>
      <td>2</td>
    </tr>
  </tbody>
</table>
</div>



Before combining sgRNAs at the gene level, it's sometimes helpful to group controls into pseudo-genes so they're easier to compare with target genes. Our annotation file maps from sgRNA sequences to gene symbols

```python
remapped_annotations = pool.group_pseudogenes(annotations=guide_annotations, pseudogene_size=4, 
                                              gene_col='Annotated Gene Symbol', 
                                              control_regex=['NO_CURRENT'])
remapped_annotations.head()
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
      <th>sgRNA Sequence</th>
      <th>Annotated Gene Symbol</th>
      <th>Annotated Gene ID</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>AAAAAAAATCCGGACAATGG</td>
      <td>SLC25A24</td>
      <td>29957</td>
    </tr>
    <tr>
      <th>1</th>
      <td>AAAAAAAGGATGGTGATCAA</td>
      <td>FASTKD3</td>
      <td>79072</td>
    </tr>
    <tr>
      <th>2</th>
      <td>AAAAAAATGACATTACTGCA</td>
      <td>BCAS2</td>
      <td>10286</td>
    </tr>
    <tr>
      <th>3</th>
      <td>AAAAAAATGTCAGTCGAGTG</td>
      <td>GPR18</td>
      <td>2841</td>
    </tr>
    <tr>
      <th>4</th>
      <td>AAAAAACACAAGCAAGACCG</td>
      <td>ZNF470</td>
      <td>388566</td>
    </tr>
  </tbody>
</table>
</div>



Using this remapped annotations file, we'll average log-fold changes for each gene, and calculate z-scores using the nonessential genes as controls. When `z_score_neg_ctls=True`, we can specify `z_score_neg_ctl_genes` as a list or regex to specify genes for our null distribution. Note that if `z_score_neg_ctl_genes=None` then all genes are used to generate the null

```python
nonessential_genes = (pd.read_table('https://raw.githubusercontent.com/gpp-rnd/genesets/master/human/non-essential-genes-Hart2014.txt', 
                                    names=['gene'])
                      .gene)
gene_lfcs = pool.average_gene_lfcs(lfcs=avg_replicate_lfc_df, annotations=remapped_annotations, gene_col='Annotated Gene Symbol',
                                   merge_on='sgRNA Sequence', z_score_neg_ctls=True, 
                                   z_score_neg_ctl_genes=nonessential_genes)
```

From our z-scores we calculate a p-value and FDR using the [Benjamini-Hochberg procedure](http://www.biostathandbook.com/multiplecomparisons.html)

```python
gene_lfcs.sort_values('z_scored_avg_lfc').head()
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
      <th>condition</th>
      <th>Annotated Gene Symbol</th>
      <th>avg_lfc</th>
      <th>n_obs</th>
      <th>z_scored_avg_lfc</th>
      <th>z_scored_avg_lfc_p_value</th>
      <th>z_scored_avg_lfc_fdr_bh</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>13317</th>
      <td>A375</td>
      <td>PSMG3</td>
      <td>-3.165679</td>
      <td>4</td>
      <td>-13.076520</td>
      <td>4.485080e-39</td>
      <td>8.684460e-35</td>
    </tr>
    <tr>
      <th>7375</th>
      <td>A375</td>
      <td>HSPA5</td>
      <td>-3.021085</td>
      <td>4</td>
      <td>-12.492789</td>
      <td>8.173727e-36</td>
      <td>7.913393e-32</td>
    </tr>
    <tr>
      <th>4927</th>
      <td>A375</td>
      <td>EIF6</td>
      <td>-3.009939</td>
      <td>4</td>
      <td>-12.447794</td>
      <td>1.437640e-35</td>
      <td>9.279010e-32</td>
    </tr>
    <tr>
      <th>14219</th>
      <td>A375</td>
      <td>RPL19</td>
      <td>-2.993390</td>
      <td>4</td>
      <td>-12.380987</td>
      <td>3.312447e-35</td>
      <td>1.603473e-31</td>
    </tr>
    <tr>
      <th>12774</th>
      <td>A375</td>
      <td>POLR2L</td>
      <td>-2.970395</td>
      <td>4</td>
      <td>-12.288152</td>
      <td>1.048764e-34</td>
      <td>4.061445e-31</td>
    </tr>
  </tbody>
</table>
</div>



Finally, to evaluate the quality this screen, we'll calculate the ROC-AUC between [essential](https://doi.org/10.1016/j.cell.2015.11.015) and [nonessential](https://doi.org/10.15252/msb.20145216) genes for each condition

```python
essential_genes = (pd.read_table('https://raw.githubusercontent.com/gpp-rnd/genesets/master/human/essential-genes-Hart2015.txt', 
                                 names=['gene'])
                   .gene)
roc_aucs = pool.get_roc_aucs(lfcs=gene_lfcs, tp_genes=essential_genes, fp_genes=nonessential_genes, 
                             gene_col='Annotated Gene Symbol', score_col='avg_lfc', group_col='condition')
print('ROC-AUC: ' + str(round(roc_aucs['ROC-AUC'].values[0], 3)))
```

    ROC-AUC: 0.976


Note that we can also use this function to calculate roc-aucs at the guide level

```python
annotated_guide_lfcs = lfc_df.merge(guide_annotations, how='inner', on='sgRNA Sequence')
roc_aucs = pool.get_roc_aucs(lfcs=annotated_guide_lfcs, tp_genes=essential_genes, fp_genes=nonessential_genes, gene_col='Annotated Gene Symbol',
                             conditions=['A375_RepA', 'A375_RepB'])
print('Rep A AUC: ' + str(round(roc_aucs.loc[roc_aucs.condition == 'A375_RepA', 'ROC-AUC'].values[0], 4)))
print('Rep B AUC: ' + str(round(roc_aucs.loc[roc_aucs.condition == 'A375_RepB', 'ROC-AUC'].values[0], 4)))
```

    Rep A AUC: 0.9185
    Rep B AUC: 0.9176

