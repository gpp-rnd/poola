# poola
> Python package for pooled screen analysis


## Install

Install from github:

`pip install git+git://github.com/gpp-rnd/poola.git#egg=poola`

## How to use

Additional packages required for this tutorial can be install using `pip install -r requirements.txt`

```python
from poola import core as pool
import pandas as pd
import seaborn as sns
import gpplot
import matplotlib.pyplot as plt
import requests
```

To demonstrate the functionality of this module we'll use read counts from [Sanson et al. 2018](https://doi.org/10.1038/s41467-018-07901-8).

```python
supp_reads = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-07901-8/MediaObjects/41467_2018_7901_MOESM4_ESM.xlsx'
read_counts = pd.read_excel(supp_reads,
                            sheet_name = 'A375_orig_tracr raw reads', 
                            header = None,
                            skiprows = 3, 
                            names = ['sgRNA Sequence', 'pDNA', 'A375_RepA', 'A375_RepB'])
guide_annotations = pd.read_excel(supp_reads,
                                  sheet_name='sgRNA annotations')
```

```python
lognorms = pool.lognorm_columns(reads_df=read_counts, columns=['pDNA', 'A375_RepA', 'A375_RepB'])
filtered_lognorms = pool.filter_pdna(lognorm_df=lognorms, pdna_cols=['pDNA'])
print('Filtered ' + str(lognorms.shape[0] - filtered_lognorms.shape[0]) + ' columns due to low pDNA abundance')
```

    Filtered 576 columns due to low pDNA abundance


```python
lfc_df = pool.calculate_lfcs(lognorm_df=filtered_lognorms, ref_col='pDNA', target_cols=['A375_RepA', 'A375_RepB'])
plt.subplots(figsize=(4,4))
gpplot.point_densityplot(data=lfc_df, x='A375_RepA', y='A375_RepB')
gpplot.add_correlation(data=lfc_df, x='A375_RepA', y='A375_RepB')
sns.despine()
```


![png](docs/images/output_8_0.png)


Since we see a strong correlation, we'll average the log-fold change of each sgRNA across replicates

```python
avg_replicate_lfc_df = pool.average_replicate_lfcs(lfcs=lfc_df, guide_col='sgRNA Sequence', condition_indices=[0])
```

Before combining sgRNAs at the gene level, it's sometimes helpful to group controls into pseudo-genes so they're easier to compare with target genes.

Using this reamapped annotations file, we'll average log-fold changes for each gene

```python
remapped_annotations = pool.group_pseudogenes(annotations=guide_annotations, pseudogene_size=4, 
                                              gene_col='Annotated Gene Symbol', 
                                              control_regex=['NO_CURRENT'])
gene_lfcs = pool.average_gene_lfcs(lfcs=avg_replicate_lfc_df, annotations=remapped_annotations, gene_col='Annotated Gene Symbol',
                                   merge_on='sgRNA Sequence', controls_to_z='NO_CURRENT')
```

Finally, to evaluate the quality this screen, we'll calculate the ROC-AUC between [essential](https://doi.org/10.1016/j.cell.2015.11.015) and [nonessential](https://doi.org/10.15252/msb.20145216) genes for each condition

```python
noness_file = "https://www.embopress.org/action/downloadSupplement?doi=10.15252%2Fmsb.20145216&file=msb145216-sup-0001-DatasetS1.xlsx"
noness_genes = (pd.read_excel(requests.get(noness_file).content, sheet_name='ReferenceSets', usecols=['Nonessential Genes (NE)'])
                .rename({'Nonessential Genes (NE)': 'gene'}, axis=1))
ess_file = 'http://tko.ccbr.utoronto.ca/Data/core-essential-genes-sym_HGNCID'
ess_genes = pd.read_table(ess_file, names=['gene', 'gene_id'])
```

    /Users/pdeweird/.local/share/virtualenvs/poola-tuJn2lJU/lib/python3.8/site-packages/openpyxl/worksheet/_reader.py:308: UserWarning: Unknown extension is not supported and will be removed
      warn(msg)


```python
roc_aucs = pool.get_roc_aucs(lfcs=gene_lfcs, tp_genes=ess_genes.gene, fp_genes=noness_genes.gene, gene_col='Annotated Gene Symbol')
assert roc_aucs['ROC-AUC'][0] > 0.9
print('ROC-AUC: ' + str(round(roc_aucs['ROC-AUC'].values[0], 3)))
```

    ROC-AUC: 0.976

