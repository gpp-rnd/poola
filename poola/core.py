# AUTOGENERATED! DO NOT EDIT! File to edit: 00_core.ipynb (unless otherwise specified).

__all__ = ['lognorm', 'lognorm_columns', 'filter_pdna', 'calculate_lfcs', 'get_condition', 'average_replicate_lfcs',
           'group_pseudogenes', 'average_gene_lfcs', 'get_roc_aucs']

# Cell
import numpy as np

def lognorm(reads):
    """
    Standardize read counts by calculating reads per million,
    adding a pseudo-count of one, and taking the log2

    reads: numpy or pandas array |
    returns: numpy or pandas array
    """
    reads_per_million = (reads/reads.sum())*(10**6)
    lognormed_reads = np.log2(reads_per_million + 1)
    return lognormed_reads


def lognorm_columns(reads_df, columns):
    """
    Calculate lognorms for specified columns

    reads_df: dataframe |
    columns: list |
    returns: lognorm dataframe
    """
    lognorm_df = reads_df.copy()
    lognorm_df[columns] = lognorm_df[columns].apply(lognorm)
    return lognorm_df

# Cell
def filter_pdna(lognorm_df, pdna_cols, z_low=-3, z_high=None):
    """
    Filter the lognorm dataframe based on the z_scored lognorms of pDNA

    lognorm_df: dataframe |
    pdna_cols: list |
    z_low: int or None, lower range of z-scores to filter |
    z_high: int or None, uppper range of z-scores to filter |
    returns: filtered dataframe
    """
    filtered_lognorms = lognorm_df.copy()
    z_scored_cols = []
    for pdna in pdna_cols:
        z_col = pdna + '_z'
        filtered_lognorms[z_col] = (filtered_lognorms[pdna] - filtered_lognorms[pdna].mean())/filtered_lognorms[pdna].std()
        z_scored_cols.append(z_col)
    if z_low:
        filtered_lognorms = filtered_lognorms[filtered_lognorms[z_scored_cols].min(axis = 1) > z_low] # allow for multiple pDNA columns
    if z_high:
        filtered_lognorms = filtered_lognorms[filtered_lognorms[z_scored_cols].max(axis = 1) < z_high]
    filtered_lognorms = filtered_lognorms.drop(z_scored_cols, axis=1)
    return filtered_lognorms

# Cell
def calculate_lfcs(lognorm_df, ref_col=None, target_cols=None, ref_map=None):
    """
    Calculate log-fold changes between reference column and target columns

    lognorm_df: dataframe |
    ref_col: str or None, if str then target_col must also be present |
    target_col: list or None |
    ref_map: dict or None, key-value pairs correspond to target-reference conditions |
    returns: dataframe with log-fold changes
    """
    lfc_df = lognorm_df.copy()
    if (ref_col is not None) & (target_cols is not None):
        lfc_df[target_cols] = lfc_df[target_cols].sub(lfc_df[ref_col], axis=0)
        lfc_df = lfc_df.drop(ref_col, axis=1)
    elif ref_map is not None:
        for target_col, ref_col in ref_map.items():
            # use lognorms in case columns double as target and ref
            lfc_df[target_col] = lognorm_df[target_col] - lognorm_df[ref_col]
        for ref_col in set(ref_map.values()):
            if ref_col not in ref_map.keys(): # not a target condition as well
                lfc_df = lfc_df.drop(ref_col, axis=1)
    else:
        raise ValueError('Either ref_col and target_cols or ref_map must be present')
    return lfc_df

# Cell
from pandas.api.types import is_numeric_dtype

def get_condition(condition_name, sep, condition_indices):
    """Split replicate name from condition name"""
    split_condition = condition_name.split(sep)
    relevant_condition_elements = [split_condition[i] for i in condition_indices]
    condition = sep.join(relevant_condition_elements)
    return condition

def average_replicate_lfcs(lfcs, guide_col, condition_indices, sep='_', lfc_cols=None,
                           condition_name='condition', lfc_name='avg_lfc'):
    """
    Average log-fold changes of sgRNAs across replicates

    lfcs: dataframe |
    guide_col: str, sgrna column name |
    condition_indices: list of int, specifies which elements to use
        for conditions after separating column names with sep |
    sep: str, separator in column names |
    lfc_cols: list or None, lfc column(s) to melt. If None use all columns that are not guide_col |
    condition_name: str, name of condition columns |
    lfc_name: str, name of new column with log-fold changes |
    returns: dataframe of average lfcs
    """
    if lfc_cols is None:
        if not (lfcs.drop(guide_col,axis=1)
                .apply(is_numeric_dtype, axis=0).all()):
            raise ValueError('If lfc_cols are not supplied then all columns except guide_col must be numeric')
    long_lfcs = (lfcs.melt(id_vars=guide_col, value_vars=lfc_cols,
                           var_name=condition_name, value_name=lfc_name)
                 .reset_index())
    conditions = long_lfcs[condition_name].unique()
    condition_map = {cond: get_condition(cond, sep, condition_indices) for cond in conditions}
    long_lfcs[condition_name] = (long_lfcs[condition_name]
                                 .replace(condition_map))
    avg_lfcs = (long_lfcs.groupby([guide_col, condition_name])
                .agg(avg_lfc = (lfc_name, 'mean'),
                     n_obs = (lfc_name, 'count'))
                .reset_index())
    avg_lfcs = avg_lfcs.rename({'avg_lfc': lfc_name}, axis=1)
    return avg_lfcs

# Cell
def group_pseudogenes(annotations, pseudogene_size,
                      gene_col, control_regex, seed=7):
    """
    Remap annotations dataframe such that control genes are grouped into
    pseudo-genes

    annotations: dataframe |
    pseudogene_size: int |
    gene_col: str |
    control_regex: list of str, regular expressions to identify groups of
        pseudogenes |
    seed: int, random seed for reproducible outputs |
    returns: dataframe of annotations with controls grouped into pseudogenes |
    """
    remapped_annotations = annotations.copy()
    genes = remapped_annotations[gene_col]
    control_remap = {}
    for regex in control_regex:
        control_genes = genes[genes.str.contains(regex)].to_list()
        np.random.seed(seed)
        np.random.shuffle(control_genes) # shuffle mutates existing variable
        n_controls = len(control_genes)
        for i in range(n_controls):
            gene = control_genes[i]
            # Use modulo to get the right number of groupings
            gene_number = i % np.ceil(n_controls / pseudogene_size)
            control_remap[gene] = regex + '_' + str(int(gene_number))
    remapped_annotations[gene_col] = remapped_annotations[gene_col].replace(control_remap)
    return remapped_annotations

# Cell
from pandas.api.types import is_list_like

def average_gene_lfcs(lfcs, annotations, gene_col, condition_col='condition',
                      lfc_col='avg_lfc', merge_on=None, lfc_merge_on=None,
                      annotation_merge_on=None, controls_to_z=None):
    """
    Average log-fold changes of sgRNAs across genes

    lfcs: dataframe |
    annotations: dataframe, mapping between sgRNAs and genes |
    gene_col: str, column which uniquely identifies a gene |
    condition_col: str, column which uniquely identifies experimental conditions |
    lfc_col: str, name of value to average |
    merge_on: str or None, name of sgRNA column to merge on. Must be present in both
        lfc and annotation dataframes if supplied,
        otherwise supply unique merge column to each |
    lfc_merge_on: str or None, name of sgRNA column to merge on |
    annotation_merge_on: str, name of sgRNA column to merge on. Present in
        annotation dataframe if supplied |
    controls_to_z: str or list of str or None, if supplied specifies
        control genes for z-scoring log-fold changes
        if string then interpreted as a regex
        otherwise interpreted as a list of gene names found in gene_col |
    returns: dataframe, lfcs and (optionally) z-scores averaged by gene
    """
    annotated_lfcs = lfcs.merge(annotations, how='inner', on=merge_on, left_on=lfc_merge_on, right_on=annotation_merge_on)
    gene_lfcs = (annotated_lfcs
                 .groupby([condition_col, gene_col])
                 .agg(avg_lfc = (lfc_col, 'mean'),
                      n_obs = (lfc_col, 'count'))
                 .reset_index())
    gene_lfcs = gene_lfcs.rename({'avg_lfc': lfc_col}, axis=1)
    if controls_to_z != None:
        if isinstance(controls_to_z, str):
            control_lfcs = annotated_lfcs[annotated_lfcs[gene_col].str.contains(controls_to_z)]
        elif is_list_like(controls_to_z):
            control_lfcs = annotated_lfcs[annotated_lfcs[gene_col].isin(controls_to_z)]
        else:
            raise ValueError('Must Supply a string or list-like object to controls_to_z')
        control_stats = (control_lfcs.groupby(condition_col)
                         .agg(ctl_mean = (lfc_col, 'mean'),
                              ctl_sd = (lfc_col, 'std'))
                         .reset_index())
        gene_lfcs = (gene_lfcs.merge(control_stats, how='inner', on=condition_col))
        gene_lfcs[lfc_col + '_z-score'] = (gene_lfcs[lfc_col] - gene_lfcs['ctl_mean'])/(gene_lfcs['ctl_sd']/np.sqrt(gene_lfcs['n_obs']))
    return gene_lfcs

# Cell
from sklearn.metrics import roc_auc_score

def get_roc_aucs(lfcs, tp_genes, fp_genes, gene_col, score_col='avg_lfc', group_col='condition'):
    """
    Calculate the ROC-AUC between true positive and false positive gene sets

    lfcs: dataframe |
    tp_genes: list-like, true positive genes |
    fp_genes: list-like, false positive genes |
    gene_col: str |
    score_col: str, column to use for ranking genes from smallest to largest |
    group_col: list, columns to use for grouping genes |
    returns: datafrme of roc_aucs
    """
    roc_df = lfcs.copy()
    roc_df = roc_df[roc_df[gene_col].isin(tp_genes) |
                    roc_df[gene_col].isin(fp_genes)]
    roc_df['tp'] = roc_df[gene_col].isin(tp_genes)
    roc_aucs = (roc_df.groupby(group_col)
                .apply(lambda df: roc_auc_score(df['tp'], -df[score_col]))
                .reset_index(name='ROC-AUC'))
    return roc_aucs