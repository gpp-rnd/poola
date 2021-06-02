# AUTOGENERATED! DO NOT EDIT! File to edit: 00_core.ipynb (unless otherwise specified).

__all__ = ['lognorm', 'lognorm_columns', 'filter_pdna', 'calculate_lfcs', 'get_condition', 'average_replicate_lfcs',
           'group_pseudogenes', 'get_annotated_lfcs', 'get_control_lfcs', 'get_neg_ctl_z_score', 'scale_lfcs',
           'annotate_guide_lfcs', 'aggregate_gene_lfcs', 'get_roc_aucs']

# Cell

import numpy as np
from pandas.api.types import is_numeric_dtype, is_list_like
from statsmodels.stats.multitest import multipletests
import scipy
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
import pandas as pd

# Cell

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
    z_high: int or None, upper range of z-scores to filter |
    returns: filtered dataframe
    """
    filtered_lognorms = lognorm_df.copy()
    z_scored_cols = []
    for pdna in pdna_cols:
        z_col = pdna + '_z'
        filtered_lognorms[z_col] = (filtered_lognorms[pdna] - filtered_lognorms[pdna].mean())/filtered_lognorms[pdna].std()
        z_scored_cols.append(z_col)
    if z_low is not None:
        filtered_lognorms = filtered_lognorms[filtered_lognorms[z_scored_cols].min(axis = 1) > z_low] # allow for multiple pDNA columns
    if z_high is not None:
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

def get_condition(condition_name, sep, condition_indices):
    """Split replicate name from condition name"""
    split_condition = condition_name.split(sep)
    relevant_condition_elements = [split_condition[i] for i in condition_indices]
    condition = sep.join(relevant_condition_elements)
    return condition

def average_replicate_lfcs(lfcs, guide_col, condition_indices=None, sep=None, condition_map=None,
                           lfc_cols=None, condition_name='condition', lfc_name='avg_lfc'):
    """Average log-fold changes of sgRNAs across replicates

    lfcs: dataframe |
    guide_col: str, sgrna column name |
    condition_indices: list of int or None, specifies which elements to use
        for conditions after separating column names with sep |
    sep: str or None, separator in column names |
    condition_map: dict or None, alternative to supplying condition_indices and sep. Keys are column
        names and values are condition names |
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
                           var_name=condition_name, value_name=lfc_name))
    if condition_map is None:
        # must supply sep and condition_indices
        if (sep is not None) and (condition_indices is not None):
            conditions = long_lfcs[condition_name].unique()
            condition_map = {cond: get_condition(cond, sep, condition_indices) for cond in conditions}
        elif ((sep is None) and (condition_indices is not None) or
              (sep is not None) and (condition_indices is None)):
            raise ValueError('Must supply sep AND condition_indices')
    if condition_map is not None:
        condition_map_df = pd.DataFrame({condition_name: condition_map.keys(),
                                         condition_name + '_new': condition_map.values()})
        long_lfcs = (long_lfcs.merge(condition_map_df, how='inner', on=condition_name)
                     .drop(condition_name, axis=1)
                     .rename({condition_name + '_new': condition_name}, axis=1))
    avg_lfcs = (long_lfcs.groupby([guide_col, condition_name])
                .agg(avg_lfc = (lfc_name, 'mean'),
                     n_obs = (lfc_name, 'count'))
                .reset_index())
    avg_lfcs = avg_lfcs.rename({'avg_lfc': lfc_name}, axis=1)
    return avg_lfcs

# Cell
def group_pseudogenes(annotations, pseudogene_size,
                      gene_col, control_regex, seed=7):
    """Remap annotations dataframe such that control genes are grouped into
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

def get_annotated_lfcs(lfcs, annotations, merge_on=None, lfc_merge_on=None, annotation_merge_on=None):
    """Merge lfcs and annotatsion

    lfcs: dataframe, log-fold changes
    annotations: dataframe, sgRNA annotations
    merge_on: str or None, name of sgRNA column to merge on. Must be present in both
        lfc and annotation dataframes if supplied,
        otherwise supply unique merge column to each |
    lfc_merge_on: str or None, name of sgRNA column to merge on |
    annotation_merge_on: str or None, name of sgRNA column to merge on
    """
    annotated_lfcs = lfcs.merge(annotations, how='inner', on=merge_on, left_on=lfc_merge_on, right_on=annotation_merge_on)
    return annotated_lfcs


def get_control_lfcs(annotated_lfcs, controls, gene_col):
    """Get lfcs from control populations

    annotated_lfcs: dataframe, merged lfcs and annotations |
    controls: str or list of str or None, if string then interpreted as a regex
        otherwise interpreted as a list of gene names found in gene_col |
    gene_col: str, column which uniquely identifies a gene |
    """
    if controls is None:
        control_lfcs = annotated_lfcs.copy()
    elif isinstance(controls, str):
        control_lfcs = annotated_lfcs[annotated_lfcs[gene_col].str.contains(controls)].copy()
    elif is_list_like(controls):
        control_lfcs = annotated_lfcs[annotated_lfcs[gene_col].isin(controls)].copy()
    else:
        raise ValueError('Must supply a string, list-like object, or None for controls')
    return control_lfcs


def get_neg_ctl_z_score(lfcs, negative_controls, gene_col, condition_col, lfc_col):
    """Z-score lfcs by negative control population

    lfcs: dataframe, merged lfcs and annotations |
    negative_controls: str or list of str or None, if string then interpreted as a regex
        otherwise interpreted as a list of gene names found in gene_col |
    gene_col: str, column which uniquely identifies a gene |
    condition_col: str, column which uniquely identifies experimental conditions |
    lfc_col: str, name of value to average |
    returns: dataframe, gene lfcs with column for z-score
    """
    # control lfcs are taken at the guide level
    control_lfcs = get_control_lfcs(lfcs, negative_controls, gene_col)
    control_stats = (control_lfcs.groupby(condition_col)
                     .agg(neg_ctl_mean = (lfc_col, 'mean'),
                          neg_ctl_sd = (lfc_col, 'std'))
                     .reset_index())
    lfcs = lfcs.copy()
    lfcs = (lfcs.merge(control_stats, how='inner', on=condition_col))
    lfcs['z_scored_' + lfc_col] = (lfcs[lfc_col] - lfcs['neg_ctl_mean'])/lfcs['neg_ctl_sd']
    lfcs = lfcs.drop(['neg_ctl_sd', 'neg_ctl_mean'], axis=1)
    return lfcs


def scale_lfcs(lfcs, negative_controls, positive_controls,
               gene_col, condition_col, lfc_col, pos_control_direction):
    """Scale lfcs between negative and positive control populations

    lfcs: dataframe, merged lfcs and annotations |
    negative_controls: str or list of str or None, if string then interpreted as a regex
        otherwise interpreted as a list of gene names found in gene_col |
    positive_controls: str or list of str of None, if string then interpreted as a regex
        otherwise interpreted as a list of gene names found in gene_col |
    gene_col: str, column which uniquely identifies a gene |
    condition_col: str, column which uniquely identifies experimental conditions |
    lfc_col: str, name of value to average |
    pos_control_direction: int, direction in which positive controls score |
    returns: dataframe, gene lfcs with column for z-score and z-score p-values
    """
    negative_control_lfcs = get_control_lfcs(lfcs, negative_controls, gene_col)
    positive_control_lfcs = get_control_lfcs(lfcs, positive_controls, gene_col)
    neg_median_lfcs = (negative_control_lfcs.groupby(condition_col)
                       .agg(neg_ctl_median = (lfc_col, 'median'))
                       .reset_index())
    pos_median_lfcs = (positive_control_lfcs.groupby(condition_col)
                       .agg(pos_ctl_median = (lfc_col, 'median'))
                       .reset_index())
    lfcs = lfcs.copy()
    lfcs = (lfcs
            .merge(neg_median_lfcs, how='inner', on=condition_col)
            .merge(pos_median_lfcs, how='inner', on=condition_col))
    lfcs['control_scaled_' + lfc_col] = (lfcs[lfc_col] - lfcs['neg_ctl_median'])/(lfcs['neg_ctl_median'] - lfcs['pos_ctl_median'])
    lfcs = lfcs.drop(['neg_ctl_median', 'pos_ctl_median'], axis=1)
    group_df_list = []
    for condition, group_df in lfcs.groupby(condition_col):
        pos_condition = positive_control_lfcs[positive_control_lfcs[condition_col] == condition].copy()
        pos_condition['positive'] = 1
        neg_condition = negative_control_lfcs[negative_control_lfcs[condition_col] == condition].copy()
        neg_condition['positive'] = 0
        condition_control_lfcs = pd.concat([pos_condition, neg_condition])
        model = LogisticRegression().fit(pos_control_direction * condition_control_lfcs[[lfc_col]], condition_control_lfcs['positive'])
        group_df['probability_positive_control'] = model.predict_proba(pos_control_direction * group_df[[lfc_col]])[:,1]
        group_df_list.append(group_df)
    lfcs = pd.concat(group_df_list)
    return lfcs

# Cell

def annotate_guide_lfcs(lfcs, annotations, gene_col, condition_col='condition', lfc_col='avg_lfc',
                             merge_on=None, lfc_merge_on=None, annotation_merge_on=None,
                             z_score_neg_ctls=False, z_score_neg_ctl_genes=None,
                             scale_ctls=False, scale_neg_ctl_genes=None, scale_pos_ctl_genes=None,
                             pos_control_direction=-1):
    """
    Join guide log-fold changes with annotations, and optionally calculate z-scores and scaled scores

    lfcs: dataframe |
    annotations: dataframe, mapping between sgRNAs and genes |
    gene_col: str, column which uniquely identifies a gene |
    condition_col: str, column which uniquely identifies experimental conditions |
    lfc_col: str, name of value to average |
    merge_on: str or None, name of sgRNA column to merge on. Must be present in both
        lfc and annotation dataframes if supplied,
        otherwise supply unique merge column to each |
    lfc_merge_on: str or None, name of sgRNA column to merge on |
    annotation_merge_on: str or None, name of sgRNA column to merge on. Present in
        annotation dataframe if supplied |
    z_score_neg_ctls: bool, z-score log-fold changes relative to negative controls |
    z_score_neg_ctl_genes: str, list of str, or None, if string then interpreted as a regex
        otherwise interpreted as a list of gene names found in gene_col. If None and z_score_neg_ctls is True,
        then uses all sgRNAs as controls |
    scale_ctls: bool, scale-lfcs by positive and negative controls |
    scale_neg_ctl_genes: str or list of str or None, if string then interpreted as a regex
        otherwise interpreted as a list of gene names found in gene_col |
    scale_pos_ctl_genes: str or list of str or None, if string then interpreted as a regex
        otherwise interpreted as a list of gene names found in gene_col |
    pos_control_direction: int, direction in which positive controls score |
    returns: dataframe, lfcs and (optionally) z-scores averaged by gene
    """
    annotated_lfcs = get_annotated_lfcs(lfcs, annotations, merge_on=merge_on, lfc_merge_on=lfc_merge_on,
                                        annotation_merge_on=annotation_merge_on)
    if z_score_neg_ctls:
        annotated_lfcs = get_neg_ctl_z_score(annotated_lfcs, z_score_neg_ctl_genes, gene_col, condition_col, lfc_col)
    if scale_ctls:
        annotated_lfcs = scale_lfcs(annotated_lfcs, scale_neg_ctl_genes, scale_pos_ctl_genes,
                                    gene_col, condition_col, lfc_col, pos_control_direction)
    return annotated_lfcs


# Cell

def aggregate_gene_lfcs(lfcs, gene_col, condition_col='condition', average_cols=None, zscore_cols=None):
    """
    Aggregate log-fold changes at the gene level

    lfcs: dataframe, scores with annotations |
    gene_col: str, column which uniquely identifies a gene |
    condition_col: str, column which uniquely identifies experimental conditions |
    average_cols: list of str, columns to average
    zscore_cols: list of str, columns to z-score --> will also caculate p-value and FDR for each column
    returns: dataframe, lfcs aggregated at the gene level
    """
    agg_mean = {col: 'mean' for col in average_cols}
    agg_sum = {col: 'sum' for col in zscore_cols}
    aggs = {condition_col: 'count', **agg_mean, **agg_sum}
    agg_lfcs = (lfcs
                .groupby([condition_col, gene_col])
                .agg(aggs)
                .rename({condition_col: 'n_guides'}, axis=1)
                .reset_index())
    for col in zscore_cols:
        agg_lfcs[col] = agg_lfcs[col]/np.sqrt(agg_lfcs['n_guides'])
        agg_lfcs[col + '_p_value'] = (agg_lfcs.groupby(condition_col)
                                      [col]
                                      .transform(lambda x: scipy.stats.norm.sf(abs(x))*2))
        agg_lfcs[col + '_fdr'] = (agg_lfcs.groupby(condition_col)
                                  [col + '_p_value']
                                  .transform(lambda x: multipletests(x, method='fdr_bh')[1]))
    return agg_lfcs

# Cell

def get_roc_aucs(lfcs, tp_genes, fp_genes, gene_col, score_col=None, condition_col=None, condition_list=None,
                 pos_control_direction=-1):
    """
    Calculate the ROC-AUC between true positive and false positive gene sets.
    Must specificy score_col (and group_col optionally) or conditions

    lfcs: dataframe |
    tp_genes: list-like or str, true positive genes |
    fp_genes: list-like or str, false positive genes |
    gene_col: str |
    score_col: str, column for ranking genes from smallest to largest.
    condition_col: str or list, columns to use for grouping genes.
    condition_list: list, columns for which to calculate ROC-AUCs
    pos_control_direction: int, direction in which positive controls score |
    returns: dataframe of roc_aucs
    """
    if condition_list is not None:
        score_col = 'lfc'
        condition_col = 'condition'
        roc_df = lfcs.melt(id_vars=gene_col, value_vars=condition_list, var_name=condition_col, value_name=score_col)
    elif score_col is not None:
        roc_df = lfcs.copy()
    else:
        raise ValueError('conditions or score_col must be specified')
    if is_list_like(tp_genes):
        tp_bool = roc_df[gene_col].isin(tp_genes)
    elif type(tp_genes) is str:
        tp_bool = roc_df[gene_col].str.contains(tp_genes)
    else:
        raise ValueError('tp_genes must be list-like or str')
    if is_list_like(fp_genes):
        fp_bool = roc_df[gene_col].isin(fp_genes)
    elif type(fp_genes) is str:
        fp_bool = roc_df[gene_col].str.contains(fp_genes)
    else:
        raise ValueError('fp_genes must be list-like or str')
    roc_df['tp'] = tp_bool
    filter_bool = (tp_bool | fp_bool)
    roc_df = roc_df[filter_bool].reset_index(drop=True)
    if condition_col is not None:
        tpr_fpr_df_list = []
        roc_auc_list = []
        for group, df in roc_df.groupby(condition_col):
            fpr, tpr, treshold = roc_curve(df['tp'], pos_control_direction * df[score_col])
            group_tpr_fpr_df = pd.DataFrame({'tpr': tpr, 'fpr': fpr, 'threshold': treshold})
            group_tpr_fpr_df[condition_col] = group
            tpr_fpr_df_list.append(group_tpr_fpr_df)
            roc_auc = auc(fpr, tpr)
            roc_auc_list.append({condition_col: group, 'ROC-AUC': roc_auc})
        roc_aucs = pd.DataFrame(roc_auc_list)
        tpr_fpr_df = (pd.concat(tpr_fpr_df_list).reset_index(drop=True))
    else:
        fpr, tpr, treshold = roc_curve(roc_df['tp'], pos_control_direction * roc_df[score_col])
        tpr_fpr_df = pd.DataFrame({'tpr': tpr, 'fpr': fpr, 'threshold': treshold})
        roc_aucs = auc(fpr, tpr)
    return roc_aucs, tpr_fpr_df