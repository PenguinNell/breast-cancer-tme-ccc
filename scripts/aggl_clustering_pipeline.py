from dataclasses import dataclass

import numpy as np
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import (
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
)


@dataclass
class ClusteringResult:
    y_pred_mapped_best: pd.Series
    df_cluster_neg: pd.DataFrame
    genefu_labels_for_non_target: pd.Series


def predict_clusters(X, linkage, metric):

    if linkage == 'ward':
        model = AgglomerativeClustering(
            n_clusters=2,
            linkage='ward',
            metric='euclidean',
        )
        return model.fit_predict(X)

    model = AgglomerativeClustering(
        n_clusters=2,
        linkage=linkage,
        metric=metric,
    )
    return model.fit_predict(X)


def best_map_labels(y_true, y_pred):

    p0, p1 = np.unique(y_pred)
    t0, t1 = np.unique(y_true)

    # var 1: p0->t0, p1->t1
    score_01 = ((y_pred == p0) & (y_true == t0)).sum() + ((y_pred == p1) & (y_true == t1)).sum()
    # var 2: p0->t1, p1->t0
    score_10 = ((y_pred == p0) & (y_true == t1)).sum() + ((y_pred == p1) & (y_true == t0)).sum()

    if score_01 >= score_10:
        return np.where(y_pred == p0, t0, t1)
    else:
        return np.where(y_pred == p0, t1, t0)


def compute_metrics(y_true, y_pred_mapped):

    scores = {
        'balanced_accuracy': balanced_accuracy_score(y_true, y_pred_mapped),
        'f1_macro': f1_score(y_true, y_pred_mapped, average='macro'),
        'mcc_scaled': (matthews_corrcoef(y_true, y_pred_mapped) + 1) / 2,
    }

    scores['weighted_score'] = (
        scores['balanced_accuracy'] * 0.35
        + scores['f1_macro'] * 0.35
        + scores['mcc_scaled'] * 0.3
    )

    return scores


def plot_metric_cm_heatmaps(df_res, cm, labels_cm, title_cm):

    metrics = [
        ('balanced_accuracy', 'Balanced accuracy', 'viridis'),
        ('f1_macro', 'F1 macro', 'magma'),
        ('mcc_scaled', 'MCC scaled [0, 1]', 'mako'),
    ]
    metric_cols = [m[0] for m in metrics]

    piv = df_res.pivot_table(index='linkage', columns='metric', values=metric_cols)

    fig, axes = plt.subplots(2, 2, figsize=(6, 5), dpi=160)
    ax1, ax2, ax3, ax_cm = axes.ravel()

    common_kws = dict(
        annot=True,
        fmt='.3f',
        vmin=0,
        vmax=1,
        cbar=False,
        annot_kws={'size': 6},
    )

    for ax, (col, title, cmap) in zip([ax1, ax2, ax3], metrics):
        sns.heatmap(piv[col], cmap=cmap, ax=ax, **common_kws)
        ax.set_title(title, fontsize=8)
        ax.set_xlabel('metric', fontsize=6)
        ax.set_ylabel('linkage', fontsize=6)
        ax.tick_params(axis='both', labelsize=6)

    sns.heatmap(
        cm,
        annot=True,
        annot_kws={'size': 6},
        fmt='d',
        cmap='Blues',
        cbar=False,
        ax=ax_cm,
    )
    ax_cm.set_title(title_cm, fontsize=8)
    ax_cm.set_xlabel('Predicted', fontsize=6)
    ax_cm.set_ylabel('True', fontsize=6)
    ax_cm.set_xticklabels(labels_cm, rotation=0)
    ax_cm.set_yticklabels(labels_cm, rotation=0)
    ax_cm.tick_params(axis='both', labelsize=6)

    fig.suptitle('Performance across metric × linkage and confusion matrix', y=1.02, fontsize=10)
    fig.tight_layout()
    plt.show()


def plot_clustermap(dfX, y_true, y_cluster_mapped, metric, linkage_method):

    true_levels = sorted(y_true.unique())
    cluster_levels = sorted(y_cluster_mapped.unique())

    pal_true = dict(zip(true_levels, sns.color_palette('Set2', len(true_levels))))
    pal_cluster = dict(zip(cluster_levels, sns.color_palette('Set3', len(cluster_levels))))

    row_colors = pd.DataFrame(
        {'true': y_true.map(pal_true), 'cluster_mapped': y_cluster_mapped.map(pal_cluster)},
        index=dfX.index,
    )

    g = sns.clustermap(
        dfX,
        row_colors=row_colors,
        cmap='vlag',
        center=0,
        figsize=(4, 5),
        yticklabels=False,
        xticklabels=True,
        col_cluster=True,
        row_cluster=True,
        metric=metric,
        method=linkage_method,
        cbar_pos=None
    )

    g.fig.suptitle(
        f'Clustermap: genefu / cluster_mapped; metric={metric}, linkage={linkage_method}',
        fontsize=10
    )

    g.fig.set_dpi(160)
    g.ax_col_dendrogram.set_visible(False)
    g.ax_heatmap.set(xlabel='', ylabel='')
    g.ax_heatmap.tick_params(axis='x', labelsize=6, rotation=90)
    g.ax_row_colors.set_axis_off()

    handles_true = [mpatches.Patch(color=pal_true[k], label=str(k)) for k in true_levels]
    handles_cluster = [mpatches.Patch(color=pal_cluster[k], label=str(k)) for k in cluster_levels]

    legend_kws = dict(
        loc='upper left',
        frameon=True,
        fontsize=6,
        title_fontsize=8,
        borderpad=0.8,
        labelspacing=0.6,
        handlelength=1.2,
        handleheight=1.2,
    )

    leg1 = g.ax_heatmap.legend(
        handles=handles_true, title='Genefu', bbox_to_anchor=(1.03, 0.98), **legend_kws
    )
    g.ax_heatmap.legend(
        handles=handles_cluster, title='Clustering', bbox_to_anchor=(1.03, 0.65), **legend_kws
    )
    g.ax_heatmap.add_artist(leg1)

    plt.show()


def cluster_breast_cancer_target_vs_rest(df_bg_genes_scaled, genes, y_subtype_genefu, target_label, metrics, linkages, plot):

    X_df = df_bg_genes_scaled.loc[:, genes].copy()
    X = X_df.to_numpy()

    not_target_label = f'Not_{target_label}'
    y_true = (y_subtype_genefu == target_label).map({True: target_label, False: not_target_label})

    rows = []
    clustering_results = {}

    for linkage in linkages:
        for metric in metrics:
            if linkage == 'ward' and metric != 'euclidean':
                continue

            y_pred_raw = predict_clusters(X, linkage, metric)
            y_pred_mapped = best_map_labels(y_true, y_pred_raw)
            scores = compute_metrics(y_true, y_pred_mapped)

            rows.append({'linkage': linkage, 'metric': metric, **scores, })
            clustering_results[(linkage, metric)] = {'y_pred_mapped': y_pred_mapped}

    results_df = pd.DataFrame(rows)
    results_df = results_df.sort_values(by='weighted_score', ascending=False)

    best_row = results_df.iloc[0]
    best_linkage = best_row['linkage']
    best_metric = best_row['metric']
    best_result = clustering_results[(best_linkage, best_metric)]
    y_pred_mapped_best = pd.Series(best_result['y_pred_mapped'], index=X_df.index)

    if plot:
        labels_cm = [target_label, not_target_label]
        cm = confusion_matrix(
            y_true.astype(str),
            y_pred_mapped_best,
            labels=labels_cm,
        )
        title_cm = (
            f'Winner: weighted_score={best_row["weighted_score"]:.3f}\n'
            f'linkage={best_linkage}, metric={best_metric}'
        )

        plot_metric_cm_heatmaps(
            df_res=results_df,
            cm=cm,
            labels_cm=labels_cm,
            title_cm=title_cm,
        )

        plot_clustermap(
            dfX=X_df,
            y_true=y_true,
            y_cluster_mapped=y_pred_mapped_best,
            metric=best_metric,
            linkage_method=best_linkage,
        )

    df_cluster_neg = df_bg_genes_scaled.loc[y_pred_mapped_best.eq(not_target_label)].copy()
    genefu_labels_for_non_target = y_subtype_genefu.loc[df_cluster_neg.index].copy()

    print('\nTop-3 combinations by general score:')
    print(results_df.head(3).to_string(index=False))

    return ClusteringResult(
        y_pred_mapped_best=y_pred_mapped_best,
        df_cluster_neg=df_cluster_neg,
        genefu_labels_for_non_target=genefu_labels_for_non_target,
    )
