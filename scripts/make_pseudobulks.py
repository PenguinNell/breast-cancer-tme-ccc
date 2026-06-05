"""
Utilities for pseudobulk construction from single-cell data.

This module aggregates raw counts by donor and cell type group.
Can optionally save the resulting count matrix and sample metadata.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData


def make_pseudobulks(
        adata: AnnData,
        donor_col: str,
        group_col: str,
        min_cells: int = 20,
        return_pb: bool = False,
        save: bool = True,
        out_dir: str | Path | None = None,
        counts_name: str = 'pseudobulks_counts_gene_by_sample.tsv',
        meta_name: str = 'pseudobulks_sample_metadata.tsv',
) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """Aggregate raw counts by donor and cell type group."""

    obs = adata.obs.loc[:, [donor_col, group_col]].copy()
    obs[donor_col] = obs[donor_col].astype(str)
    obs[group_col] = obs[group_col].astype(str)
    obs['sample_id'] = obs[donor_col] + '__' + obs[group_col]

    sample_ids = obs['sample_id'].to_numpy()
    uniq, inv = np.unique(sample_ids, return_inverse=True)

    n_samples = len(uniq)
    n_genes = adata.raw.shape[1]

    pb = np.zeros((n_samples, n_genes), dtype=np.int64)

    for i in range(n_samples):
        rows = np.where(inv == i)[0]
        pb[i, :] = np.asarray(adata.raw.X[rows].sum(axis=0)).ravel().astype(np.int64)

    meta_pb = (
        obs.groupby('sample_id', observed=True)
        .agg({
            donor_col: 'first',
            group_col: 'first',
        })
        .reindex(uniq)
        .reset_index()
    )

    meta_pb['n_cells'] = (
        obs.groupby('sample_id', observed=True)
        .size()
        .reindex(uniq)
        .to_numpy()
    )

    keep = meta_pb['n_cells'].ge(min_cells).to_numpy()
    pb = pb[keep, :]
    pb_df = pd.DataFrame(
        pb,
        index=meta_pb['sample_id'],
        columns=adata.raw.var_names
    )
    meta_pb = meta_pb.loc[keep].copy().reset_index(drop=True)

    if save:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        pb_df.to_csv(out_dir / counts_name, sep='\t')
        meta_pb.to_csv(out_dir / meta_name, sep='\t', index=False)

    if return_pb:
        return pb_df, meta_pb
