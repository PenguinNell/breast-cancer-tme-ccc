"""
LIANA and Tensor-cell2cell pipeline for cell-cell communication analysis.

This module runs donor-level LIANA, applies subtype-level and interaction filters,
and builds a communication tensor for Tensor-cell2cell.
"""

import argparse
import gc
import pickle
import warnings
from pathlib import Path

import pandas as pd

warnings.filterwarnings('ignore')


def get_device(use_gpu: bool) -> str:
    """Return the compute device."""

    if not use_gpu:
        return 'cpu'

    import torch
    import tensorly as tl

    if torch.cuda.is_available():
        tl.set_backend('pytorch')
        return 'cuda'
    return 'cpu'


def parse_interactions(mode: str) -> set[str]:
    """Parse interaction-mode string."""

    valid_base = {'tt', 't_nt', 'nt_t', 'nt_nt'}
    mode = mode.strip()

    if mode == 'all':
        return valid_base

    if mode.startswith('all_except_'):
        excluded = mode.replace('all_except_', '')
        if excluded not in valid_base:
            raise ValueError(f'Unknown interaction mode: {mode}')
        return valid_base - {excluded}

    parts = {part.strip() for part in mode.split(',') if part.strip()}
    if not parts <= valid_base:
        raise ValueError(f'Unknown interaction mode: {mode}')
    return parts


def make_param_suffix(
        min_cells: int,
        expr_prop: float,
        min_donors_prop: float,
        interaction_mode_tag: str,
) -> str:
    """Build the parameter suffix for output files."""

    return '__'.join([
        f'mincell_{min_cells}',
        f'exprprop_{str(expr_prop).replace(".", "p")}',
        f'mindonorsprop_{str(min_donors_prop).replace(".", "p")}',
        f'interactmode_{interaction_mode_tag}',
    ])


def filter_interactions(
        df: pd.DataFrame,
        interaction_mode: set[str],
        tumor_label: str,
) -> pd.DataFrame:
    """Filter interactions by sender and receiver type."""

    source_is_tumor = df['source'] == tumor_label
    target_is_tumor = df['target'] == tumor_label

    masks = {
        'tt': source_is_tumor & target_is_tumor,
        't_nt': source_is_tumor & ~target_is_tumor,
        'nt_t': ~source_is_tumor & target_is_tumor,
        'nt_nt': ~source_is_tumor & ~target_is_tumor,
    }

    final_mask = masks[next(iter(interaction_mode))].copy()
    for key in interaction_mode:
        final_mask |= masks[key]

    return df.loc[final_mask].copy()


def run_liana(
        adata_path: str | Path,
        outdir: Path,
        min_cells: int,
        expr_prop: float,
        min_donors_prop: float,
        interaction_mode: set[str],
        subtype_col: str,
        sample_col: str,
        celltype_col: str,
        param_suffix: str,
        tumor_label: str,
        interaction_mode_tag: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run LIANA and apply post-processing filters."""

    import scanpy as sc # too large to load at startup
    import liana as li # too large to load at startup

    adata = sc.read_h5ad(adata_path)
    context_df = adata.obs[[sample_col, subtype_col]].drop_duplicates()

    res = li.mt.rank_aggregate.by_sample(
        adata,
        groupby=celltype_col,
        resource_name='consensus',
        sample_key=sample_col,
        use_raw=False,
        expr_prop=expr_prop,
        min_cells=min_cells,
        verbose=True,
        inplace=False,
    )

    del adata
    gc.collect()

    res = res.merge(context_df, on=sample_col, how='left')
    res.to_csv(outdir / f'liana_res_full__{param_suffix}.csv', index=False)

    donors_per_subtype = context_df.groupby(subtype_col)[sample_col].nunique()

    interaction_coverage = (
        res[[sample_col, subtype_col, 'source', 'target', 'ligand_complex', 'receptor_complex']]
        .drop_duplicates()
        .groupby([subtype_col, 'source', 'target', 'ligand_complex', 'receptor_complex'])[sample_col]
        .nunique()
        .reset_index(name='n_donors_with_interaction')
    )

    interaction_coverage['n_donors_subtype'] = interaction_coverage[subtype_col].map(donors_per_subtype).astype(int)
    interaction_coverage['donor_prop'] = (
            interaction_coverage['n_donors_with_interaction'] / interaction_coverage['n_donors_subtype']
    )

    interactions_to_keep = interaction_coverage.loc[
        interaction_coverage['donor_prop'] >= min_donors_prop,
        [subtype_col, 'source', 'target', 'ligand_complex', 'receptor_complex']
    ]

    res_min_donors = res.merge(
        interactions_to_keep,
        on=[subtype_col, 'source', 'target', 'ligand_complex', 'receptor_complex'],
        how='inner'
    )

    res_min_donors.to_csv(outdir / f'liana_res_min_donors_filtered__{param_suffix}.csv', index=False)

    liana_res_selected = filter_interactions(
        df=res_min_donors,
        interaction_mode=interaction_mode,
        tumor_label=tumor_label,
    )

    selected_result_name = f'liana_res_filtered_by_{interaction_mode_tag}__{param_suffix}.csv'
    liana_res_selected.to_csv(outdir / selected_result_name, index=False)

    summary = {
        'Rows in full LIANA result': len(res),
        'Rows after min-donors filter': len(res_min_donors),
        'Rows after interaction-mode filter': len(liana_res_selected),
        'Unique LR pairs': len(liana_res_selected[['ligand_complex', 'receptor_complex']].drop_duplicates()),
        'Unique source-target-LR combinations': len(
            liana_res_selected[['source', 'target', 'ligand_complex', 'receptor_complex']].drop_duplicates()
        ),
    }
    print('\n'.join(f'{k}: {v}' for k, v in summary.items()))

    return liana_res_selected, context_df


def run_tensor(
        liana_res_filt: pd.DataFrame,
        context_df: pd.DataFrame,
        sample_col: str,
        subtype_col: str,
        outdir: Path,
        param_suffix: str,
        device: str,
) -> None:
    """Build and decompose the communication tensor."""

    import liana as li # too large to load at startup
    import cell2cell as c2c # too large to load at startup

    context_map = context_df.set_index(sample_col)[subtype_col].to_dict()

    print('\nBuilding a Tensor')
    tensor = li.multi.to_tensor_c2c(
        liana_res=liana_res_filt,
        sample_key=sample_col,
        score_key='magnitude_rank',
        non_expressed_fill=0,
        how='outer',
    )
    print(f'Tensor shape: {tensor.tensor.shape}')

    print('\nBuild Metadata')
    tensor_meta = c2c.tensor.generate_tensor_metadata(
        interaction_tensor=tensor,
        metadata_dicts=[context_map, None, None, None],
        fill_with_order_elements=True,
    )

    print('\nRunning Tensor-cell2cell')
    print(f'Device: {device}')
    tensor = c2c.analysis.run_tensor_cell2cell_pipeline(
        tensor,
        tensor_meta,
        copy_tensor=True,
        rank=None,
        tf_optimization='regular',
        random_state=0,
        device=device,
        elbow_metric='error',
        smooth_elbow=False,
        upper_rank=20,
        tf_init='random',
        tf_svd='numpy_svd',
        cmaps=None,
        sample_col='Element',
        group_col='Category',
        output_fig=False,
    )

    with open(outdir / f'tensor__{param_suffix}.pkl', 'wb') as f:
        pickle.dump(tensor, f)

    with open(outdir / f'tensor_meta__{param_suffix}.pkl', 'wb') as f:
        pickle.dump(tensor_meta, f)


def run_pipeline(
        adata_path: str | Path,
        outdir: str | Path = "liana_tensor_cell2cell_result",
        min_cells: int = 5,
        expr_prop: float = 0.1,
        min_donors_prop: float = 0.3,
        tumor_label: str = "malignant cell",
        subtype_col: str = "subtype",
        sample_col: str = "sample",
        celltype_col: str = "celltype",
        interaction_mode: str = "all",
        use_gpu: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the LIANA and Tensor-cell2cell pipeline."""

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    interaction_mode_raw = interaction_mode.strip()
    interaction_mode_set = parse_interactions(interaction_mode_raw)
    interaction_mode_tag = interaction_mode_raw.replace(',', '+')

    device = get_device(use_gpu)

    param_suffix = make_param_suffix(
        min_cells=min_cells,
        expr_prop=expr_prop,
        min_donors_prop=min_donors_prop,
        interaction_mode_tag=interaction_mode_tag,
    )

    print('\n=== Starting CCC pipeline ===')
    print(f'Params for LIANA: min_cells={min_cells}, expr_prop={expr_prop}')
    print(f'Filters: min_donors_prop={min_donors_prop}, interaction_mode="{interaction_mode}"')
    print('Running LIANA')

    liana_res, context_df = run_liana(
        adata_path=adata_path,
        outdir=outdir,
        min_cells=min_cells,
        expr_prop=expr_prop,
        min_donors_prop=min_donors_prop,
        interaction_mode=interaction_mode_set,
        subtype_col=subtype_col,
        sample_col=sample_col,
        celltype_col=celltype_col,
        param_suffix=param_suffix,
        tumor_label=tumor_label,
        interaction_mode_tag=interaction_mode_tag,
    )

    run_tensor(
        liana_res_filt=liana_res,
        context_df=context_df,
        sample_col=sample_col,
        subtype_col=subtype_col,
        outdir=outdir,
        param_suffix=param_suffix,
        device=device,
    )

    print('Done')

    return liana_res, context_df


def main() -> None:
    """Parse CLI arguments and run the pipeline."""

    parser = argparse.ArgumentParser()

    parser.add_argument('--adata-path', required=True)
    parser.add_argument('--outdir', default='liana_tensor_cell2cell_result')
    parser.add_argument('--min-cells', default=5, type=int)
    parser.add_argument('--expr-prop', default=0.1, type=float)
    parser.add_argument('--min-donors-prop', default=0.3, type=float)
    parser.add_argument('--tumor-label', required=True)
    parser.add_argument('--subtype-col', required=True)
    parser.add_argument('--sample-col', required=True)
    parser.add_argument('--celltype-col', required=True)
    parser.add_argument(
        '--interaction-mode',
        default='all',
        help='all, all_except_tt, all_except_t_nt, all_except_nt_t, all_except_nt_nt, or comma-separated subset of tt,t_nt,nt_t,nt_nt'
    )

    gpu_group = parser.add_mutually_exclusive_group()
    gpu_group.add_argument('--use-gpu', dest='use_gpu', action='store_true')
    gpu_group.add_argument('--no-use-gpu', dest='use_gpu', action='store_false')
    parser.set_defaults(use_gpu=False)

    args = parser.parse_args()

    run_pipeline(
        adata_path=args.adata_path,
        outdir=args.outdir,
        min_cells=args.min_cells,
        expr_prop=args.expr_prop,
        min_donors_prop=args.min_donors_prop,
        tumor_label=args.tumor_label,
        subtype_col=args.subtype_col,
        sample_col=args.sample_col,
        celltype_col=args.celltype_col,
        interaction_mode=args.interaction_mode,
        use_gpu=args.use_gpu,
    )


if __name__ == '__main__':
    main()
