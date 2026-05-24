# breast-cancer-tme-ccc





## Results

### Part 1: Breast cancer subtype classification from tumor-cell pseudobulks

The goal is to compare two subtype assignment strategies:

- **PAM50 classification with `genefu`**, the standard approach;
- a **custom RNA-seq-oriented classifier** based on sequential hierarchical clustering with published subtype-discriminative gene panels.

The main limitation of `genefu` in this setting is that its PAM50 reference centroids were derived from **microarray** data, whereas this project works with **RNA-seq* data.

#### Analysis notebook

The first part of the analysis is available here: [Subtype classification notebook](notebooks/01_breast_cancer_subtype_classification.ipynb)

#### What is done in this notebook

- malignant cells are extracted from the full dataset;
- donor-level tumor pseudobulks are constructed;
- pseudobulks are normalized with **TMM + logCPM**;
- donor samples are classified with **`genefu` PAM50**;
- an alternative classifier is built using published genes from a sequence of binary boosted trees;
- these gene panels are used in a **series of agglomerative clustering steps**:
  - Basal vs rest
  - ERBB2 vs rest
  - LumB vs LumA

#### Main result

The custom clustering-based workflow produces a four-class solution:

- **Basal**: 32
- **ERBB2**: 22
- **LumA**: 35
- **LumB**: 36

Overall, the clustering-based labels show strong agreement with `genefu` for the **Basal** group, moderate agreement for **ERBB2**, and a less stable split between **LumA** and **LumB**. This is expected, since LumA and LumB are biologically closer. At the same time, the heatmap suggests that the clustering-based approach separates luminal samples more clearly than `genefu`.

#### Comparison with `genefu`

![Clustering vs genefu confusion matrix](pics/clustering.png)

### Part 2a: Data preparation for cell-cell communication analysis

The goal of this step is to prepare the breast cancer scRNA-seq atlas for downstream **cell-cell communication analysis**.

This notebook prepares the breast cancer scRNA-seq atlas for downstream **LIANA + Tensor-cell2cell** analysis and generates **cell-type-specific pseudobulks** for subsequent **edgeR-based validation**.

#### Analysis notebook

The preparation step is available here: [CCC data preparation notebook](notebooks/02a_prep_data_pam50_for_liana.ipynb)

#### What is done in this notebook

- the breast cancer atlas is filtered to retain donors with sufficient **tumor** and **non-tumor** cell counts;
- **PAM50 subtype labels** from Part 1 are mapped to donor cells;
- subtype structure is inspected on UMAPs to check for **batch effects** and **donor-specific tumor states**;
- fine cell annotations are collapsed into **16 major cell types** for communication analysis;
- subtype differences in **major cell type composition** are tested;
- the expression matrix is reformatted to keep only **HGNC-mapped genes** and collapse duplicated symbols;
- filtered AnnData objects are prepared for **LIANA**;
- **cell-type-specific donor pseudobulks** are generated for downstream **edgeR differential expression** and validation of the **LIANA + Tensor-cell2cell** results.

#### Main result

1. After filtering, **121 of 138 donors** were retained for downstream analysis. The subtype labels are well distributed across donors and batches.

2. The original atlas contained **29 cell types**, many of them represented by very small numbers of cells per donor.  
For the communication analysis, these were collapsed into **16 major groups**.

3. TME composition differences across PAM50 subtypes: significant subtype-dependent differences were detected for **9 of 15 non-malignant major cell types**. This indicates that PAM50 subtypes differ not only in tumor-cell expression states, but also in their TME composition.

4. Preparation of HGNC-compatible gene annotations
Because LIANA expects **HGNC symbols**, the expression matrix was reformatted:
    - genes lacking HGNC symbols were removed;
    - duplicated HGNC entries were collapsed by retaining the gene with the highest total expression.

    The resulting filtered objects were saved to disk and will be used as input for the downstream LIANA + Tensor-cell2cell analysis in Part 2b.

5. In addition, **donor-level pseudobulks stratified by cell type** were generated.

    These pseudobulks are intended for:
    - downstream differential expression analysis with **edgeR**;
    - validation of subtype-associated communication patterns detected by **LIANA + Tensor-cell2cell**;
    - follow-up testing of whether ligand, receptor, or pathway-related signals highlighted in the communication analysis are also reflected at the expression level within specific cell compartments.

### Part 2b: Cell-cell communication analysis across PAM50 breast cancer subtypes

The goal of this part is to identify **cell-cell communication programs** associated with breast cancer PAM50 subtypes and to determine which cell types and ligand-receptor interactions contribute most strongly to these subtype-specific patterns.

This analysis combines:
- **LIANA**, to infer donor-level ligand-receptor interactions;
- **Tensor-cell2cell**, to decompose these interactions into latent communication programs shared across donors.

The workflow focuses specifically on **tumor-microenvironment crosstalk**, restricting the analysis to:
- **tumor -> TME**
- **TME -> tumor**

#### Analysis notebook

The second part of the analysis is available here: [Cell-cell communication notebook](notebooks/02b_cell_cell_communication_analysis.ipynb)

#### What is done in this notebook

- donor-level cell-cell communication is inferred with **LIANA** using the `consensus` resource;
- interactions are retained only if they are present within each PAM50 subtype in at least 40% of donors;
- the analysis is restricted to tumor-TME interactions;
- LIANA results are converted into a 4D tensor:
  - donors (context)
  - ligand-receptor pairs
  - sender cell types
  - receiver cell types
- the tensor is decomposed with **Tensor-cell2cell** to identify latent communication factors;
- donor factor loadings are compared across PAM50 subtypes;
- factor-derived scores are assigned back to individual LIANA interactions;
- interactions are additionally filtered using **CellChatDB** to keep curated and pathway-annotated pairs;
- subtype differences are analyzed separately for:
  - **TME -> Tumor**
  - **Tumor -> TME**

#### Main result

Tensor-cell2cell identified **one dominant communication factor**, suggesting that the main variation in tumor-TME signaling can be summarized by a single latent program.

This program differs significantly across PAM50 subtypes:
- it is **strongest in Basal-like**
- also enriched in **ERBB2+**
- weaker in **Luminal A** and **Luminal B**

![Tensor score across subtypes](pics/tensor_score_by_subtypes.png)

After filtering to curated interactions, the same subtype trend remains significant both for both TME -> Tumor and Tumor -> TME signaling.

#### Biological interpretation
The inferred factor combines two major components:
- a **stromal-to-tumor axis**, dominated by Fibroblast/Stromal, Endothelial, and Mural cells;
- a **tumor-to-immune axis**, dominated by malignant-cell signaling toward Macrophage, Dendritic_cell, Treg, NK.

Top LR-interactions:
| TME -> Tumor (Stromal) | Tumor -> TME (Immune) |
|---|---|
| COL1A1 \ COL1A2 ^ CD44** | MIF ^ CD74_CXCR4 |
| FN1 ^ CD44 | MIF ^ CD44_CD74 |
| APP ^ CD74 | APP ^ CD74 |

Overall, the results point to a communication program combining **stromal remodeling** and **immune-modulatory signaling**, with higher activity in the more aggressive **Basal-like** and **ERBB2+** tumors and lower activity in the **Luminal** subtypes.

### CCC pipeline

This repository includes a standalone **LIANA + Tensor-cell2cell** CLI pipeline (`ccc_pipeline.py`) for reproducible cell-cell communication analysis from an annotated `.h5ad` object.

The pipeline is designed for studying communication programs in the **tumor microenvironment (TME)**, with explicit support for separating **tumor** and **non-tumor** cell compartments.

#### Why use it

The pipeline was designed to make CCC analysis easy to rerun with different:

- LIANA filtering parameters,
- subtype recurrence thresholds,
- tumor/non-tumor interaction constraints,
- CPU/GPU settings.

#### Input

The script expects an annotated single-cell `h5ad` object with:

- a donor/sample column,
- a cell type annotation column,
- a subtype annotation column.

#### Workflow

Starting from the input `h5ad`, the pipeline:

1. loads the dataset;
2. runs **LIANA** separately for each donor/sample using the **consensus** resource;
3. adds subtype metadata to donor-level LIANA results;
4. applies subtype-level recurrence filtering across donors;
5. filters interactions by communication direction;
6. builds a 4D communication tensor;
7. runs **Tensor-cell2cell** decomposition with elbow-based rank selection;
8. saves all main outputs to disk.

#### Interaction modes

The `--interaction-mode` argument controls which communication directions are retained:

- `tt` - tumor -> tumor
- `t_nt` - tumor -> non-tumor
- `nt_t` - non-tumor -> tumor
- `nt_nt` - non-tumor -> non-tumor

Multiple modes can be combined, for example:

```bash
--interaction-mode 't_nt,nt_t'
```

#### Example run

```bash
python ccc_pipeline.py \
    --adata-path ../data/data_hgnc.h5ad \
    --outdir ../data/liana_tensor_cell2cell_result \
    --min-cells 20 \
    --expr-prop 0.2 \
    --min-donors-prop 0.4 \
    --tumor-label 'Malignant cell' \
    --subtype-col PAM50_subtype \
    --sample-col donor_id \
    --celltype-col cell_type_major \
    --interaction-mode 't_nt,nt_t' \
    --use-gpu
```

#### Main CLI arguments

**Input / output**
- `--adata-path` - path to the input `.h5ad` file
- `--outdir` - directory for all output files

**LIANA parameters**
- `--min-cells` - minimum number of cells per cell type
- `--expr-prop` - minimum expression proportion threshold

**Metadata columns**
- `--subtype-col` - subtype label column
- `--sample-col` - donor/sample identifier column
- `--celltype-col` - cell type annotation column
- `--tumor-label` - label used to define tumor cells

**Filtering**
- `--min-donors-prop` - minimum donor proportion within each subtype required to retain an interaction
- `--interaction-mode` - interaction classes to retain

**Hardware**
- `--use-gpu` / `--no-use-gpu` - whether to run tensor decomposition on GPU if available

For the full list of options:

```bash
python ccc_pipeline.py --help
```

#### Output

The pipeline writes multiple intermediate and final files to the output directory, including:

- full LIANA results;
- LIANA results after subtype recurrence filtering;
- LIANA results after interaction-direction filtering;
- serialized **tensor** object;
- serialized **tensor metadata** object.

Output file names automatically include the parameter settings used in the run, making it easier to track and compare analyses.
