# breast-cancer-tme-ccc





## Results

### Part 1. Breast cancer subtype classification from tumor-cell pseudobulks

**What for:**  
This step assigns PAM50 subtypes to donors and compares a standard PAM50 classifier with an RNA-seq-oriented alternative better suited to pseudobulk expression profiles.

**Pipeline:**  
Malignant cells were extracted from the atlas and aggregated into donor-level tumor pseudobulks. After TMM + logCPM normalization, subtypes were assigned in two ways: (1) PAM50 classification with `genefu`, and (2) a custom clustering-based workflow using published subtype-discriminative gene panels to separate Basal, ERBB2, LumA, and LumB tumors.

**Results:**  
The custom workflow produced a four-class solution across donors: **Basal-like (32), ERBB2+ (22), Luminal A (35), Luminal B (36)**.  
Agreement with `genefu` was strongest for **Basal-like**, moderate for **ERBB2+**, and less stable for the **Luminal A / Luminal B** split, which is expected given their biological similarity. At the same time, the clustering-based approach gave a clearer separation of luminal samples in the expression heatmap.

![Clustering vs genefu confusion matrix](imgs/clustering.png)

**Takeaway:**  
Subtype labels were robust enough for downstream analysis, and the custom RNA-seq-oriented approach provided a biologically interpretable alternative to microarray-based PAM50 classification.

**Notebook:**
More details can be found here: [Subtype classification notebook](notebooks/01_breast_cancer_subtype_classification.ipynb)

### Part 2. Tumor–microenvironment communication across PAM50 subtypes

**What for:**  
This part tests whether breast cancer subtypes differ not only in tumor-cell state, but also in how tumor and microenvironment compartments communicate.

**Pipeline:**  
After filtering donors and merging cell annotations into major cell types, subtype labels were mapped back to single cells and the dataset was prepared for communication analysis. Donor-level ligand-receptor interactions were inferred with **LIANA**, restricted to **tumor -> TME** and **TME -> tumor** signaling, and then decomposed with **Tensor-cell2cell** to identify shared communication programs across donors. In parallel, cell-type-specific donor pseudobulks were generated for later expression-based validation.

**Results:**  
After filtering, **121 of 138 donors** were retained. Subtype-dependent differences in TME composition were detected in **9 of 15** non-malignant major cell types, indicating that PAM50 subtypes differ not only transcriptionally, but also in microenvironment structure.

Tensor-cell2cell identified **one dominant communication program** that captured most subtype-related variation. This factor was strongest in **Basal-like**, elevated in **ERBB2+**, and weaker in **Luminal A/B** tumors. The same trend remained significant after restricting to curated interactions.

![Tensor score across subtypes](imgs/tensor_score_by_subtypes.png)

Biologically, the program combined:
- a **stromal-to-tumor axis** involving Fibroblast/Stromal, Endothelial, and Mural cells;
- a **tumor-to-immune axis** involving signaling toward Macrophages, Dendritic cells, Treg, and NK cells.

Top interactions included **COL1A1/COL1A2 → CD44**, **FN1 → CD44**, and **MIF → CD74/CXCR4**.

![Communication scheme](imgs/communication_scheme.png)

**Takeaway:**  
Breast cancer subtypes differ along a major tumor-TME communication gradient, with stronger stromal remodeling and immune-modulatory signaling in the more aggressive **Basal-like** and **ERBB2+** tumors.

**Notebook:**
This part was divided into two steps: data preparation and cell-cell communication analysis. Details can be found here:
- [CCC data preparation notebook](notebooks/02a_prep_data_pam50_for_liana.ipynb)
- [Cell-cell communication notebook](notebooks/02b_cell_cell_communication_analysis.ipynb)

### Part 3. Expression-level validation of subtype-associated communication patterns

**What for:**  
This step evaluates whether the communication patterns identified in Part 2 are also reflected in donor-level expression changes within specific cell compartments.

**Pipeline.**  
Cell-type-specific donor pseudobulks were analyzed in **R** using **edgeR** for differential expression across PAM50 subtypes. In addition, receptor-centered **Reactome** pathway scores were computed to test whether signaling programs highlighted by the CCC analysis were also altered at the transcriptional level.

**Results:**  
Subtype-associated expression differences were detected not only in malignant cells, but also in several TME compartments. In the cell types most strongly involved in the communication signal, the same broad pattern was observed as in Part 2: expression programs were generally higher in **Basal-like** and **ERBB2+** tumors and lower in **Luminal** subtypes.

Pathway analysis provided targeted support for this result. Significant subtype-associated differences were detected in:
- **Malignant cells**: stromal and inflammatory programs;
- **Macrophages**: inflammatory program.

These findings are consistent with the dominant CCC program, which linked a **stromal TME -> tumor axis** with a **tumor -> immune axis**.

![DGE heatmaps](imgs/dge.png)

**Takeaway:**  
The communication gradient identified by LIANA + Tensor-cell2cell is supported by independent expression- and pathway-level evidence, suggesting that PAM50 subtypes differ in both tumor state and tumor-microenvironment interaction programs.

## Reproducible CCC pipeline

This repository includes a standalone CLI pipeline, `ccc_pipeline.py`, for reproducible **LIANA + Tensor-cell2cell** analysis from an annotated `.h5ad` object.

The pipeline runs donor-level **LIANA** inference, applies subtype-aware recurrence filtering, restricts interactions by communication direction, builds a 4D communication tensor, and performs **Tensor-cell2cell** decomposition. It is designed for flexible re-analysis of tumor microenvironment communication under different filtering settings.

**Input requirements:**
- donor/sample annotation
- cell type annotation
- subtype annotation

**Supported interaction modes:**
- `tt` - tumor -> tumor
- `t_nt` - tumor -> non-tumor
- `nt_t` - non-tumor -> tumor
- `nt_nt` - non-tumor -> non-tumor

**Example**
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

All main intermediate and final outputs are saved to the specified output directory.  

For the full list of options:
```bash
python ccc_pipeline.py --help
```

