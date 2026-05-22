# breast-cancer-tme-ccc





## Results

### Part 1: breast cancer subtype classification from tumor-cell pseudobulks

The goal is to compare two subtype assignment strategies:

- **PAM50 classification with `genefu`**, the standard approach;
- a **custom RNA-seq-oriented classifier** based on sequential hierarchical clustering with published subtype-discriminative gene panels.

The main limitation of `genefu` in this setting is that its PAM50 reference centroids were derived from **microarray** data, whereas this project works with **RNA-seq* data.

#### Analysis notebook

The first part of the analysis is available here:

- [Subtype classification notebook](notebooks/Nbreast_cancer_subtype_classification.ipynb)

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
