# Genomic Prediction of Complex Traits in Rice: Leveraging Structural Variation and Applying Deep Learning Methods

**Abstract**  
Prediction of complex traits is challenging for plant breeders due to the high number of quantitative trait loci (QTLs) affecting them and because of the relevance of non-additive genetic effects for some traits. Therefore, it is important to fine-tune this process as much as possible. We present a case study of four traits in *Oryza sativa*, carried out with a dataset of 738 accessions representing the whole genetic diversity of the species. We address the problem of the high number of QTLs by introducing structural variants (SVs) as an additional source of genetic variation with respect to single-nucleotide polymorphisms (SNPs). The SVs dataset includes transposon insertion polymorphisms (TIPs), deletions, tandem duplications and inversions. The problem of the non-additive genetic effects is tackled by implementing deep learning methods, which are based on artificial neural network (ANN) models. These models can implement complex non-linear transformations of an input, so they should be able to capture the effects of dominance and epistasis. Results show that SVs can substantially improve the prediction of certain traits, as seen by the 7.40% improvement for leaf senescence. Deep learning methods proved to be competitive with linear models for regression tasks. They also performed well for classification, in which the linear models struggled.

This repository contains the scripts used to develop the project. Source data and results are not included because of the large size of the files.

| Folder | Content                                                                                       |
| ------ | --------------------------------------------------------------------------------------------- |
| R      | Data processing, exploratory data analysis, genomic variance inference, GP with linear models |
| Python | Genomic prediction with ANNs                                                                  |

## R

| Folder                    | Content                                    |
| ------------------------- | ------------------------------------------ |
| data                      | Data processing, exploratory data analysis |
| genetic-variance_gaussian | genomic variance inference                 |
| GP-linear                 | GP with linear models                      |

### data

| Script                         | Description |
| ------------------------------ | ----------- |
| 01_accessions_phenotypes.R     |             |
| 02_SV-info.R                   |             |
| 03_geno-matrix_SNP.R           |             |
| 04_geno-matrix_TIP.R           |             |
| 05_geno-matrix_SV.R            |             |
| 06_CNV-processing.sh           |             |
| 07_matrices_combined-markers.R |             |
| 08_partitions.R                |             |
| 09a_GWAS_param-grid.R          |             |
| 09b_GWAS.R                     |             |
| 09c_run_GWAS.sh                |             |
| 10_top-markers_tables.R        |             |
| 11_top-markers_geno-matrices.R |             |
| 12_kernels-from-markers.R      |             |
| 13_PCA.R                       |             |
| 14_kernel-PCA.R                |             |
| 15_trait-not-NA-count.R        |             |
| s01_plot_TIP-freq.R            |             |
| s02_plot_CNV-stats.R           |             |
| s03a_top-markers_plots.R       |             |
| s03b_top-markers_plots.R       |             |
| s04_trait-plots.R              |             |

### genetic-variance_gaussian

| Script                    | Description |
| ------------------------- | ----------- |
| 01_h2-estimate_gaussian.R |             |
| 02_merge-results_plot.R   |             |

### GP-linear

| Script                    | Description |
| ------------------------- | ----------- |
| 01a_do-param-grid.R       |             |
| 01b_GP-BGLR.R             |             |
| 01c_run-GP.sh             |             |
| 02a_plot-BayesC-Htune.R   |             |
| 02b_plot.R                |             |
| 02c_plot_alt.R            |             |
| 03_compute-accuracy-AUC.R |             |
| bash_instruction.txt      |             |


## Python

| Script                        | Description |
| ----------------------------- | ----------- |
| 01_parameters.ipynb           |             |
| hypermodels.py                |             |
| 02_run_hypermodel_final.ipynb |             |
| 02_run_hypermodel_final.py    |             |
| 03_run_prediction_final.ipynb |             |
| 03_run_prediction_final.py    |             |
| 04_run_python_script.sh       |             |
| 04_bash_instruction.txt       |             |
| 05_get_results.ipynb          |             |
| 06_plot_results.R             |             |
