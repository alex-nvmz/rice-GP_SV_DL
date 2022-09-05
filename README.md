# Genomic Prediction of Complex Traits in Rice: Leveraging Structural Variation and Applying Deep Learning Methods

**Abstract**  
Prediction of complex traits is challenging for plant breeders due to the high number of quantitative trait loci (QTLs) affecting them and because of the relevance of non-additive genetic effects for some traits. Therefore, it is important to fine-tune this process as much as possible. We present a case study of four traits in *Oryza sativa*, carried out with a dataset of 738 accessions representing the whole genetic diversity of the species. We address the problem of the high number of QTLs by introducing structural variants (SVs) as an additional source of genetic variation with respect to single-nucleotide polymorphisms (SNPs). The SVs dataset includes transposon insertion polymorphisms (TIPs), deletions, tandem duplications and inversions. The problem of the non-additive genetic effects is tackled by implementing deep learning methods, which are based on artificial neural network (ANN) models. These models can implement complex non-linear transformations of an input, so they should be able to capture the effects of dominance and epistasis. Results show that SVs can substantially improve the prediction of certain traits, as seen by the 7.40% improvement for leaf senescence. Deep learning methods proved to be competitive with linear models for regression tasks. They also performed well for classification, in which the linear models performed sub-optimally.

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

| Script                         | Description                                                               |
| ------------------------------ | ------------------------------------------------------------------------- |
| 01_accessions_phenotypes.R     | Process original files with phenotypic and rice population information    |
| 02_SV-info.R                   | Summarize information from the original SV files                          |
| 03_geno-matrix_SNP.R           | Turn SNPs .bed file into .csv genotype matrix                             |
| 04_geno-matrix_TIP.R           | Generate TIP genotype matrices and filter by MAF                          |
| 05_geno-matrix_SV.R            | Generate SV genotype matrices from the information files                  |
| 06_CNV-processing.sh           | Compute missing rate of the CNVs (data eventually discarded)              |
| 07_matrices_combined-markers.R | Concatenation of genotype matrices. Creates MITE-DTX and RLX-RIX matrices |
| 08_partitions.R                | Establishes the data partitions for the 10-fold CV                        |
| 09a_GWAS_param-grid.R          | Creates grid where every row has the parameters for a GWAS run            |
| 09b_GWAS.R                     | GWAS core script (1 trait and 1 partition at a time)                      |
| 09c_run_GWAS.sh                | Runs all GWAS jobs with the option to parallelize                         |
| 10_top-markers_tables.R        | Process GWAS output to generate tables with the top 10k markers           |
| 11_top-markers_geno-matrices.R | Generate sorted genotype matrices with the top 10k markers                |
| 12_kernels-from-markers.R      | Compute genomic relationship matrices with AGHmatrix                      |
| 13_PCA.R                       | Perform PCA of the genotype matrices as a exploratory analysis            |
| 14_kernel-PCA.R                | Obtain kerPC inputs for the ANNs                                          |
| s01_plot_TIP-freq.R            | Plot TIP allelic frequencies before filtering by MAF                      |
| s02_plot_CNV-stats.R           | Plot CNV missing rate                                                     |
| s03_top-markers_plots.R        | Plot marker type composition of the top 10k marker genotype matrices      |
| s04_trait-plots.R              | Exploratory analysis of the trait variables                               |

### genetic-variance_gaussian

| Script                    | Description                                           |
| ------------------------- | ----------------------------------------------------- |
| 01_h2-estimate_gaussian.R | Run Bayesian RKHS models to estimate genomic variance |
| 02_merge-results_plot.R   | Merge the results into a table and plot               |

### GP-linear

| Script                    | Description                                                                            |
| ------------------------- | -------------------------------------------------------------------------------------- |
| 01a_do-param-grid.R       | Create grids where every row has the parameters for an independent RKHS or Bayes C run |
| 01b_GP-BGLR.R             | Core BGLR script for RKHS or Bayes C                                                   |
| 01c_run-GP.sh             | Runs all GP jobs with the option to parallelize                                        |
| 01c_bash_instruction.txt  | Bash commands used with 01c_run-GP.sh                                                  |
| 02_merge-results_plot.R   | Merge results into a table and plot                                                    |
| 03_compute-accuracy-AUC.R | Compute accuracy and AUC metric a posteriori using the .RData objects of each analysis |

## Python

| Script                        | Description                                                                                                                                                     |
| ----------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| environment.yml               | Specifies the conda environment used. To replicate, use command: `conda env create -f environment.yml`                                                          |
| 01_parameters.ipynb           | Creates grids where every row has the parameters for an independent ANN run                                                                                     |
| hypermodels.py                | Defines subclasses of the keras_tuner HyperModel(): hyperMLP() and hyperCNN(). They specify the hyperparameter space and have methods to build and train models |
| 02_run_hypermodel_final.ipynb | Run hyperparameter search and early prediction results. Interactive version                                                                                     |
| 02_run_hypermodel_final.py    | Run hyperparameter search and early prediction results                                                                                                          |
| 03_run_prediction_final.ipynb | Run prediction with a given set of hyperparameters. Interactive version                                                                                         |
| 03_run_prediction_final.py    | Run prediction with a given set of hyperparameters                                                                                                              |
| 04_run_python_script.sh       | Run either 02_run_hypermodel_final.py or 03_run_prediction_final.py one job at a time                                                                           |
| 04_bash_instruction.txt       | Bash commands used with 04_run_python_script.sh                                                                                                                 |
| 05_get_results.ipynb          | Gather all results into a table                                                                                                                                 |
| 06_plot_results.R             | Make final plots with linear models + ANNs results                                                                                                              |
