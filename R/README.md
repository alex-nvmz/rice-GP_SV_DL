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
