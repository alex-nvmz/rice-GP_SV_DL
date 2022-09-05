### GP-linear

| Script                    | Description                                                                            |
| ------------------------- | -------------------------------------------------------------------------------------- |
| 01a_do-param-grid.R       | Create grids where every row has the parameters for an independent RKHS or Bayes C run |
| 01b_GP-BGLR.R             | Core BGLR script for RKHS or Bayes C                                                   |
| 01c_run-GP.sh             | Runs all GP jobs with the option to parallelize                                        |
| 01c_bash_instruction.txt  | Bash commands used with 01c_run-GP.sh                                                  |
| 02_merge-results_plot.R   | Merge results into a table and plot                                                    |
| 03_compute-accuracy-AUC.R | Compute accuracy and AUC metric a posteriori using the .RData objects of each analysis |