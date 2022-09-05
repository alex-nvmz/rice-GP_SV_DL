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