# %%
import sys
import os
import shutil
import json
import ast

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
plt.style.use("seaborn-whitegrid")
import seaborn as sns
sns.set_style("whitegrid")

import tensorflow as tf
print("tensorflow:", tf.__version__)
from tensorflow import keras
print("keras:", keras.__version__)
import keras_tuner as kt
print("keras_tuner:", kt.__version__)

from keras import layers
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix

# Check if tensorflow can use your GPU
tf.config.list_physical_devices('GPU')

# %%
# Do prediction with the best hpars of current partition or AroAdm partition

# which_hps = part
which_hps = "AroAdm"
print(which_hps)

# %%
# Random_state seed
RS = 1234

# %%
# 1. Obtain parameters for the different analyses ----------------------------------------
par_file = sys.argv[1]
row_number = int(sys.argv[2])

# par_file = "AroAdm_kerPC.csv"
# par_file = "AroAdm_top-markers-10k.csv"
# row_number = 2

# Load grid
par_grid = pd.read_csv(os.path.join("parameters", par_file))
par_grid

# %%
# Extract parameters =====================================================================

# Which row of the grid
pars = par_grid.iloc[row_number]
print(pars, "\n")

# Model type: top-markers-10k or kerPC
model_type = pars["model_type"]

# Common pars
architecture = pars["architecture"]
m_set = pars["m_set"]
traits = pars["traits"]
part = pars["part"]
input_type = pars["input_type"]
output_type = pars["output_type"]

# Output folder
out_folder = os.path.join("models", f"{model_type}_{architecture}_{m_set}_{input_type}", traits, part)
if not os.path.isdir(out_folder):
    os.makedirs(out_folder)
print(out_folder)

# Save traits as list
# If multi-output, change literal string to list
if output_type == "multi-output":
    traits = ast.literal_eval(traits)
else:
    traits = [traits]
print(traits)

# Target variable type
# Save as list, just like traits
target_dict = {
    "culm.diameter.1st.internode": "binary",
    "leaf.senescence": "binary",
    "grain.weight": "continuous",
    "time.to.flowering.from.sowing": "continuous"
}

target_types = []
for trait_i in traits:
    if trait_i not in target_dict.keys():
        print("Define variable type of the target in the dictionary")
        exit(1)
    target_types.append(target_dict[trait_i])

print(target_types)

# %%
# 3. Prepare data ------------------------------------------------------------------------

# Features ===============================================================================
# Generate X_list, even if single-input

# Top-markers-10k
if model_type == "top-markers-10k":
    X_list = []
    X_scaler_list = []
    
    # For every trait, 10k most associated markers
    for i in range(len(traits)):
        X_name = f"geno_top-10000-{m_set}_{part}_{traits[i]}.csv"
        # Load and remove accession column
        X_df = pd.read_csv(os.path.join("data", model_type, X_name))
        X_df = X_df.iloc[:, 1:]
        
        X = X_df.to_numpy()
        # Scale using all data because genomic prediction
        scaler_X_all = StandardScaler()
        scaler_X_all.fit(X)
        X = scaler_X_all.transform(X)
        
        X_list.append(X)
        X_scaler_list.append(scaler_X_all)
        
# kerPC
elif model_type == "kerPC":
    # Save initial list with kerPCs based on marker set
    X_list_initial = []
    
    # Define marker sets
    if m_set == "SNP":
        marker_sets = ["SNP"]
    elif m_set == "all":
        marker_sets = ["DEL", "DUP", "INV", "MITE-DTX", "RLX-RIX", "SNP"]
    
    # Load matrices
    for i in range(len(marker_sets)):
        X_name = f"PC_{marker_sets[i]}.csv"
        X_df = pd.read_csv(os.path.join("data", model_type, X_name))
        X_df = X_df.iloc[:, 1:]
        
        X = X_df.to_numpy()
        X_list_initial.append(X)
    
    # Process depending on input-type
    X_list = []
    X_scaler_list = []
    
    # If single-input, concatenate matrices of list
    if input_type == "single-input":
        for i in range(len(X_list_initial)):
            if i == 0:
                X = X_list_initial[i]
                continue
            
            X = np.concatenate((X, X_list_initial[i]), axis = 1)
        
        # Scale using all data because genomic prediction
        scaler_X_all = StandardScaler()
        scaler_X_all.fit(X)
        X = scaler_X_all.transform(X)
        
        X_list.append(X)
        X_scaler_list.append(scaler_X_all)
    # Multi-input
    else:
        for i in range(len(X_list_initial)):
            X = X_list_initial[i]
            
            # Scale using all data because genomic prediction
            scaler_X_all = StandardScaler()
            scaler_X_all.fit(X)
            X = scaler_X_all.transform(X)
            
            X_list.append(X)
            X_scaler_list.append(scaler_X_all)

print(len(X_list))
print(X_list[0].shape)

# %%
# Target =================================================================================

Y_df = pd.read_csv(os.path.join("data", "pheno_original.csv"))
parts = pd.read_csv(os.path.join("data", "partitions.csv"))

# Obtain train and test masks ############################################################
# Single-output
if output_type == "single-output":
    y_df = Y_df[traits[0]]
    y = y_df.to_numpy().reshape(-1, 1)

    train_mask = (parts[part] == "train") & pd.notna(y_df)
    test_mask = (parts[part] == "test") & pd.notna(y_df)
# Multi-output
else:
    # Obtain combined mask of missing values for the different traits (final_NA_mask)
    NA_mask_list = []
    for i in range(len(traits)):
        y_df = Y_df[traits[i]]
        
        NA_mask = pd.notna(y_df)
        NA_mask_list.append(NA_mask)
    
    for i in range(len(NA_mask_list)):
        if i == 0:
            final_NA_mask = NA_mask_list[i]
            continue
        final_NA_mask = final_NA_mask & NA_mask_list[i]
        # print(initial_mask_list[i].sum())
        print(final_NA_mask.sum())
    
    train_mask = (parts[part] == "train") & final_NA_mask
    test_mask = (parts[part] == "test") & final_NA_mask
    
print(train_mask.sum(), test_mask.sum())

# Train / test split + train_sub / validation split ######################################

# Dicts of y splits for each trait (output)
y_train = {}; y_test = {}; y_train_sub = {}; y_val = {}
# y_train_scaler_dict = {}

for i in range(len(traits)):
    y_df = Y_df[traits[i]]
    y = y_df.to_numpy().reshape(-1, 1)
    
    key = traits[i]
    
    # 1. Target split
    
    y_train_i = y[train_mask, ]
    y_test_i = y[test_mask, ]
    
    # 2. Target process

    # Scale continuous targets (just with train data)
    # We no longer scale targets
    if target_types[i] == "continuous":
        # scaler_y_train = StandardScaler()
        # scaler_y_train.fit(y_train_i)
        # y_train_i = scaler_y_train.transform(y_train_i)
        # y_test_i = scaler_y_train.transform(y_test_i)
        
        # y_train_scaler_dict[key] = scaler_y_train
        pass
    # Recode binary targets from 1/2 to 0/1
    elif target_types[i] == "binary":
        y_train_i = np.where(y_train_i == 2, 1, 0)
        y_test_i = np.where(y_test_i == 2, 1, 0)
    else:
        print("Unknown target variable type")
        exit(1)
    
    y_train[key] = y_train_i
    y_test[key] = y_test_i
    
    # 3. Features split
    
    X_train = []
    X_test = []
    
    X_train_sub = []
    X_val = []
    
    # For every input
    for j in range(len(X_list)):
        X_train.append(X_list[j][train_mask, ])
        X_test.append(X_list[j][test_mask, ])
        
        # RS seed is important in this case to have the same splits
        X_train_sub_i, X_val_i, y_train_sub_i, y_val_i = train_test_split(
            X_train[j], y_train_i,
            test_size = 0.2,
            random_state = RS
        )
        X_train_sub.append(X_train_sub_i)
        X_val.append(X_val_i)
    
    y_train_sub[key] = y_train_sub_i
    y_val[key] = y_val_i

# X_lists and y_dicts
print(len(X_train), len(y_train.keys()))

print(X_train[0].shape, X_test[0].shape, X_train_sub[0].shape, X_val[0].shape)
print(y_train[traits[0]].shape, y_test[traits[0]].shape,
      y_train_sub[traits[0]].shape, y_val[traits[0]].shape)

# %%
# 4. Load hypermodels and create corresponding object ------------------------------------

from hypermodels import hyper_MLP
from hypermodels import hyper_CNN

# %%
input_dim = (X_train[0].shape[1], )
n_inputs = len(X_train)
n_outputs = len(traits)
print(input_dim, traits, target_types, n_inputs, n_outputs)


if architecture == "MLP":
    hypermodel = hyper_MLP(
        input_dim = input_dim, targets = traits, target_types = target_types,
        n_inputs = n_inputs, n_outputs = n_outputs
    )
elif architecture == "CNN":
    hypermodel = hyper_CNN(
        input_dim = input_dim, targets = traits, target_types = target_types,
        n_inputs = n_inputs, n_outputs = n_outputs
    )

# %%
# model = hypermodel.build(kt.HyperParameters())
# model.summary()
# history = model.fit(x = X_train, y = y_train, epochs = 10, validation_data = (X_val, y_val))
# model.evaluate(x = X_test, y = y_test, return_dict = True)

# %%
# 5. Prediction --------------------------------------------------------------------------

# Load best hpars ========================================================================

# Do prediction with the best hpars of current partition or AroAdm partition
print(which_hps)

trait_path = os.path.dirname(out_folder)
hps_path = os.path.join(trait_path, which_hps, "top_hps")
print(hps_path)

# Load dictionary of best hpars
with open(os.path.join(hps_path, f"top_hps_{1}.json"), "r") as infile:
    hps_dict = json.load(infile)
print(hps_dict)

# Create hpar object from dictionary
best_hps = kt.HyperParameters()
best_hps.values = hps_dict

# Results path: different folder than hypermodel results #################################
pred_path = os.path.join(trait_path, part, f"results_hps_{which_hps}")
if not os.path.isdir(pred_path):
    os.mkdir(pred_path)
##########################################################################################

# %%
plt.style.use("seaborn-whitegrid")
sns.set_style("whitegrid")

# Do prediction n_iter times =============================================================
n_iter = 5

# Results list (of Series)
res_list = []
# Prediction dictionary (by trait) which contains lists (of dataframes)
prediction_dict = {}
for i in range(len(traits)):
    prediction_dict[traits[i]] = []


# For each iteration
for iter in range(n_iter):

    # Create model from hypermodel class #################################################
    model = hypermodel.build(best_hps)
    # model.summary()
    
    # Training with early stopping. Use validation set ###################################
    callbacks = [
        keras.callbacks.EarlyStopping(
            monitor = "val_loss",
            patience = 10,
            min_delta = 0.001,
            restore_best_weights = True
        )
    ]

    history = hypermodel.fit(
        best_hps, model,
        x = X_train_sub, y = y_train_sub,
        epochs = 50,
        callbacks = callbacks,
        validation_data = (X_val, y_val),
        verbose = 2
    )

    # Training plot
    plt.plot(history.history["loss"], label = "loss")
    plt.plot(history.history["val_loss"], label = "val_loss")
    plt.legend(loc = "upper right")
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.savefig(os.path.join(pred_path, f"training_loss_plot_{iter}.png"), dpi = 300)
    # plt.show()
    plt.clf()

    # Save prediction results ############################################################
    
    # Compute metrics and prediction
    eval_result = model.evaluate(x = X_test, y = y_test, return_dict = True)
    pred_result = model.predict(X_test)

    # Series of results
    res = pd.Series(eval_result)

    # Dataset size
    size_train = y_train[traits[0]].shape[0]
    size_test = y_test[traits[0]].shape[0]
    size_total = size_train + size_test

    res["size_train"] = size_train
    res["size_test"] = size_test
    res["size_total"] = size_total

    # Compute extra metrics, depending on target type

    for i in range(len(target_types)):
        y_test_i = y_test[traits[i]]
        
        if target_types[i] == "continuous":
            # Correlation
            if output_type == "single-output":
                pred_result_i = pred_result
            else:
                pred_result_i = pred_result[i]
            
            y_hat = pred_result_i
            cor = np.corrcoef(y_test_i, y_hat, rowvar = False)[0, 1]
            res[f"{traits[i]}_cor"] = cor
            
        elif target_types[i] == "binary":
            # AUC
            if output_type == "single-output":
                pred_result_i = pred_result
            else:
                pred_result_i = pred_result[i]
            m = keras.metrics.AUC()
            m.update_state(y_test_i, pred_result_i)
            AUC = m.result().numpy()
            res[f"{traits[i]}_AUC"] = AUC
            
        # Save prediction just in case
        prediction = pd.DataFrame({
            "original_index": test_mask.index[test_mask == True],
            "prediction": pred_result_i.squeeze() 
        })
        # Prediction to dict (in list)
        prediction_dict[traits[i]].append(prediction)
    
    # Results to list
    res_list.append(res)

# %%
# Compute average results from res_list ==================================================

# Join res in dataframe
res_df = pd.DataFrame(columns = res_list[0].index)

for i in range(len(res_list)):
    series = res_list[i]
    res_df.loc[i, ] = series

# Compute mean and sd
RES = pd.DataFrame({
    "mean": res_df.mean(),
    "sd": res_df.std()
})
print(RES)

RES.to_csv(os.path.join(pred_path, "results.csv"), index = True, header = True)

# %%
# Compute average prediction from prediction_dict ========================================

# For each trait
for i in range(len(traits)):
    # Get list, and join predictions into dataframe
    prediction_list = prediction_dict[traits[i]]
    prediction_df = pd.DataFrame(columns = prediction_list[0]["original_index"])

    for j in range(len(prediction_list)):
        series = prediction_list[j]
        series.index = series["original_index"]
        series = series["prediction"]
        prediction_df.loc[j, ] = series
    
    # Compute mean and sd
    PREDICTION = pd.DataFrame({
        "original_index": prediction_df.columns,
        "mean": prediction_df.mean(),
        "sd": prediction_df.std()
    })
    # print(PREDICTION)
    
    PREDICTION.to_csv(os.path.join(pred_path, f"prediction_{traits[i]}.csv"),
                      index = False, header = True)
    
    
    # Make plots of predictions
    y_test_i = y_test[traits[i]]
        
    if target_types[i] == "continuous":
        
        pred_result_i = PREDICTION["mean"]
        y_hat = pred_result_i
        
        # Correlation plot
        plot = sns.scatterplot(x = y_test_i.squeeze(), y = y_hat.squeeze())
        plot.set(xlabel = "y_test", ylabel = "y_hat", title = traits[i])
        plt.savefig(os.path.join(pred_path, f"cor_plot_{traits[i]}.png"), dpi = 300)
        # plt.show()
        plt.clf()
    
    elif target_types[i] == "binary":
        
        pred_result_i = PREDICTION["mean"]
        
        # Confusion matrix plot
        y_hat = np.where(pred_result_i > 0.5, 1, 0)
        conf = confusion_matrix(y_test_i, y_hat)
        plot = sns.heatmap(conf, annot = True, cmap = "Blues")
        plot.set(xlabel = "y_hat", ylabel = "y_test", title = traits[i])
        plt.savefig(os.path.join(pred_path, f"conf_matrix_{traits[i]}.png"), dpi = 300)
        # plt.show()
        plt.clf()

# %%



