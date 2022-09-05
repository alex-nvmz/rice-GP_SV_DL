# Correlation metric evaluated over the entire dataset, not the average between mini-batches
# class PearsonCorrelation(keras.metrics.Metric):
#     def __init__(self, **kwargs):
#         super().__init__(name = "correlation", **kwargs)

# No checkpoint for Hpar search
    # keras.callbacks.ModelCheckpoint(
    #     filepath = os.path.join(out_folder, trait, "checkpoints", "checkpoint_{epoch}"),
    #     save_freq = "epoch"
    # )
    
    
# Old ------------------------------------------------------------------------------------
def build_model(hp):
    model = keras.Sequential()
    model.add(layers.Flatten())
    # Tune the number of layers.
    for i in range(hp.Int("num_layers", 1, 3)):
        model.add(
            layers.Dense(
                # Tune number of units separately.
                units=hp.Int(f"units_{i}", min_value=32, max_value=512, step=32),
                activation=hp.Choice("activation", ["relu", "tanh"]),
            )
        )
    if hp.Boolean("dropout"):
        model.add(layers.Dropout(rate=0.25))
    model.add(layers.Dense(10, activation="softmax"))
    learning_rate = hp.Float("lr", min_value=1e-4, max_value=1e-2, sampling="log")
    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=learning_rate),
        loss="categorical_crossentropy",
        metrics=["accuracy"],
    )
    return model


def build_model(hp):
    
    inputs = layers.Flatten()
    x = layers.Dense(
        units = hp.Int("units", min_value = 10, max_value = 100, step = 10),
        activation = "relu"
    )(inputs)
    if hp.Boolean("dropout"):
        x = layers.Dropout(rate = 0.25)
    
    outputs = layers.Dense(1, activation = "linear")(x)
    
    model = keras.Model(inputs = inputs, outputs = outputs)
    
    learning_rate = hp.Float("lr", min_value = 1e-4, max_value = 1e-2, sampling = "log")
    model.compile(optimizer = keras.optimizers.SGD(learning_rate = learning_rate),
                  loss = "mse",
                  metrics = ["mse"])
    return model

# Default MLP

n_features = X_train.shape[1]

inputs = keras.Input(shape = (n_features))
x = layers.Dense(64, activation = "relu")(inputs)
x = layers.Dense(32, activation = "softplus")(x)

outputs = layers.Dense(1, activation = "linear")(x)

model = keras.Model(inputs = inputs, outputs = outputs)

model.summary()

model.input_shape

class hyper_MLP(kt.HyperModel):
    def __init__(self, input_dim):
        self.input_dim = input_dim
    
    def build(self, hp):
        inputs = keras.Input(shape = (self.input_dim))
        x = layers.Dense(
            units = hp.Int("units", min_value = 10, max_value = 60, step = 10),
            activation = "relu"
        )(inputs)
        x = layers.Dense(32, activation = "softplus")(x)
        outputs = layers.Dense(1, activation = "linear")(x)
        
        model = keras.Model(inputs = inputs, outputs = outputs)
        model.compile(optimizer = keras.optimizers.Adam(),
                      loss = "mse",
                      metrics = ["accuracy", "mse"])
        return model
    
    def fit(self, hp, model, *args, **kwargs):
        return model.fit(
            *args,
            **kwargs
        )
        
# ----------------------------------------------------------------------------------------

# 3. Prepare data ------------------------------------------------------------------------

# Features
if model == "top-markers":
    X_name = f"geno_top-{n_mark}-{m_set}_{part}_{trait}.csv"
    # Load and remove accession column
    X_df = pd.read_csv(os.path.join("data", model, X_name))
    X_df = X_df.iloc[:, 1:]
    
    X = X_df.to_numpy()
    # Scale using all data because genomic prediction
    scaler_X_all = StandardScaler()
    scaler_X_all.fit(X)
    X = scaler_X_all.transform(X)
    
    print(X[0:5, 0:5])
    
elif model == "kerPC":
    # Save dictionaries with kerPCs for each marker set and scaler objects
    X_dict = {}
    X_scaler_dict = {}
    
    # Define marker sets
    if m_set == "SNP":
        marker_sets = ["SNP"]
    elif m_set == "all":
        marker_sets = ["DEL", "DUP", "INV", "MITE-DTX", "RLX-RIX", "SNP"]
    
    # Load matrices, scale and save in dict
    for marker_set in marker_sets:
        X_name = f"PC_{marker_set}.csv"
        X_df = pd.read_csv(os.path.join("data", model, X_name))
        X_df = X_df.iloc[:, 1:]
        
        X = X_df.to_numpy()
        # Scale using all data because genomic prediction
        scaler_X_all = StandardScaler()
        scaler_X_all.fit(X)
        X = scaler_X_all.transform(X)

        X_dict[marker_set] = X
        X_scaler_dict[marker_set] = scaler_X_all
            
    # If single-input, concatenate matrices of dictionary
    if input == "single-input":
        flag = True
        for marker_set in X_dict.keys():
            if flag:
                X = X_dict[marker_set]
                flag = False
                continue
            
            X_tmp = X_dict[marker_set]
            X = np.concatenate((X, X_tmp), axis = 1)
    elif input == "multi-input":
        print("Multi-input model, use X_dict")
        
# ----------------------------------------------------------------------------------------

# Target
Y_df = pd.read_csv(os.path.join("data", "pheno_original.csv"))
parts = pd.read_csv(os.path.join("data", "partitions.csv"))

y_df = Y_df[trait]
y = y_df.to_numpy().reshape(-1, 1)

# Train / test split
train_mask = (parts[part] == "train") & pd.notna(y_df)
test_mask = (parts[part] == "test") & pd.notna(y_df)

y_train = y[train_mask, ]
y_test = y[test_mask, ]

# Scale continuous targets (just with train data)
# We no longer scale targets
if target_type == "continuous":
    # scaler_y_train = StandardScaler()
    # scaler_y_train.fit(y_train)
    # y_train = scaler_y_train.transform(y_train)
    # y_test = scaler_y_train.transform(y_test)
    pass
# Recode binary targets from 1/2 to 0/1
elif target_type == "binary":
    y_train = np.where(y_train == 2, 1, 0)
    y_test = np.where(y_test == 2, 1, 0)
else:
    print("Unknown target variable type")
    exit(1)


# Split features

# Multi-input case
if model == "kerPC":
    if input == "multi-input":
        X_dict_train = {}
        X_dict_test = {}
        
        X_dict_train_sub = {}
        X_dict_val = {}
        for marker_set in X_dict.keys():
            X_dict_train[marker_set] = X_dict[marker_set][train_mask, ]
            X_dict_test[marker_set] = X_dict[marker_set][test_mask, ]
            
            # RS seed is important in this case
            X_train_sub, X_val, y_train_sub, y_val = train_test_split(
                X_dict_train[marker_set], y_train,
                test_size = 0.2,
                random_state = RS
            )
            X_dict_train_sub[marker_set] = X_train_sub
            X_dict_val[marker_set] = X_val
# Single-input
else:
    X_train = X[train_mask, ]
    X_test = X[test_mask, ]

    # Create validation set
    # For Hpar tuning and to decide best n epochs for prediction
    X_train_sub, X_val, y_train_sub, y_val = train_test_split(X_train, y_train,
                                                              test_size = 0.2,
                                                              random_state = RS)
    
# ----------------------------------------------------------------------------------------

# if output_type == "multi-output":
#     objective = []
#     for i in range(len(target_type)):
#         if target_type[i] == "continuous":
#             objective.append(f"val_mse")
#         elif target_type[i] == "binary":
#             objective.append(f"val_binary_crossentropy")

# else:
#     if target_type == "continuous":
#         objective = "val_mse"
#     elif target_type == "binary":
#         objective = "val_binary_crossentropy"

# ----
# Dicts of y splits for each trait (output)
    # y_list = []
    # for i in range(len(trait)):
    #     y_df = Y_df[trait[i]]
    #     y = y_df.to_numpy().reshape(-1, 1)
        
        # y_list.append(
        #     y_process_and_Xy_split(
        #         y = y, target_type_i = target_type[i], train_mask = train_mask, test_mask = test_mask
        #     )[4:]
        # )
    # Single-input, so only one split for X
    # X_train, X_test, X_train_sub, X_val = y_process_and_Xy_split(
    #     y = y, target_type_i = target_type[i], train_mask = train_mask, test_mask = test_mask
    # )[:4]

# ---
# # y_dict of output: y
# # INDEX 0: train, 1: test, 2: train_sub, 3: val
# def get_y_dict(index, y_list):
#     y_dict = {}
#     for i in range(len(y_list)):
#         y_dict[f"output_{i}"] = y_list[i][index]
#     return y_dict

# ----------------------------------------------------------------------------------------
# To load dictionary ###
with open(os.path.join(hps_path, f"top_hps_{1}.json"), "r") as infile:
    hps_dict = json.load(infile)
print(hps_dict)

# To create hp object to build hypermodel ###
in_hps = kt.HyperParameters()
in_hps.values = hps_dict
model = hypermodel.build(in_hps)
model.summary()

# To reload tuning results ###
in_tuner = kt.Hyperband(
    hypermodel,
    directory = out_folder,
    project_name = "tuning",
    overwrite = False
)
in_tuner.results_summary(1)