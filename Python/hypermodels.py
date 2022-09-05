import tensorflow as tf
from tensorflow import keras
import keras_tuner as kt
from keras import layers

class hyper_MLP(kt.HyperModel):
    def __init__(self, input_dim, targets, target_types, n_inputs, n_outputs):
        self.input_dim = input_dim
        self.targets = targets
        self.target_types = target_types
        self.n_inputs = n_inputs
        self.n_outputs = n_outputs
    
    def build(self, hp):
        # Same activation function for all dense layers
        activation = hp.Choice("activation", values = ["relu", "tanh", "linear", "softplus"])
        # Same regularization for all dense layers
        reg_type = hp.Choice("regularization_type", values = ["L2", "L1"])
        reg_rate = hp.Choice("regularization_rate", values = [0.0, 0.1, 0.01, 0.001])
        
        if reg_type == "L2":
            kernel_regularizer = keras.regularizers.L2(l2 = reg_rate)
        elif reg_type == "L1":
            kernel_regularizer = keras.regularizers.L1(l1 = reg_rate)
        
        # Input layer ####################################################################
        
        if self.n_inputs == 1:
            inputs = layers.Input(shape = self.input_dim)
            x = inputs
        # Multi-input
        else:
            inputs = []
            # First dense layer for each input
            first_multi_layer = []
            
            first_units = hp.Choice("first_units", values = [10, 20, 40, 80, 160])
            for i in range(self.n_inputs):
                inputs.append(layers.Input(shape = self.input_dim))
                
                first_multi_layer.append(
                    layers.Dense(units = first_units, activation = activation,
                                 kernel_regularizer = kernel_regularizer)(inputs[i])
                )
            
            # Concatenate
            x = layers.Concatenate(axis = 1)(first_multi_layer)
        
        # Body: Dense layers #############################################################
        for i in range(hp.Int("n_layers", min_value = 1, max_value = 5, step = 1)):
            x = layers.Dense(
                    units = hp.Choice(f"units_{i}", values = [10, 20, 40, 80, 160]),
                    activation = activation, kernel_regularizer = kernel_regularizer
                )(x)
        
        # Dropout layer ##################################################################
        x = layers.Dropout(
            rate = hp.Float("dropout_rate", min_value = 0.0, max_value = 0.3, step = 0.05)
        )(x)
        
        
        # Output layer ###################################################################
        
        # Get dictionaries of losses, metrics and activations for each output
        losses = {}
        metrics_dict = {}
        out_activation_dict = {}
        
        for i in range(self.n_outputs):
            if self.target_types[i] == "continuous":
                    out_activation = "linear"
                    loss = "mse"
                    metrics = ["mse"]
            elif self.target_types[i] == "binary":
                out_activation = "sigmoid"
                loss = "binary_crossentropy"
                metrics = ["binary_crossentropy", "accuracy"]
                
            losses[self.targets[i]] = loss
            metrics_dict[self.targets[i]] = metrics
            out_activation_dict[self.targets[i]] = out_activation
        
        # Add layers
        if self.n_outputs == 1:
            outputs = layers.Dense(
                units = 1,
                activation = out_activation_dict[self.targets[0]],
                name = self.targets[0]
            )(x)
        # Multi-output
        else:
            # Dense layer for each output
            preout_multi_layer = []
            outputs = []
            
            last_units = hp.Choice("last_units", values = [10, 20, 40, 80, 160])
            
            for i in range(self.n_outputs):
                preout_multi_layer.append(
                    layers.Dense(units = last_units, activation = activation,
                                 kernel_regularizer = kernel_regularizer)(x)
                )
                
                outputs.append(
                    layers.Dense(
                        units = 1,
                        activation = out_activation_dict[self.targets[i]],
                        name = self.targets[i]
                    )(preout_multi_layer[i])
                )
        
        ##################################################################################
        
        # Create model
        model = keras.Model(inputs = inputs, outputs = outputs)
        
        # Compile
        # optimizer = hp.Choice("optimizer", values = ["Adam", "RMSprop", "SGD"])
        optimizer = "Adam"
        model.compile(
            optimizer = optimizer,
            loss = losses,
            metrics = metrics_dict
        )
        
        return model
        
    def fit(self, hp, model, *args, **kwargs):
        return model.fit(
            *args,
            **kwargs
        )

# ----------------------------------------------------------------------------------------

class hyper_CNN(kt.HyperModel):
    def __init__(self, input_dim, targets, target_types, n_inputs, n_outputs):
        self.input_dim = input_dim
        self.targets = targets
        self.target_types = target_types
        self.n_inputs = n_inputs
        self.n_outputs = n_outputs
    
    def build(self, hp):
        # Same activation function for all dense layers
        activation = hp.Choice("activation", values = ["relu", "tanh", "linear", "softplus"])
        # Same regularization for all dense and convolutional layers
        reg_type = hp.Choice("regularization_type", values = ["L2", "L1"])
        reg_rate = hp.Choice("regularization_rate", values = [0.0, 0.1, 0.01, 0.001])
        
        if reg_type == "L2":
            kernel_regularizer = keras.regularizers.L2(l2 = reg_rate)
        elif reg_type == "L1":
            kernel_regularizer = keras.regularizers.L1(l1 = reg_rate)
        
        # Input layer ####################################################################
        
        if self.n_inputs == 1:
            inputs = layers.Input(shape = self.input_dim)
            
            # Conv1D input: batch (steps, channels) -> None (10000, 1)
            x = layers.Reshape(target_shape = (self.input_dim[0], 1))(inputs)
            x = layers.Conv1D(
                filters = hp.Choice("n_filters", values = [16, 32, 64, 128]),
                kernel_size = 3,
                activation = activation,
                kernel_regularizer = kernel_regularizer
            )(x)
            # Output: batch (new_steps, filters)
            
            # Pooling
            x = layers.MaxPooling1D(pool_size = 3)(x)
            
            # Flatten output to 1 dim
            x = layers.Flatten()(x)
            
        # Multi-input
        else:
            inputs = []
            pre_concat = []
            
            for i in range(self.n_inputs):
                inputs.append(layers.Input(shape = self.input_dim))
                
                # Conv1D input: batch (steps, channels) -> None (10000, 1)
                pre_concat.append(
                    layers.Reshape(target_shape = (self.input_dim[0], 1))(inputs[i])
                )
                pre_concat[i] = layers.Conv1D(
                    filters = hp.Choice("n_filters", values = [16, 32, 64, 128]),
                    kernel_size = 3,
                    activation = activation,
                    kernel_regularizer = kernel_regularizer
                    )(pre_concat[i])
                    # Output: batch (new_steps, filters)
                
                # Pooling
                pre_concat[i] = layers.MaxPooling1D(pool_size = 3)(pre_concat[i])
                
                # Flatten output to 1 dim
                pre_concat[i] = layers.Flatten()(pre_concat[i])
                
                # Because multi-input, add additional dense layer before concatenating
                pre_concat[i] = layers.Dense(
                    units = hp.Choice("first_units", values = [10, 20, 40, 80, 160]),
                    activation = activation,
                    kernel_regularizer = kernel_regularizer
                )(pre_concat[i])
            
            # Concatenate 
            x = layers.Concatenate(axis = 1)(pre_concat)
        
        # Body: Dense layers #############################################################
        for i in range(hp.Int("n_layers", min_value = 1, max_value = 5, step = 1)):
            x = layers.Dense(
                    units = hp.Choice(f"units_{i}", values = [10, 20, 40, 80, 160]),
                    activation = activation,
                    kernel_regularizer = kernel_regularizer
                )(x)
        
        # Dropout layer ##################################################################
        x = layers.Dropout(
            rate = hp.Float("dropout_rate", min_value = 0.0, max_value = 0.3, step = 0.05)
        )(x)
        
        
        # Output layer ###################################################################
        
        # Get dictionaries of losses, metrics and activations for each output
        losses = {}
        metrics_dict = {}
        out_activation_dict = {}
        
        for i in range(self.n_outputs):
            if self.target_types[i] == "continuous":
                    out_activation = "linear"
                    loss = "mse"
                    metrics = ["mse"]
            elif self.target_types[i] == "binary":
                out_activation = "sigmoid"
                loss = "binary_crossentropy"
                metrics = ["binary_crossentropy", "accuracy"]
                
            losses[self.targets[i]] = loss
            metrics_dict[self.targets[i]] = metrics
            out_activation_dict[self.targets[i]] = out_activation
        
        # Add layers
        if self.n_outputs == 1:
            outputs = layers.Dense(
                units = 1,
                activation = out_activation_dict[self.targets[0]],
                name = self.targets[0]
            )(x)
        # Multi-output
        else:
            # Dense layer for each output
            preout_multi_layer = []
            outputs = []
            
            last_units = hp.Choice("last_units", values = [10, 20, 40, 80, 160])
            
            for i in range(self.n_outputs):
                preout_multi_layer.append(
                    layers.Dense(units = last_units, activation = activation,
                                 kernel_regularizer = kernel_regularizer)(x)
                )
                
                outputs.append(
                    layers.Dense(
                        units = 1,
                        activation = out_activation_dict[self.targets[i]],
                        name = self.targets[i]
                    )(preout_multi_layer[i])
                )
        
        ##################################################################################
        
        # Create model
        model = keras.Model(inputs = inputs, outputs = outputs)
        
        # Compile
        # optimizer = hp.Choice("optimizer", values = ["Adam", "RMSprop", "SGD"])
        optimizer = "Adam"
        model.compile(
            optimizer = optimizer,
            loss = losses,
            metrics = metrics_dict
        )
        
        return model
        
    def fit(self, hp, model, *args, **kwargs):
        return model.fit(
            *args,
            **kwargs
        )