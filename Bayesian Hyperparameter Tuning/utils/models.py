# This file creates and returns ML models

import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, BatchNormalization, Dropout, Dense, Activation
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras import activations

def mlp_model(input_shape,
              n_layers,
              n_neurons=None,
              activation='relu', # Nothing else is supported now
              dropout=False,
              dropout_rate=0.2,
              batchnorm=False,
              opt=Adam(learning_rate=0.001),):
    '''
    Creates and return an MLP model based on the input parameters
    
    Inputs:
    input_shape = shape of the input
    n_layers = # layers
    n_neurons = list of # neurons in each layer, uses 64 neurons in each layer by default
    activation = activation function for hidden layers
    opt = optimizer with all necessary parameters included
    dropout = flag of whether to use dropout
    dropout_rate = dropout rate to use if using dropout, ignored otherwise
    batchnorm = flag of whether to use batch normalization
    
    Outputs:
    model = an MLP model using keras
    '''
    inputs = Input(shape=input_shape)
    if n_neurons == None:
        n_neurons = [64]*n_layers
    
    layer = inputs
    for n in n_neurons:
        # FC layer
        layer = Dense(n, activation)(layer)
        # batch norm
        if batchnorm:
            layer = BatchNormalization()(layer)
        # Activation
        layer = Activation(activations.relu)(layer)
        # dropout
        if dropout:
            layer = Dropout(dropout_rate)(layer)
    
    outputs = Dense(1, activation='sigmoid')(layer)
    
    model = Model(inputs, outputs)
    model.compile(optimizer=opt, loss='binary_crossentropy', metrics=[tf.keras.metrics.AUC(), 'acc'])
    
    return model