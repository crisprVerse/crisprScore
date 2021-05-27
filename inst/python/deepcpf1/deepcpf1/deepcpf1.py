from numpy import *
import sys  
import os
import numpy as np

from os import environ
stderr = sys.stderr
sys.stderr = open(os.devnull, 'w')
from keras import backend as K
sys.stderr = stderr
from keras.models import Model
from keras.layers import Input
from keras.layers.merge import Multiply
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.layers.convolutional import Convolution1D, AveragePooling1D


dir_path = os.path.dirname( os.path.realpath( __file__ ) )
weights_path = os.path.join( dir_path, 'weights/Seq_deepCpf1_weights.h5' )


# user defined function to change keras backend
def set_keras_backend(backend):
    if K.backend() != backend:
       environ['KERAS_BACKEND'] = backend
       stderr = sys.stderr
       sys.stderr = open(os.devnull, 'w')
       reload(K)
       sys.stderr = stderr
       assert K.backend() == backend

# call the function with "theano"


set_keras_backend("theano")

def deepcpf1(sequences):
    if K.backend() != "theano":
        print "ERROR: Not using the theano backend. Check the requirements."
    Seq_deepCpf1_Input_SEQ = Input(shape=(34,4))
    Seq_deepCpf1_C1 = Convolution1D(80, 5, activation='relu')(Seq_deepCpf1_Input_SEQ)
    Seq_deepCpf1_P1 = AveragePooling1D(2)(Seq_deepCpf1_C1)
    Seq_deepCpf1_F = Flatten()(Seq_deepCpf1_P1)
    Seq_deepCpf1_DO1= Dropout(0.3)(Seq_deepCpf1_F)
    Seq_deepCpf1_D1 = Dense(80, activation='relu')(Seq_deepCpf1_DO1)
    Seq_deepCpf1_DO2= Dropout(0.3)(Seq_deepCpf1_D1)
    Seq_deepCpf1_D2 = Dense(40, activation='relu')(Seq_deepCpf1_DO2)
    Seq_deepCpf1_DO3= Dropout(0.3)(Seq_deepCpf1_D2)
    Seq_deepCpf1_D3 = Dense(40, activation='relu')(Seq_deepCpf1_DO3)
    Seq_deepCpf1_DO4= Dropout(0.3)(Seq_deepCpf1_D3)
    Seq_deepCpf1_Output = Dense(1, activation='linear')(Seq_deepCpf1_DO4)
    Seq_deepCpf1 = Model(inputs=[Seq_deepCpf1_Input_SEQ], outputs=[Seq_deepCpf1_Output])
    #Seq_deepCpf1.load_weights('weights/Seq_deepCpf1_weights.h5')
    Seq_deepCpf1.load_weights(weights_path)
    SEQ=preprocess_sequences(sequences)
    Seq_deepCpf1_SCORE = Seq_deepCpf1.predict([SEQ], batch_size=50, verbose=0)
    return Seq_deepCpf1_SCORE

def preprocess_sequences(sequences):
    data_n = sequences.shape[0]
    SEQ = zeros((data_n, 34, 4), dtype=int)
    
    for l in range(data_n):
        seq  = sequences[l]
        for i in range(34):
            if seq[i] in "Aa":
                SEQ[l, i, 0] = 1
            elif seq[i] in "Cc":
                SEQ[l, i, 1] = 1
            elif seq[i] in "Gg":
                SEQ[l, i, 2] = 1
            elif seq[i] in "Tt":
                SEQ[l, i, 3] = 1
    return SEQ