from sgrna_modeler import features as fe
from sklearn.model_selection import train_test_split
from sklearn import ensemble
from tensorflow import keras as k
import pandas as pd
import os
from joblib import load
import sgrna_modeler.enzymes as en

def curr_path():
    return os.path.dirname(__file__)

def get_deepcpf1_weights():
    path = os.path.join(curr_path(), 'data/saved_models/Seq_deepCpf1_weights_tf.h5')
    return path

def get_enpam_gb():
    path = os.path.join(curr_path(), 'data/saved_models/enPAM_GB.joblib')
    return path

def build_kim2018(input_shape=(34, 4)):
    """
    Build a convolutional neural network

    From:
    Kim, Hui Kwon, et al. "Deep learning improves prediction of CRISPR–Cpf1 guide RNA activity." \
    Nature biotechnology 36.3 (2018): 239.

    :param input_shape: guide length by nts (4)
    :type input_shape: tuple
    :return: CNN architecture
    :rtype: keras Model object
    """
    """Build a Convolutional neural network model from Kim 2018

    Parmeters
    ---------
    input_shape: tuple, optional (default (34, 4)
        shape of the first layer of the model

    Returns
    -------
    model: keras model object
    """
    Input_SEQ = k.layers.Input(shape=input_shape)
    C1 = k.layers.Convolution1D(80, 5, activation='relu')(Input_SEQ)
    P1 = k.layers.AveragePooling1D(2)(C1)
    F = k.layers.Flatten()(P1)
    DO1 = k.layers.Dropout(0.3)(F)
    D1 = k.layers.Dense(80, activation='relu')(DO1)
    DO2 = k.layers.Dropout(0.3)(D1)
    D2 = k.layers.Dense(40, activation='relu')(DO2)
    DO3 = k.layers.Dropout(0.3)(D2)
    D3 = k.layers.Dense(40, activation='relu')(DO3)
    DO4 = k.layers.Dropout(0.3)(D3)
    Output = k.layers.Dense(1, activation='linear')(DO4)
    model = k.models.Model(inputs = Input_SEQ, outputs = Output)
    return model

class KerasSgrnaModel(object):
    """This class is for creating, training, and predicting guide activity with a Keras model

    :param random_state: set random state in train/test split for reproducibility
    :type random_stat: int
    :param val_frac: amount of data to use for early stopping
    :type val_frac: float
    :param base_arc: base architecture to build neural network, defaults to build_kim2018
    :type base_arc: function, which takes an input shape and returns a keras model

    :Example:

    >>> from sgrna_modeler import datasets as da
    >>> from sgrna_modeler import models as sg
    >>> train_data = da.load_kim_2018_train()
    >>> train_model = sg.KerasSgrnaModel()
    >>> train_model.fit(train_data)
    >>> test_data = da.load_kim_2018_test()
    >>> test_predictions = train_model.predict(test_data)
    """
    def __init__(self, random_state = 7, val_frac = 0.1, base_arc = None):
        """Constructor
        """
        self.base_name = 'Keras_CNN'
        self.val_frac = val_frac
        self.random_state = random_state
        if base_arc is None:
            self.base_arc = build_kim2018
        else:
            self.base_arc = base_arc
        self.train_dataset = None
        self.enzyme = None
        self.model = None
        self.model_history = None
        self.train_name = None

    def load_weights(self, weights, enzyme, name):
        """Load previously trained weights

        :param enzyme: cas9 or cas12a
        :type enyme: dict
        :param weights: filepath to weights
        :type weights: str
        :param name: name of the model
        :type name:str
        """
        if weights is None:
            weights = get_deepcpf1_weights()
            self.train_name = 'Seq-DeepCpf1'
            self.enzyme = en.cas12a
        else:
            self.train_name = name
            self.enzyme = enzyme
        model = self.base_arc(input_shape = (self.enzyme['context_length'],4))
        model.load_weights(weights)
        self.model = model
        return self

    def fit(self, train_dataset):
        """ Fit a model to the training data

        :param train_dataset: training data
        :type train_dataset: :class:`sgrna_modeler.datasets.ActivityData`
        :return: self
        """
        self.train_dataset = train_dataset
        self.train_name = train_dataset.name
        self.enzyme = train_dataset.enzyme
        train_val_x, y = train_dataset.get_xy()
        encoded_train_val_x = fe.encode_seqs(train_val_x)
        train_x, val_x, train_y, val_y = train_test_split(encoded_train_val_x, y, test_size=self.val_frac,
                                                          random_state=self.random_state)
        model = self.base_arc(input_shape = (self.enzyme['context_length'],4))
        model.compile(optimizer='RMSprop',loss='mse',metrics=['mae'])
        self.model_history = model.fit(train_x, train_y, epochs = 200,
                                       validation_data = (val_x, val_y),
                                       callbacks = [k.callbacks.EarlyStopping(patience=20,restore_best_weights=True),
                                                    k.callbacks.History()],
                                       verbose = 0)
        self.model = model
        return self

    def predict(self, test_dataset):
        """Predict activity of test data

        :param test_dataset: testing data
        :type test_dataset: :class:`sgrna_modeler.datasets.ActivityData`
        :return: dataframe of predictions and other meta information
        :rtype: pandas dataframe
        """
        x, y = test_dataset.get_xy()
        encoded_x = fe.encode_seqs(x)
        predictions = self.model.predict(encoded_x)
        out_data = pd.DataFrame({'kmer': x, 'y': y})
        if test_dataset.group_column:
            out_data['group'] = test_dataset.data[test_dataset.group_column]
        else:
            out_data['group'] = ''
        out_data['prediction'] = predictions
        out_data['model'] = self.base_name
        out_data['training_data'] = self.train_name
        out_data['test_data'] = test_dataset.name
        return out_data

    def predict_seqs(self, seqs):
        """ Predict from sequences

        :param seqs: sequences to predict
        :return: numeric vector of predcitions
        """
        featurized_x = fe.encode_seqs(seqs)
        predictions = self.model.predict(featurized_x).flatten()
        return predictions

class SklearnSgrnaModel(object):
    """scikit-learn gradient boosting for modeling sgRNA activity

    :param random_state: set random state in train/test split for reproducibility
    :type random_state: int
    :param val_frac: amount of data to use for early stopping
    :type val_frac: float
    :param model: base model
    :type model: sklearn GradientBoostingRegressor
    :param features: features to model
    :type features: list

    :Example:
    >>> from sgrna_modeler import datasets as da
    >>> from sgrna_modeler import models as sg
    >>> train_model = sg.SklearnSgrnaModel()
    >>> rs2_data = da.load_doench_2016()
    >>> train_model.fit(rs2_data)
    """
    def __init__(self, random_state = 7, val_frac = 0.1, model = None, features = None):
        """Constructor
        """
        self.base_name = 'Sklearn_GB'
        self.val_frac = val_frac
        self.random_state = random_state
        if model is None:
            # Gradient boosted model
            self.model = ensemble.GradientBoostingRegressor(n_iter_no_change=20,
                                                            validation_fraction = self.val_frac,
                                                            random_state=self.random_state)
        else:
            self.model = model
        if features is None:
            # Default features for RuleSet2
            self.features = ['Pos. Ind. 1mer', 'Pos. Ind. 2mer', 'Pos. Dep. 1mer', 'Pos. Dep. 2mer', 'GC content', 'Tm']
        else:
            self.features = features
        self.enzyme = None
        self.train_dataset = None
        self.train_name = None

    def load_model(self, model, enzyme, name):
        """Load previously trained model

        :param enzyme: cas9 or cas12a
        :type enyme: dict
        :param model: filepath to trained model
        :type model: str (*.joblib)
        :param name: name of the model
        :type name:str
        """
        self.enzyme = enzyme
        self.model = load(model)
        self.train_name = name
        return self

    def fit(self, train_dataset):
        """ Fit a model to the training data

        :param train_dataset: training data
        :type train_dataset: :class:`sgrna_modeler.datasets.ActivityData`
        :return: self
        """
        self.train_name = train_dataset.name
        self.enzyme = train_dataset.enzyme
        train_val_x, y = train_dataset.get_xy()
        featurized_train_val_x = fe.featurize_guides(train_val_x, features=self.features,
                                                     guide_start = self.enzyme['guide_start'],
                                                     guide_length = self.enzyme['guide_length'])
        self.model.fit(featurized_train_val_x, y)
        return self

    def predict(self, test_dataset):
        """Predict activity of test data

        :param test_dataset: testing data
        :type test_dataset: :class:`sgrna_modeler.datasets.ActivityData`
        :return: dataframe of predictions and other meta information
        :rtype: pandas dataframe
        """
        x, y = test_dataset.get_xy()
        featurized_x = fe.featurize_guides(x, features=self.features,
                                           guide_start=test_dataset.enzyme['guide_start'],
                                           guide_length=test_dataset.enzyme['guide_length'])
        predictions = self.model.predict(featurized_x)
        out_data = pd.DataFrame({'kmer': x, 'y': y})
        if test_dataset.group_column:
            out_data['group'] = test_dataset.data[test_dataset.group_column]
        else:
            out_data['group'] = ''
        out_data['prediction'] = predictions
        out_data['model'] = self.base_name
        out_data['training_data'] = self.train_name
        out_data['test_data'] = test_dataset.name
        return out_data

    def predict_seqs(self, seqs):
        """ Predict from sequences

        :param seqs: sequences to predict
        :return: numeric vector of predcitions
        """
        featurized_x = fe.featurize_guides(seqs, features=self.features,
                                           guide_start=self.enzyme['guide_start'],
                                           guide_length=self.enzyme['guide_length'])
        predictions = self.model.predict(featurized_x)
        return predictions

