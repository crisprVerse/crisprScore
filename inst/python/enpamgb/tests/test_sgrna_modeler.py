#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `sgrna_modeler` package."""

import pandas as pd
from pandas.util.testing import assert_frame_equal
import os

import numpy as np
from sgrna_modeler import features as ft
from sgrna_modeler import datasets as da
from sgrna_modeler import models as sg
from sgrna_modeler import mutagenesis as mu
import sgrna_modeler.enzymes as en
from scipy.stats import stats

def curr_path():
    return os.path.dirname(__file__)

def test_encoding(seqs=None):
    if seqs is None:
        seqs = ['ACT', 'ACT', 'GCT', 'AAT', 'AAA']
    encoded = ft.encode_seqs(seqs)
    np.testing.assert_array_equal(encoded[3, :, :],
                                  np.array([[1, 0, 0, 0],
                                            [1, 0, 0, 0],
                                            [0, 0, 0, 1]]))

def test_context_order():
    assert ft.get_context_order(4) == ['1', '2', '3', '4']

def test_guide_seq():
    # SpCas9
    assert ft.get_guide_sequence('CGTCCCCATCCACGGCCTTCACCCGGGCAG', 5, 20) == 'CCCATCCACGGCCTTCACCC'
    # AsCas12a
    assert ft.get_guide_sequence('CCAGTTTGAACTCTCGCCCATCACCTATCAGTGC', 9, 20) == 'AACTCTCGCCCATCACCTAT'

def test_featurization():
    features = {'Pos. Ind. 1mer',
             'Pos. Ind. 2mer',
             'Pos. Ind. 3mer',
             'Pos. Ind. Zipper',
             'Pos. Dep. 1mer',
             'Pos. Dep. 2mer',
             'Pos. Dep. 3mer',
             'Pos. Dep. Zipper',
             'Pos. Ind. Rep.',
             'GC content',
             'Tm',
             'Cas9 PAM',
             'Physio',
             'OOF Mutation Rate',
             'Double Zipper'}

    kmers = pd.Series(['ACTGGTGGG'])
    one_hot = ft.featurize_guides(kmers, features, guide_start=3, guide_length=6).iloc[0]
    assert (one_hot['3T'] == 1)
    assert (one_hot['TGG'] == 0.5)
    assert (one_hot['GC content'] == 2/3)
    assert (one_hot['bendability'] != 0)
    assert (one_hot['Tm, context'] != 0)
    assert (one_hot['1AC'] == 1)
    assert (one_hot['7GGG'] == 1)
    assert (one_hot['G'] == 2/3)
    assert (one_hot['TG'] == 0.4)
    assert (one_hot['1T'] == 0)

def test_data_load():
    # fname = os.path.join(os.path.dirname(__file__), 'kim_2018_test.csv')
    # data = pd.read_csv(fname)
    # cas12a = en.Enzyme(guide_start=9, guide_length=23,
    #                    pam_start=5, pams = ['TTTN'], context_length=34)
    # kim_2018_train = da.Activity_Data(data, cas12a, 'Context Sequence',
    #                                   'Indel frequency')
    # x, y = kim_2018_train.get_xy()
    # assert x.shape == y.shape
    # assert kim_2018_train.enzyme.context_length == 34
    doench = da.load_doench_2016()
    x, y = doench.get_xy()
    assert np.isscalar(y[0])
    assert len(x[0]) == 30

    meyers_train = da.load_meyers_2017_train()
    x, y = meyers_train.get_xy()
    assert np.isscalar(y[0])
    assert len(x[0]) == 30

    meyers_test = da.load_meyers_2017_test()
    x, y = meyers_test.get_xy()
    assert np.isscalar(y[0])
    assert len(x[0]) == 30

    kim_2019_train = da.load_kim_2019_train()
    x, y = kim_2019_train.get_xy()
    assert np.isscalar(y[0])
    assert len(x[0]) == 30

    kim_2019_test = da.load_kim_2019_test()
    x, y = kim_2019_test.get_xy()
    assert np.isscalar(y[0])
    assert len(x[0]) == 30

    kim_2018_train = da.load_kim_2018_train()
    x, y = kim_2018_train.get_xy()
    assert np.isscalar(y[0])
    assert len(x[0]) == 34

    kim_2018_test = da.load_kim_2018_test()
    x, y = kim_2018_test.get_xy()
    assert np.isscalar(y[0])
    assert len(x[0]) == 34


def test_model_keras():
    # train a model
    train_data = da.load_kim_2018_train()
    train_model = sg.KerasSgrnaModel()
    train_model.fit(train_data)
    test_data = da.load_kim_2018_test()
    test_predictions = train_model.predict(test_data)
    np.testing.assert_almost_equal(stats.pearsonr(test_predictions.y, test_predictions.prediction)[0],0.75, decimal=2)

    # load a model
    load_model = sg.KerasSgrnaModel()
    load_model.load_weights(sg.get_deepcpf1_weights(), en.cas12a, 'Seq-DeepCpf1')
    load_predictions = load_model.predict(test_data)
    np.testing.assert_almost_equal(stats.pearsonr(load_predictions.y, load_predictions.prediction)[0], 0.75, decimal=2)

    # predict seqs
    x, y = test_data.get_xy()
    seq_predictions = load_model.predict_seqs(x)
    np.testing.assert_array_equal(seq_predictions, load_predictions['prediction'])

def load_kim_2019_rs2():
    data = pd.read_csv(os.path.join(curr_path(), 'test_data/kim_2019_rs2_predictions.csv.zip'))
    data_class = da.ActivityData(data = data, enzyme = en.cas9, kmer_column='Input String',
                                 activity_column='on-target',
                                 name = 'km_2019_rs2')
    return data_class

def test_model_sklearn():
    # train a model, compare with rs2
    train_model = sg.SklearnSgrnaModel()
    rs2_data = da.load_doench_2016()
    train_model.fit(rs2_data)
    kim_2019_rs2_predictions = load_kim_2019_rs2()
    train_predictions = train_model.predict(kim_2019_rs2_predictions)
    _, kim_rs2_y = kim_2019_rs2_predictions.get_xy()
    assert stats.pearsonr(train_predictions['prediction'], kim_rs2_y)[0] > 0.9

    # load a model
    enpam_gb = sg.get_enpam_gb()
    load_model = sg.SklearnSgrnaModel()
    load_model.load_model(enpam_gb, en.cas12a, 'enPAM_GB')
    kim_2018 = da.load_kim_2018_test()
    x, y = kim_2018.get_xy()
    predictions = load_model.predict_seqs(x)
    np.testing.assert_almost_equal(stats.pearsonr(y, predictions)[0], 0.57, decimal = 2)


def test_mutagenesis():
    load_model = sg.KerasSgrnaModel()
    load_model.load_weights(sg.get_deepcpf1_weights(), en.cas12a, 'Seq-DeepCpf1')
    deltas = mu.mutagenize_model(load_model, 500)
    delta_summaries = (deltas.groupby(by = ['nt', 'position'])
                       .agg({'delta': 'mean'})
                       .reset_index())
    np.array_equal(delta_summaries.nsmallest(1, 'delta')[['nt', 'position']].values,
                   pd.DataFrame({'nt': ['T'], 'position': [8]}).values)

