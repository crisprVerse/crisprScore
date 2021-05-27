import numpy as np
import pandas as pd
from Bio.SeqUtils import MeltingTemp
import re
import os.path

nt_codes = {'A':[1,0,0,0],
            'C':[0,1,0,0],
            'G':[0,0,1,0],
            'T':[0,0,0,1]}

def encode_seqs(seqs):
    # 3d array with samples x position x nt
    encoded_seqs = np.array([[nt_codes.get(x) for x in seq] for seq in seqs])
    return encoded_seqs

## Feature Engineering
def get_frac_g_or_c(dict, guide_sequence):
    """Get gc content
    :param context: sequence
    :param len: length of guide
    :param guide_start: position that guide starts
    :return: curr_dict with gc content
    """
    g_count = guide_sequence.count('G')
    c_count = guide_sequence.count('C')
    gc_frac = (g_count + c_count)/len(guide_sequence)
    dict['GC content'] = gc_frac
    return dict

def get_one_nt_counts(dict, guide, nts):
    for nt in nts:
        nt_frac = guide.count(nt)/len(guide)
        dict[nt] = nt_frac
    return dict

def get_two_nt_counts(dict, guide, nts):
    for nt1 in nts:
        for nt2 in nts:
            two_mer = nt1 + nt2
            nts_counts = guide.count(two_mer)
            nts_frac = nts_counts/(len(guide) - 1)
            dict[nt1 + nt2] = nts_frac
    return dict

def get_three_nt_counts(dict, guide, nts):
    for nt1 in nts:
        for nt2 in nts:
            for nt3 in nts:
                k_mer = nt1 + nt2 + nt3
                nts_counts = guide.count(k_mer)
                nts_frac = nts_counts/(len(guide) - 2)
                dict[nt1 + nt2 + nt3] = nts_frac
    return dict

def get_one_nt_pos(dict, context_sequence, nts, context_order):
    for i in range(len(context_order)):
        curr_nt = context_sequence[i]
        for nt in nts:
            key = context_order[i] + nt
            if curr_nt == nt:
                dict[key] = 1
            else:
                dict[key] = 0
    return dict

def get_two_nt_pos(dict, context_sequence, nts, context_order):
    for i in range(len(context_order) - 1):
        curr_nts = context_sequence[i:i+2]
        for nt1 in nts:
            for nt2 in nts:
                match_nts = nt1+nt2
                key = context_order[i] + match_nts
                if curr_nts == match_nts:
                    dict[key] = 1
                else:
                    dict[key] = 0
    return dict

def get_three_nt_pos(dict, context_sequence, nts, context_order):
    for i in range(len(context_order) - 2):
        curr_nts = context_sequence[i:i+3]
        for nt1 in nts:
            for nt2 in nts:
                for nt3 in nts:
                    match_nts = nt1+nt2+nt3
                    key = context_order[i] + match_nts
                    if curr_nts == match_nts:
                        dict[key] = 1
                    else:
                        dict[key] = 0
    return dict

def get_thermo(dict, guide_sequence, context_sequence):
    # Use Biopython to get thermo info. from context and guides
    dict['Tm, context'] = MeltingTemp.Tm_NN(context_sequence)
    third = len(guide_sequence)//3
    dict['Tm, start'] = MeltingTemp.Tm_NN(guide_sequence[0:third])
    dict['Tm, mid'] = MeltingTemp.Tm_NN(guide_sequence[third:2*third])
    dict['Tm, end'] = MeltingTemp.Tm_NN(guide_sequence[2*third:])
    return dict

def get_context_order(k):
    """
    :param k: length of kmer
    :return: list of characters of each nt position
    """
    context_order = [str(x + 1) for x in range(k)]
    return context_order

def get_guide_sequence(context, guide_start, guide_length):
    return context[(guide_start-1):(guide_start + guide_length -1)]

def get_physiochemical(curr_dict, guide, nts, physiochemical_data):
    nt_counts = np.array([guide.count(nt1+nt2) for nt1 in nts for nt2 in nts])
    no_dinucs = physiochemical_data.drop(columns=['Dinucleotides'])
    numeric_physio = np.transpose(np.array(no_dinucs))
    physio_sum = np.sum(nt_counts*numeric_physio, axis = 1)
    curr_dict = {**curr_dict, **dict(zip(no_dinucs.keys(), physio_sum))}
    return curr_dict

def get_zipper_pos(dict, context_sequence, nts, context_order):
    for i in range(len(context_order) - 2):
        curr_nts = context_sequence[i] + context_sequence[i+2]
        for nt1 in nts:
            for nt2 in nts:
                match_nts = nt1+nt2
                key = context_order[i] + nt1 + 'N' + nt2
                if curr_nts == match_nts:
                    dict[key] = 1
                else:
                    dict[key] = 0
    return dict

def get_zipper_counts(dict, guide, nts):
    for nt1 in nts:
        for nt2 in nts:
            regex = '(' + nt1 + '(A|C|T|G)' + nt2 + ')'
            nts_counts = len(re.findall(regex, guide))
            nts_frac = nts_counts/(len(guide) - 2)
            dict[nt1 + 'N' + nt2] = nts_frac
    return dict

def get_rep_counts(dict, context, nts, length):
    for nt in nts:
        k_mer = nt*length
        nts_counts = context.count(k_mer)
        nts_frac = nts_counts/(len(context) - 2)
        dict[nt + '*' + str(length)] = nts_frac
    return dict

def get_double_zipper(dict, context_sequence, nts, context_order):
    for i in range(len(context_order) - 3):
        curr_nts = context_sequence[i] + context_sequence[i+3]
        for nt1 in nts:
            for nt2 in nts:
                match_nts = nt1+nt2
                key = context_order[i] + nt1 + 'NN' + nt2
                if curr_nts == match_nts:
                    dict[key] = 1
                else:
                    dict[key] = 0
    return dict


def featurize_guides(kmers, features = None,
                     guide_start = 9, guide_length = 20):
    """Take guides and encodes for modeling
    :param kmers: vector of context sequences
    :param features: boolean dictionary of which feature types to inlcude:
        ['Pos. Ind. 1mer', 'Pos. Ind. 2mer', 'Pos. Ind. 3mer',
        'Pos. Ind. Zipper','Pos. Dep. 1mer', 'Pos. Dep. 2mer',
        'Pos. Dep. 3mer', 'Pos. Dep. Zipper', 'Pos. Ind. Rep.', 'GC content',
        'Tm', 'Physio', 'Double Zipper']
    :param pam_start: int
    :param pam_end: int
    :param guide_start: int
    :param guide_end: int
    :return: featurized matrix
    """
    if features == None:
        # RS2
        features = ['Pos. Ind. 1mer', 'Pos. Ind. 2mer',
                    'Pos. Dep. 1mer', 'Pos. Dep. 2mer',
                    'GC content', 'Tm']
    possible_feats = {'Pos. Ind. 1mer', 'Pos. Ind. 2mer', 'Pos. Ind. 3mer',
                      'Pos. Ind. Zipper','Pos. Dep. 1mer', 'Pos. Dep. 2mer',
                      'Pos. Dep. 3mer', 'Pos. Dep. Zipper', 'Pos. Ind. Rep.',
                      'GC content', 'Tm', 'Physio', 'Double Zipper'}
    if not set(features).issubset(possible_feats):
        diff = features - possible_feats
        assert ValueError(str(diff) + 'Are not currently supported as features')
    current_path = os.path.abspath(os.path.dirname(__file__))
    physio_path = os.path.join(current_path, 'data/features/physiochem.csv.zip')
    physiochemical_data = pd.read_csv(physio_path)
    k = len(kmers[0])
    context_order = get_context_order(k)
    nts = ['A', 'C', 'T', 'G']
    feature_dict_list = []
    for i in range(len(kmers)):
        curr_dict = {}
        context = kmers[i]
        guide_sequence = get_guide_sequence(context, guide_start, guide_length)
        if 'GC content' in features:
            curr_dict = get_frac_g_or_c(curr_dict, guide_sequence)
        if 'Pos. Ind. 1mer' in features:
            curr_dict = get_one_nt_counts(curr_dict, guide_sequence, nts)
        if 'Pos. Ind. 2mer' in features:
            curr_dict = get_two_nt_counts(curr_dict, guide_sequence, nts)
        if 'Pos. Ind. 3mer' in features:
            curr_dict = get_three_nt_counts(curr_dict, guide_sequence, nts)
        if 'Pos. Dep. 1mer' in features:
            curr_dict = get_one_nt_pos(curr_dict, context, nts, context_order)
        if 'Pos. Dep. 2mer' in features:
            curr_dict = get_two_nt_pos(curr_dict, context, nts, context_order)
        if 'Pos. Dep. 3mer' in features:
            curr_dict = get_three_nt_pos(curr_dict, context, nts, context_order)
        if 'Tm' in features:
            curr_dict = get_thermo(curr_dict, guide_sequence, context)
        if 'Physio' in features:
            curr_dict = get_physiochemical(curr_dict, guide_sequence, nts, physiochemical_data)
        if 'Pos. Dep. Zipper' in features:
            curr_dict = get_zipper_pos(curr_dict, context, nts, context_order)
        if 'Pos. Ind. Zipper' in features:
            curr_dict = get_zipper_counts(curr_dict, guide_sequence, nts)
        if 'Pos. Ind. Rep.' in features:
            curr_dict = get_rep_counts(curr_dict, context, nts, 4)
        if 'Double Zipper' in features:
            curr_dict = get_double_zipper(curr_dict, context, nts, context_order)
        feature_dict_list.append(curr_dict)
    feature_matrix = pd.DataFrame(feature_dict_list)
    return feature_matrix
