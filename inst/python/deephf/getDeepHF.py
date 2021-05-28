import sys
import os
from os.path import dirname, abspath, join
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
from deephf.training_util import *
from deephf.prediction_util import *

#enzyme choices: 'wt_u6', 'wt_t7', 'esp' or 'hf'
def getDeepHF(sequences, enzyme='wt_u6', model_file=''):
	results = get_scores(sequences, enzyme, model_file)
	return results

#sequences = np.array(["ACGTGTGACTACCGGCGGCGCGG",
#	"GGAAGTCTGGAGTCTCCAGGTGG"])
#results = getDeepWt(sequences)