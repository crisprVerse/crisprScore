import os
import numpy as np
import warnings
with warnings.catch_warnings():	
	warnings.simplefilter('ignore')
	from deepcpf1.deepcpf1 import *
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

def getDeepCpf1(sequences):	
	results = deepcpf1(sequences)
	return results