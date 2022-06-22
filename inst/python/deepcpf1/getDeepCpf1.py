#sys.argv[1] should be the path of the file with input sequences
#sys.argv[2] should be the path of the file where to save sequence
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

sequences = np.loadtxt(sys.argv[1], dtype="U34", ndmin=1)
scores = getDeepCpf1(sequences)
np.savetxt(sys.argv[2], scores)