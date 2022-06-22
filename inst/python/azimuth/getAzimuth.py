#sys.argv[1] should be the path of the file with input sequences
#sys.argv[2] should be the path of the file where to save sequence
import sys
import azimuth.model_comparison
import numpy as np 

def getAzimuth(sequences):
	predictions = azimuth.model_comparison.predict(sequences, None, None)
	return predictions

sequences = np.loadtxt(sys.argv[1], dtype="U34", ndmin=1)
scores = getAzimuth(sequences)
np.savetxt(sys.argv[2], scores)
