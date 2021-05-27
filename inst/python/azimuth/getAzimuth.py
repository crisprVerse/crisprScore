import sys
import azimuth.model_comparison
import numpy as np 

def getAzimuth(sequences):
	predictions = azimuth.model_comparison.predict(sequences, None, None)
	return predictions