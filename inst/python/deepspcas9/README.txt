1. System Requirements:
	Ubuntu	16.04
	Python	2.7.12
	Python Packages:
		numpy 1.14.5
		scipy 1.1.0

	Tensorflow and dependencies:
		Tensorflow  1.4.1
		CUDA	    8.0.61
		cuDNN	    5.1.10

2. Installation Guide (required time, <120 minutes):
	Package Installation:
		pip install numpy==1.14.5
		pip install scipy==1.1.0
		pip install tensorflow==1.4.1


3. Demo Instructions (required time, <1 min):

Input1: ./dataset/        # List of Target Sequence(s)
	File format:
	Target number   30 bp target sequence (4 bp + 20 bp protospacer + PAM + 3 bp)
	  1   TAAGAGAGTGGTAATAGAAGTGCCAGGTAT
	  2   CCCTCATGGTGCAGCTAAAGGCCCAGGAGC

Input2: ./DeepCas9_Final/ # Pre-trained Weight Files

Output: RANK_final_DeepCas9_Final.txt
	Predicted activity score for sequence 1 and 2:
	67.5565185546875, 56.930904388427734   

Run script:
	python ./DeepCas9_TestCode.py

Modification for personalized runs:

	<DeepCas9_TestCode.py>
	## System Paths ##
	path                 = './dataset/'
	parameters           = {'0': 'sample.txt'}

	## Run Parameters ##
	TEST_NUM_SET         = [0] # List can be expanded in case of multiple test parameters
	best_model_path_list = ['./DeepCas9_Final/']

sample.txt can be replaced or modified to include target sequence of interest

 


