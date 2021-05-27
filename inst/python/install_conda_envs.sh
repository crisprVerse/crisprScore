# Setting up conda env for azimuth:
conda create -n azimuth python=2.7
conda activate azimuth
pip install ./azimuth

# Setting up conda env for deepcpf1:
conda create -n deepcpf1 python=2.7
conda activate deepcpf1
pip install ./deepcpf1

# Setting up conda env for lindel:
conda create -n lindel python=3.6
conda activate lindel
pip install ./lindel

# Setting up conda env for deepcas9
conda create -n deepcas9 python=3.6
conda activate deepcas9
pip install ./deepcas9
#Notes:
#1. Tensorflow has to be used as the backend for KERAS
#2. ViennaRNA has to be installed, and the path added to .getPaths()
###########################################################################

# Setting up conda env for deepcas9
conda create -n enpamgb python=3.6
conda activate enpamgb2
pip install ./enpamgb
###########################################################################