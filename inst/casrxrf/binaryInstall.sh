# The ViennaRNA and RNAhybrid packages must be installed
# so that the following binaries are available to crisprScore:
# - viennarna/bin/RNAfold
# - viennarna/bin/RNAplfold
# - RNAhybrid/bin/RNAhybrid

# For RNAhybrid, I downloaded the source code from
# https://bibiserv.cebitec.uni-bielefeld.de/rnahybrid/
#
# I unzipped the folder, entered the folder, and compiled the code from source using 
# ./configure --prefix=/gstore/data/omni/crispr/programs/RNAhybrid-2.1.2
# make
# make install


# For ViennaRNA, I downlaoded the source code from their website
# https://www.tbi.univie.ac.at/RNA/

# I unzipped the folder, entered the folder, and compiled the code from source using 
# ./configure --prefix=/gstore/data/omni/crispr/programs/ViennaRNA-2.4.9 --without-gsl
# make
# make install
