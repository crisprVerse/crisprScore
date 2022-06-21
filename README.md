# crisprScore

## Overview

crisprScore Provides R wrappers of several on-target and off-target scoring methods for CRISPR guide RNAs (gRNAs).
    The following nucleases are supported: SpCas9, AsCas12a, enAsCas12a, and RfxCas13d (CasRx).
    The available on-target cutting efficiency scoring methods are RuleSet1, Azimuth, DeepHF, DeepCpf1, 
    enPAM+GB, and CRISPRscan. Both the CFD and MIT scoring methods are available for off-target
    specificity prediction. The package also provides a Lindel-derived score to predict the probability
    of a gRNA to produce indels inducing a frameshift for the Cas9 nuclease.
    Note that DeepHF, DeepCpf1 and enPAM+GB are not available on Windows machines. 

Our work is described in a recent bioRxiv preprint: ["A comprehensive Bioconductor ecosystem for the design of CRISPR guide RNAs across nucleases and technologies"](https://www.biorxiv.org/content/10.1101/2022.04.21.488824v2)


## Software requirements

### OS Requirements

This package is supported for macOS, Linux and Windows machines.
Some functionalities are not supported for Windows machines.

Packages were developed and tested on R version 4.2.

### R Dependencies 

- crisprScoreData: https://github.com/Jfortin1/crisprScoreData


## Installation

`crisprScore` and its dependencies can be installed by typing the following commands inside of an R session:

```r
install.packages("devtools")
library(devtools)
install_github("Jfortin1/crisprScoreData")
install_github("Jfortin1/crisprScore")
```

## Demo 

A reproducible and comprehensive workflow can be found here:
https://github.com/Jfortin1/crisprScore/blob/master/vignettes/crisprScore.Rmd


## License 

This project is covered under the MIT License.

