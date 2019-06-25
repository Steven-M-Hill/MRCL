# Manifold Regularized Causal Learning (MRCL)
R code for Manifold Regularized Causal Learning (MRCL) and scripts to run the analysis presented in:
#### Hill, Oates, Blythe &amp; Mukherjee (2019). Causal Learning via Manifold Regularization. arXiv preprint [arXiv:1612.05678](https://arxiv.org/abs/1612.05678).

### MRCL source code
The file [`mrcl.R`](./code/mrcl.R) in the [code](code) directory contains the source code for running MRCL.  
Functions are documented within this file.  
This code was co-authored with [Umberto Noè](https://github.com/unoe). 

##### Required installs:
```
install.packages('kernlab')
```

### Manuscript analysis scripts

R scripts and functions are provided in the [code](code) directory to reproduce the results in Section 3 of the manuscript.  

To run the analysis scripts, your working directory in R must contain the [code](code) directory and the [data](data) directory.  
Details of the origin of the three data files in the [data](data) directory are provided below.

##### Required installs:
```
install.packages('pROC', 'glmnet', 'caret', 'pcalg', 'dplyr', 
                 'kernlab', 'doParallel', 'doRNG', 'ggplot2', 'RColorBrewer')
```

##### Dataset D1: Yeast Gene Expression (Section 3.2 of manuscript)

<!--note that code is provided that obtains two of the files ([yeastData.RData](data/yeastData.RData) and [cellLineData.RData](data/cellLineData.RData)) from their source and performs post-processing.
'R.matlab',--> 
