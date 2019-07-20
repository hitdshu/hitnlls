# hitnllls
A nonlinear least square(NLLS) solver. Fomulate the NLLS as graph optimization. For the name, hit here refers to Harbin Institute of Technology to honor the six wonderful years spent there and nlls is short for nonlinear least square.

## Abstract
This solver uses Gauss-Newton type algorithms for NLLS problems. We define node to be variables waiting for optimization and factor to be observation over which jacobian and error are evaluated. Without any dependence on other libraries, we choose to implement dense matrix and sparse matrix operations used in later computation. As for the code integration, we use catch to support unit tests. This repo will update continuously (and slowly...). 