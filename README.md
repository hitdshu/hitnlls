# hitnllls
A nonlinear least square(NLLS) solver. Fomulate the NLLS as graph optimization. For the name, hit here refers to Harbin Institute of Technology to honor the six wonderful years spent there and nlls is short for nonlinear least square.

## Abstract
This solver uses Gauss-Newton type algorithms for NLLS problems. We define node to be variables waiting for optimization and factor to be observation over which jacobian and error are evaluated. Without any dependence on other libraries, we choose to implement dense matrix and sparse matrix operations used in later computation. As for the code integration, we use catch to support unit tests. This repo will update continuously (and slowly...). 

## Implementation
- matrix: folder for dense and sparse matrix operations. We mainly implement Cholesky/LUP/Inverse operation for dense matrix and sparse cholesky for sparse matrix.
- common: folder for type register currently.
- factor: folder for error and jacobian evaluation of several factors. We mainly implement some linear factors and SE2/SE3 related factors. Included also are robust kernels and several types of cameras.
- graph: folder for optimizable graph which organizes nodes and factors. Currently only some naive graphs are available.
- node: folder for optimizable variables, mainly SE2/SE3 related types.
- solver: folder for optimization algorithms. We implement Gauss-Newton/LM/Preconditioned Conjugate Gradient methods.
- utils: folder for utilities. Now only timer class is available.
- test: folder for test cases. Have fun with them.

# hitcadmm
In directory cadmm. We recently(2020/2/14) developed a cone ADMM algorithm for convex cone optimization problem(LP/SOCP/SDP). Check files there if you feel interested. 
