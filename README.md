# hitnlls
A nonlinear least square(NLLS) solver. Fomulate the NLLS as graph optimization. For the name, HIT here refers to Harbin Institute of Technology to honor the six wonderful years spent there and NLLS is short for nonlinear least square. If you want to see our old version of HITNLLS, check branch v1.

## Abstract
This solver uses Gauss-Newton type algorithms for NLLS problems. We define node to be variables waiting for optimization and factor to be observation over which jacobians and errors are evaluated. Without any dependence on other libraries, we choose to implement dense matrix and sparse matrix operations used in later computation. We've updated this new version of HITNLLS, which extensively uses template metaprogramming and CRTP to make the code more structured, more efficient (and more difficult to read). If you want to see old version of HITNLLS, check the other branch of the repo.

## Implementation

### admm
We developed a cone ADMM algorithm for convex cone optimization problem(LP/SOCP/SDP). Check files there if you feel interested.

### common
This folder defines the register and timer, which are needed by other modules.

### geometry
Our customized implementation of several lie groups(SO2/SO3/SE2/SE3) can be utilized for representation and optimization of poses in 2d and 3d space. The jet class here is for auto differentiation just as ceres. The lie groups defined here can be combined with jet to make auto differentiation in mainfolds much simpler.

### matrix
We use CRTP and expression templates to enable lazy evaluation and avoid the generation of temporaray variabls. Our api of dense matrix operations resembles those of Eigen. We mainly implement cholesky/qr/inverse/evd/svd/lup operations for dense matrix. 

### nlls
This folder contains optimizable graph which organizes nodes and factors and optimization algorithms. Currently only some simple graphs are available. We implement Gauss-Newton/LM algorithms in this folder. This folder also has our implementation of sparse block matrix and several operations on it.

### type
Factors, vertices and cameras are defined here. We use variable tuples to avoid dynamic_cast in the evaluation of factors. The framework can support hand diff/auto diff/numeric diff now. Inherits Factor/FactorAutodiff/FactorNumeDiff according to your own purpose.
