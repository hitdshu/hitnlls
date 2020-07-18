# hitnlls
A nonlinear least square(NLLS) solver. Fomulate the NLLS as graph optimization. For the name, HIT here refers to Harbin Institute of Technology to honor the six wonderful years spent there and NLLS is short for nonlinear least square. 

## Abstract
This solver uses Gauss-Newton type algorithms for NLLS problems. We define node to be variables waiting for optimization and factor to be observation over which jacobians and errors are evaluated. Without any dependence on other libraries, we choose to implement dense matrix and sparse matrix operations used in later computation. We've updated this new version of HITNLLS, which extensively uses template metaprogramming and CRTP to make the code more structured, more efficient (and more difficult to read). If you want to see old version of HITNLLS, check the other branch of the repo. Currently I'm adding unit tests and examples.

## Usage

### install
```terminal
$ git clone https://github.com/hitdshu/hitnlls.git
$ cd hitnlls && mkdir build && cd build
$ cmake ..
$ make
$ sudo make install
```

### usage in cmake
```cmake
find_package(hitnlls REQUIRED)
include_directories(
${hitnlls_INCLUDE_DIRS}
)
target_link_libraries(Test 
${hitnlls_LIBRARIES})
```

## Implementation

### admm
We developed a cone ADMM algorithm for convex cone optimization problem(LP/SOCP/SDP). Check files there if you feel interested.

### camera
Several camera models are defined here. Our implementation of cameras are integrated naturally into our optimization framework.

### common
This folder defines the register and timer, which are needed by other modules. Included also is a light weight logger.

### geometry
Our customized implementation of several lie groups(SO2/SO3/SE2/SE3) can be utilized for representation and optimization of poses in 2d and 3d space. The jet class here is for auto differentiation just as ceres. The lie groups defined here can be combined with jet to make auto differentiation in mainfolds much simpler.

### matrix
We use CRTP and expression templates to enable lazy evaluation and avoid the generation of temporaray variabls. Our api of dense matrix operations resembles those of Eigen. We mainly implement cholesky/qr/inverse/evd/svd/lup operations for dense matrix. 

### ils
This folder contains optimizable problem which organizes nodes and factors and optimization algorithms. We implement Gauss-Newton/LM/Doglog strategies and DENSE_QR/DENSE_CHOLESKY/SPARSE_CHOLESKY/DENSE_SCHUR algorithms in this folder.

### type
Factors and vertices are defined here. We use variable tuples to avoid dynamic_cast in the evaluation of factors. The framework can support hand diff/auto diff/numeric diff now. Inherits Factor/FactorAutodiff/FactorNumeDiff according to your own purpose.
