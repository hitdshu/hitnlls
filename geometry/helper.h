#pragma once

#include "matrix/dense.h"

namespace hitnlls {
namespace geometry {

template <typename T>
matrix::Matrix<T, 3, 3> Skew(const matrix::Matrix<T, 3, 1> &v) {
    matrix::Matrix<T, 3, 3> result;
    result << 0, -v[2], v[1], 
        v[2], 0, -v[0], 
        -v[1], v[0], 0;
    return result;
}

} // namespace geometry
} // namespace hitnlls