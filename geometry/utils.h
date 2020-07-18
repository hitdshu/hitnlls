#pragma once

#include "../matrix/dense.h"

namespace nlls {

template <typename T>
Matrix<T, 2, 2> Skew2(const Matrix<T, 1, 1> &v) {
    Matrix<T, 2, 2> result;
    result << T(0), -v[0], 
        v[0], T(0);
    return result;
}

template <typename T>
Matrix<T, 2, 2> Skew2(const T &v) {
    Matrix<T, 2, 2> result;
    result << T(0), -v, 
        v, T(0);
    return result;
}

template <typename T>
Matrix<T, 3, 3> Skew3(const Matrix<T, 3, 1> &v) {
    Matrix<T, 3, 3> result;
    result << T(0), -v[2], v[1], 
        v[2], T(0), -v[0], 
        -v[1], v[0], T(0);
    return result;
}

} // namespace nlls