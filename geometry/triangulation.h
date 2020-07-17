#pragma once

#include "../matrix/dense.h"

#include "utils.h"

namespace nlls {

template <typename T>
Matrix<T, 3, 1> TriangulateDLT(const Matrix<T, 4, 4> &tc1w, const Matrix<T, 3, 1> &ray1, const Matrix<T, 4, 4> &tc2w, const Matrix<T, 3, 1> &ray2) {
    Matrix<T, 3, 1> ray1n = ray1 / ray1[2];
    Matrix<T, 3, 1> ray2n = ray2 / ray2[2];
    Matrix<T, 4, 4> A;
    A.Row(0) = tc1w.Row(0) - ray1n[0] * tc1w.Row(2);
    A.Row(1) = tc1w.Row(1) - ray1n[1] * tc1w.Row(2);
    A.Row(2) = tc2w.Row(0) - ray2n[0] * tc2w.Row(2);
    A.Row(3) = tc2w.Row(1) - ray2n[1] * tc2w.Row(2);
    SVD<Matrix<T, 4, 4>> svd(A);
    Matrix<T, 4, 1> v_last = svd.V().Col(3);
    Matrix<T, 3, 1> result;
    result = v_last.Block(0,0,3,1) / v_last[3];
    return result;
}

template <typename T>
Matrix<T, 3, 1> TriangulateMLD(const Matrix<T, 4, 4> &tc1w, const Matrix<T, 3, 1> &ray1, const Matrix<T, 4, 4> &tc2w, const Matrix<T, 3, 1> &ray2) {
    Matrix<T, 3, 1> ray1n = ray1.Normalized();
    Matrix<T, 3, 1> ray2n = ray2.Normalized();
    Matrix<T, 4, 4> twc1 = tc1w.Inverse();
    Matrix<T, 4, 4> twc2 = tc2w.Inverse();
    ray1n = twc1.Block(0, 0, 3, 3) * ray1n;
    ray2n = twc2.Block(0, 0, 3, 3) * ray2n;
    Matrix<T, 3, 3> ray1skew = Skew3<T>(ray1n);
    Matrix<T, 3, 3> ray2skew = Skew3<T>(ray2n);
    Matrix<T, 6, 3> A;
    Matrix<T, 6, 1> b;
    A.Block(0, 0, 3, 3) = ray1skew;
    A.Block(3, 0, 3, 3) = ray2skew;
    b.Block(0, 0, 3, 1) = ray1skew * twc1.Block(0, 3, 3, 1);
    b.Block(3, 0, 3, 1) = ray2skew * twc2.Block(0, 3, 3, 1);
    QR<Matrix<T, 6, 3>> qr(A);
    return qr.Solve(b);
}

} // namespace nlls