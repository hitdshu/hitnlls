#pragma once

#include <algorithm>

#include "variable.h"

namespace nlls {

template <int N>
Matrix<float, N, N> Vec2Mat(const Matrix<float, N * N, 1> &v) {
    Matrix<float, N, N> mat;
    for (int cidx = 0; cidx < N; ++cidx) {
        for (int ridx = 0; ridx < N; ++ridx) {
            mat(ridx, cidx) = v[ridx * N + cidx];
        }
    }
    return mat;
}

MatrixXf VecX2MatX(const VectorXf &v) {
    using std::round;
    int n = round(sqrt(v.Size()));
    MatrixXf mat(n, n);
    for (int cidx = 0; cidx < n; ++cidx) {
        for (int ridx = 0; ridx < n; ++ridx) {
            mat(ridx, cidx) = v[ridx * n + cidx];
        }
    }
    return mat;
}

template <int N>
Matrix<float, N * N, 1> Mat2Vec(const Matrix<float, N, N> &mat) {
    Matrix<float, N * N, 1> v;
    for (int ridx = 0; ridx < N; ++ridx) {
        for (int cidx = 0; cidx < N; ++cidx) {
            v[ridx * N + cidx] = mat(ridx, cidx);
        }
    }
    return v;
}

VectorXf MatX2VecX(const MatrixXf &mat) {
    int n = mat.Rows();
    VectorXf v(mat.Size());
    for (int ridx = 0; ridx < n; ++ridx) {
        for (int cidx = 0; cidx < n; ++cidx) {
            v[ridx * n + cidx] = mat(ridx, cidx);
        }
    }
    return v;
}

template <int N>
class Sdc : public VariableImp<N * N, Matrix<float, N, N>> {
public:
    using BaseType = VariableImp<N * N, Matrix<float, N, N>>;
    explicit Sdc() : BaseType() {}
    explicit Sdc(const Matrix<float, N, N> &val) : BaseType() { BaseType::SetValue(val); }

    virtual void Project() override {
        using std::max;
        Matrix<float, N, N> val = BaseType::GetValue();
        EVD<Matrix<float, N, N>> solver(val);
        Matrix<float, N, 1> eigen_vals = solver.V();
        Matrix<float, N, N> eigen_vecs = solver.U();
        for (int idx = 0; idx < N; ++idx) {
            eigen_vals[idx] = max<float>(0.0, eigen_vals[idx]);
        }
        Matrix<float, N, N> eigen_vals_mat;
        eigen_vals_mat.SetIdentity();
        for (int idx = 0; idx < N; ++idx) {
            eigen_vals_mat(idx, idx) = eigen_vals[idx];
        }
        val = eigen_vecs * eigen_vals_mat * eigen_vecs.Transpose();
        BaseType::SetValue(val);
    }
    virtual void SetVector(const VectorXf &v) override { if (BaseType::CheckDim(v)) { BaseType::SetValue(Vec2Mat<N>(v)); } }
    virtual VectorXf GetVector() const override { return Mat2Vec<N>(BaseType::GetValue()); }
};

class SdcX : public VariableImpX<MatrixXf> {
public:
    using BaseType = VariableImpX<MatrixXf>;
    explicit SdcX(const MatrixXf &val) : BaseType(val.Size()) { BaseType::SetValue(val); }

    virtual void Project() override;
    virtual void SetVector(const VectorXf &v) override { if (BaseType::CheckDim(v)) { BaseType::SetValue(VecX2MatX(v)); } }
    virtual VectorXf GetVector() const override { return MatX2VecX(BaseType::GetValue()); }
};

} // namespace nlls