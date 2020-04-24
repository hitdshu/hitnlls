#pragma once

#include <algorithm>

#include "admm/variable.h"

namespace hitcadmm {

template <int ndim>
hitnlls::matrix::Matrix<double, ndim, ndim> Vec2Mat(const hitnlls::matrix::Matrix<double, ndim * ndim, 1> &v) {
    hitnlls::matrix::Matrix<double, ndim, ndim> mat;
    for (int cidx = 0; cidx < ndim; ++cidx) {
        for (int ridx = 0; ridx < ndim; ++ridx) {
            mat(ridx, cidx) = v[cidx * ndim + ridx];
        }
    }
    return mat;
}

template <int ndim>
hitnlls::matrix::Matrix<double, ndim * ndim, 1> Mat2Vec(const hitnlls::matrix::Matrix<double, ndim, ndim> &mat) {
    hitnlls::matrix::Matrix<double, ndim * ndim, 1> v;
    for (int cidx = 0; cidx < ndim; ++cidx) {
        for (int ridx = 0; ridx < ndim; ++ridx) {
            v[cidx * ndim + ridx] = mat(ridx, cidx);
        }
    }
    return v;
}

template <int ndim>
class Sdc : public VariableImp<ndim * ndim, hitnlls::matrix::Matrix<double, ndim, ndim>> {
public:
    typedef std::shared_ptr<Sdc> Ptr;

    explicit Sdc() : VariableImp<ndim * ndim, hitnlls::matrix::Matrix<double, ndim, ndim>>() {}
    explicit Sdc(const hitnlls::matrix::Matrix<double, ndim, ndim> &val) : VariableImp<ndim * ndim, hitnlls::matrix::Matrix<double, ndim, ndim>>() { this->SetValue(val); }

    virtual void Project() override {
        hitnlls::matrix::Matrix<double, ndim, ndim> val = this->GetValue();
        hitnlls::matrix::EVD<hitnlls::matrix::Matrix<double, ndim, ndim>> solver(val);
        hitnlls::matrix::Matrix<double, ndim, 1> eigen_vals = solver.V();
        hitnlls::matrix::Matrix<double, ndim, ndim> eigen_vecs = solver.U();
        for (int idx = 0; idx < ndim; ++idx) {
            eigen_vals[idx] = std::max<double>(0, eigen_vals[idx]);
        }
        hitnlls::matrix::Matrix<double, ndim, ndim> eigen_vals_mat;
        eigen_vals_mat.SetIdentity();
        for (int idx = 0; idx < ndim; ++idx) {
            eigen_vals_mat(idx, idx) = eigen_vals[idx];
        }
        val = eigen_vecs * eigen_vals_mat * eigen_vecs.Transpose();
        this->SetValue(val);
    }

    virtual void SetVector(const hitnlls::matrix::VectorXd &v) override {
        if (this->CheckDim(v)) {
            this->SetValue(Vec2Mat<ndim>(v));
        }
    }

    virtual hitnlls::matrix::VectorXd GetVector() const override {
        return Mat2Vec<ndim>(this->GetValue()); 
    }
};

} // namespace hitcadmm