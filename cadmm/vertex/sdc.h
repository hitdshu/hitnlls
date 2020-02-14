#pragma once

#include <algorithm>
#include "vertex/vertex_imp.h"

namespace hitcadmm {

template <int ndim>
Eigen::Matrix<double, ndim, ndim> Vec2Mat(const Eigen::Matrix<double, ndim * ndim, 1> &v) {
    Eigen::Matrix<double, ndim, ndim> mat;
    for (int cidx = 0; cidx < ndim; ++cidx) {
        for (int ridx = 0; ridx < ndim; ++ridx) {
            mat(ridx, cidx) = v[cidx * ndim + ridx];
        }
    }
    return mat;
}

template <int ndim>
Eigen::Matrix<double, ndim * ndim, 1> Mat2Vec(const Eigen::Matrix<double, ndim, ndim> &mat) {
    Eigen::Matrix<double, ndim * ndim, 1> v;
    for (int cidx = 0; cidx < ndim; ++cidx) {
        for (int ridx = 0; ridx < ndim; ++ridx) {
            v[cidx * ndim + ridx] = mat(ridx, cidx);
        }
    }
    return v;
}

template <int ndim>
class Sdc : public VertexImp<ndim * ndim, Eigen::Matrix<double, ndim, ndim>> {
public:
    typedef std::shared_ptr<Sdc> Ptr;

    explicit Sdc() : VertexImp<ndim * ndim, Eigen::Matrix<double, ndim, ndim>>() {}
    explicit Sdc(const Eigen::Matrix<double, ndim, ndim> &val) : VertexImp<ndim * ndim, Eigen::Matrix<double, ndim, ndim>>() { this->SetValue(val); }

    virtual void Project() override {
        Eigen::Matrix<double, ndim, ndim> val = this->GetValue();
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, ndim, ndim>> solver(val);
        Eigen::Matrix<double, ndim, 1> eigen_vals = solver.eigenvalues();
        Eigen::Matrix<double, ndim, ndim> eigen_vecs = solver.eigenvectors();
        for (int idx = 0; idx < ndim; ++idx) {
            eigen_vals[idx] = std::max<double>(0, eigen_vals[idx]);
        }
        val = eigen_vecs * eigen_vals.asDiagonal() * eigen_vecs.transpose();
        this->SetValue(val);
    }

    virtual void SetVector(const Eigen::VectorXd &v) override {
        if (this->CheckDim(v)) {
            this->SetValue(Vec2Mat<ndim>(v));
        }
    }

    virtual Eigen::VectorXd GetVector() const override {
        return Mat2Vec<ndim>(this->GetValue()); 
    }
};

} // namespace hitcadmm