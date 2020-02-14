#pragma once

#include <Eigen/Dense>
#include "vertex/vertex_imp.h"

namespace hitcadmm {

template <int ndim>
class Var : public VertexImp<ndim, Eigen::Matrix<double, ndim, 1> > {
public:
    typedef std::shared_ptr<Var> Ptr;

    explicit Var() : VertexImp<ndim, Eigen::Matrix<double, ndim, 1>>() {}
    explicit Var(const Eigen::Matrix<double, ndim, 1> &val) : VertexImp<ndim, Eigen::Matrix<double, ndim, 1>>() { this->SetValue(val); }

    virtual void Project() override { return; }
    virtual void SetVector(const Eigen::VectorXd &v) override { if (this->CheckDim(v)) this->SetValue(v); }
    virtual Eigen::VectorXd GetVector() const override { return this->GetValue(); }
};

} // namespace hitcadmm