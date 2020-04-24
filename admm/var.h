#pragma once

#include "admm/variable.h"

namespace hitcadmm {

template <int ndim>
class Var : public VariableImp<ndim, hitnlls::matrix::Matrix<double, ndim, 1> > {
public:
    typedef std::shared_ptr<Var> Ptr;

    explicit Var() : VariableImp<ndim, hitnlls::matrix::Matrix<double, ndim, 1>>() {}
    explicit Var(const hitnlls::matrix::Matrix<double, ndim, 1> &val) : VariableImp<ndim, hitnlls::matrix::Matrix<double, ndim, 1>>() { this->SetValue(val); }

    virtual void Project() override { return; }
    virtual void SetVector(const hitnlls::matrix::VectorXd &v) override { if (this->CheckDim(v)) this->SetValue(v); }
    virtual hitnlls::matrix::VectorXd GetVector() const override { return this->GetValue(); }
};

} // namespace hitcadmm