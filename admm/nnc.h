#pragma once

#include <algorithm>

#include "admm/variable.h"

namespace hitcadmm {

template <int ndim>
class Nnc : public VariableImp<ndim, hitnlls::matrix::Matrix<double, ndim, 1>> {
public:
    typedef std::shared_ptr<Nnc> Ptr;

    explicit Nnc() : VariableImp<ndim, hitnlls::matrix::Matrix<double, ndim, 1>>() {}
    explicit Nnc(const hitnlls::matrix::Matrix<double, ndim, 1> &val) : VariableImp<ndim, hitnlls::matrix::Matrix<double, ndim, 1>>() { this->SetValue(val); }

    virtual void Project() override {
        hitnlls::matrix::Matrix<double, ndim, 1> val = this->GetValue();
        for (int idx = 0; idx < ndim; ++idx) {
            val[idx] = std::max<double>(0, val[idx]);
        }
        this->SetValue(val);
    }

    virtual void SetVector(const hitnlls::matrix::VectorXd &v) override { if (this->CheckDim(v)) this->SetValue(v); }
    virtual hitnlls::matrix::VectorXd GetVector() const override { return this->GetValue(); }
};

} // namespace hitcadmm