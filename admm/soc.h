#pragma once

#include <algorithm>

#include "admm/variable.h"

namespace hitcadmm {

template <int ndim>
class Soc : public VariableImp<ndim, hitnlls::matrix::Matrix<double, ndim, 1>> {
public:
    typedef std::shared_ptr<Soc> Ptr;

    explicit Soc() : VariableImp<ndim, hitnlls::matrix::Matrix<double, ndim, 1>>() {}
    explicit Soc(const hitnlls::matrix::Matrix<double, ndim, 1> &val) : VariableImp<ndim, hitnlls::matrix::Matrix<double, ndim, 1>>() { this->SetValue(val); }

    virtual void Project() override {
        hitnlls::matrix::Matrix<double, ndim, 1> val = this->GetValue();
        hitnlls::matrix::Matrix<double, ndim - 1, 1> v = val.Block(0, 0, ndim - 1, 1);
        double t = val[ndim - 1];
        double vn = v.Norm();
        if (t <= -vn) {
            val.SetZero();
        } else if (t >= vn) {
        } else {
            double a = (vn + t) / 2;
            val.Block(0, 0, ndim - 1, 1) *= a / vn;
            val[ndim - 1] = a;
        }
        this->SetValue(val);
    }

    virtual void SetVector(const hitnlls::matrix::VectorXd &v) override { if (this->CheckDim(v)) this->SetValue(v); }
    virtual hitnlls::matrix::VectorXd GetVector() const override { return this->GetValue(); }
};

} // namespace hitcadmm