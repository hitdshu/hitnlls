#pragma once

#include <algorithm>
#include "vertex/vertex_imp.h"

namespace hitcadmm {

template <int ndim>
class Soc : public VertexImp<ndim, Eigen::Matrix<double, ndim, 1>> {
public:
    typedef std::shared_ptr<Soc> Ptr;

    explicit Soc() : VertexImp<ndim, Eigen::Matrix<double, ndim, 1>>() {}
    explicit Soc(const Eigen::Matrix<double, ndim, 1> &val) : VertexImp<ndim, Eigen::Matrix<double, ndim, 1>>() { this->SetValue(val); }

    virtual void Project() override {
        Eigen::Matrix<double, ndim, 1> val = this->GetValue();
        Eigen::Matrix<double, ndim - 1, 1> v = val.block(0, 0, ndim - 1, 1);
        double t = val[ndim - 1];
        double vn = v.norm();
        if (t <= -vn) {
            val.setZero();
        } else if (t >= vn) {
        } else {
            double a = (vn + t) / 2;
            val.block(0, 0, ndim - 1, 1) *= a / vn;
            val[ndim - 1] = a;
        }
        this->SetValue(val);
    }

    virtual void SetVector(const Eigen::VectorXd &v) override { if (this->CheckDim(v)) this->SetValue(v); }
    virtual Eigen::VectorXd GetVector() const override { return this->GetValue(); }
};

} // namespace hitcadmm