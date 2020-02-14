#pragma once

#include <algorithm>
#include "vertex/vertex_imp.h"

namespace hitcadmm {

template <int ndim>
class Nnc : public VertexImp<ndim, Eigen::Matrix<double, ndim, 1>> {
public:
    typedef std::shared_ptr<Nnc> Ptr;

    explicit Nnc() : VertexImp<ndim, Eigen::Matrix<double, ndim, 1>>() {}
    explicit Nnc(const Eigen::Matrix<double, ndim, 1> &val) : VertexImp<ndim, Eigen::Matrix<double, ndim, 1>>() { this->SetValue(val); }

    virtual void Project() override {
        Eigen::Matrix<double, ndim, 1> val = this->GetValue();
        for (int idx = 0; idx < ndim; ++idx) {
            val[idx] = std::max<double>(0, val[idx]);
        }
        this->SetValue(val);
    }

    virtual void SetVector(const Eigen::VectorXd &v) override { if (this->CheckDim(v)) this->SetValue(v); }
    virtual Eigen::VectorXd GetVector() const override { return this->GetValue(); }
};

} // namespace hitcadmm