#pragma once

#include <algorithm>

#include "variable.h"

namespace nlls {

template <int N>
class Soc : public VariableImp<N, Matrix<float, N, 1>> {
public:
    using BaseType = VariableImp<N, Matrix<float, N, 1>>;
    explicit Soc() : BaseType() {}
    explicit Soc(const Matrix<float, N, 1> &val) : BaseType() { BaseType::SetValue(val); }

    virtual void Project() override {
        Matrix<float, N, 1> val = BaseType::GetValue();
        Matrix<float, N - 1, 1> v = val.Block(0, 0, N - 1, 1);
        float t = val[N - 1];
        float vn = v.Norm();
        if (t <= -vn) {
            val.SetZero();
        } else if (t >= vn) {
        } else {
            float a = (vn + t) / 2;
            val.Block(0, 0, N - 1, 1) *= a / vn;
            val[N - 1] = a;
        }
        BaseType::SetValue(val);
    }
    virtual void SetVector(const VectorXf &v) override { if (BaseType::CheckDim(v)) { BaseType::SetValue(v); } }
    virtual VectorXf GetVector() const override { return BaseType::GetValue(); }
};

class SocX : public VariableImpX<VectorXf> {
public:
    using BaseType = VariableImpX<VectorXf>;
    explicit SocX(const VectorXf &val) : BaseType(val.Size()) { BaseType::SetValue(val); }

    virtual void Project() override;
    virtual void SetVector(const VectorXf &v) override { if (BaseType::CheckDim(v)) { BaseType::SetValue(v); } }
    virtual VectorXf GetVector() const override { return BaseType::GetValue(); }
};

} // namespace nlls