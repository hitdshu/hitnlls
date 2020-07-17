#pragma once

#include <algorithm>

#include "variable.h"

namespace nlls {

template <int N>
class Nnc : public VariableImp<N, Matrix<float, N, 1>> {
public:
    using BaseType = VariableImp<N, Matrix<float, N, 1>>;
    explicit Nnc() : BaseType() {}
    explicit Nnc(const Matrix<float, N, 1> &val) : BaseType() { BaseType::SetValue(val); }

    virtual void Project() override {
        using std::max;
        Matrix<float, N, 1> val = BaseType::GetValue();
        for (int idx = 0; idx < N; ++idx) {
            val[idx] = max<float>(0, val[idx]);
        }
        BaseType::SetValue(val);
    }
    virtual void SetVector(const VectorXf &v) override { if (BaseType::CheckDim(v)) { BaseType::SetValue(v); } }
    virtual VectorXf GetVector() const override { return BaseType::GetValue(); }
};

class NncX : public VariableImpX<VectorXf> {
public:
    using BaseType = VariableImpX<VectorXf>;
    explicit NncX(const VectorXf &val) : BaseType(val.Size()) { BaseType::SetValue(val); }

    virtual void Project() override;
    virtual void SetVector(const VectorXf &v) override { if (BaseType::CheckDim(v)) { BaseType::SetValue(v); } }
    virtual VectorXf GetVector() const override { return BaseType::GetValue(); }
};

} // namespace nlls