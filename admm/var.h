#pragma once

#include "variable.h"

namespace nlls {

template <int N>
class Var : public VariableImp<N, Matrix<float, N, 1>> {
public:
    using BaseType = VariableImp<N, Matrix<float, N, 1>>;
    explicit Var() : BaseType() {}
    explicit Var(const Matrix<float, N, 1> &val) : BaseType() { BaseType::SetValue(val); }

    virtual void Project() override { return; }
    virtual void SetVector(const VectorXf &v) override { if (BaseType::CheckDim(v)) { BaseType::SetValue(v); } }
    virtual const VectorXf &GetVector() const override { return BaseType::GetValue(); }
};

class VarX : public VariableImpX<VectorXf> {
public:
    using BaseType = VariableImpX<VectorXf>;
    explicit VarX(const VectorXf &val) : BaseType(val.Size()) { BaseType::SetValue(val); }

    virtual void Project() override { return; }
    virtual void SetVector(const VectorXf &v) override { if (BaseType::CheckDim(v)) { BaseType::SetValue(v); } }
    virtual const VectorXf &GetVector() const override { return BaseType::GetValue(); }
};

} // namespace nlls