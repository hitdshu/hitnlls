#pragma once

#include "../ils/factor.h"

#include "vertex_rn.h"

namespace nlls {

template <int M, int N>
class FactorLinear : public Factor<M, VertexRn<N>, VertexRn<M>> {
public:
    using ThisType = FactorLinear;
    using BaseType = Factor<M, VertexRn<N>, VertexRn<M>>;
    using NType = VertexRn<N>;
    using MType = VertexRn<M>;
    using JacobianMatrix = typename BaseType::JacobianMatrix;
    explicit FactorLinear() = default;
    explicit FactorLinear(const Matrix<float, M, N> &A, const Matrix<float, M, 1> &b) : A_(A) { BaseType::SetMeasurement(b); }

    void SetA(const Matrix<float, M, N> &A) { A_ = A; }

    virtual void ErrorAndJacobian(bool chi_only = false) override {
        auto &vt = BaseType::GetVertexTuple();
        NType *nv = vt.template GetVertex<0>();
        MType *mv = vt.template GetVertex<1>();
        BaseType::residual_ = BaseType::sqrt_info_ * (BaseType::measure_ - A_ * nv->GetEstimate() + mv->GetEstimate());
        BaseType::chi_ = pow(BaseType::residual_.Norm(), 2);
        if (chi_only) {
            return;
        }
        BaseType::jacobians_.Block(0, 0, M, N) = BaseType::sqrt_info_ * A_;
        BaseType::jacobians_.Block(0, N, M, M) = -BaseType::sqrt_info_;
    }

private:
    Matrix<float, M, N> A_;
};

} // namespace nlls