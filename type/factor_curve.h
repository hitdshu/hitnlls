#pragma once

#include "../ils/factor.h"

#include "vertex_rn.h"

namespace nlls {

class FactorCurveExp : public FactorAutoDiff<1, VertexR2Jet> {
public:
    using BaseType = FactorAutoDiff<1, VertexR2Jet>;
    using VertexTupleType = typename BaseType::VertexTupleType;
    using EvalVectorTypeJet = typename BaseType::EvalVectorTypeJet;
    explicit FactorCurveExp(float x = 0) : x_(x) {}
    void SetX(float x) { x_ = x; }
    virtual EvalVectorTypeJet operator()(VertexTupleType &vars) override;

private:
    float x_;
};

} // namespace nlls