#pragma once

#include "type/factor.h"
#include "type/vertex_pt2.h"
#include "type/vertex_se2.h"

namespace hitnlls {

class FactorPt2SE2 : public FactorAutoDiff<2, VertexPt2Jet, VertexSE2Jet> {
public:
    using ThisType = FactorPt2SE2;
    using BaseType = FactorAutoDiff<2, VertexPt2Jet, VertexSE2Jet>;
    using VertexTupleType = typename BaseType::VertexTupleType;
    using EvalVectorTypeJet = typename BaseType::EvalVectorTypeJet;

    virtual EvalVectorTypeJet operator()(VertexTupleType &vars) override final;
};

} // namespace hitnlls