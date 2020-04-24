#pragma once

#include "type/factor.h"
#include "type/vertex_pt3.h"
#include "type/vertex_se3.h"
#include "type/camera_pinhole.h"

namespace hitnlls {

class FactorPt3SE3 : public FactorAutoDiff<2, VertexPt3Jet, VertexSE3Jet> {
public:
    using ThisType = FactorPt3SE3;
    using BaseType = FactorAutoDiff<2, VertexPt3Jet, VertexSE3Jet>;
    using VertexTupleType = typename BaseType::VertexTupleType;
    using EvalVectorTypeJet = typename BaseType::EvalVectorTypeJet;

    void SetCamera(CameraPinhole *camera) { camera_ = camera; }
    virtual EvalVectorTypeJet operator()(VertexTupleType &vars) override final;

private:
    CameraPinhole *camera_;
};

} // namespace hitnlls