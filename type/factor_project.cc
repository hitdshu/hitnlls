#include "factor_project.h"

namespace nlls {

void FactorProjectPt3SE3::SetCameraParam(const Vector4f &intri) {
    fx_ = intri[0];
    fy_ = intri[1];
    cx_ = intri[2];
    cy_ = intri[3];
}

FactorProjectPt3SE3::EvalVectorTypeJet FactorProjectPt3SE3::operator()(FactorProjectPt3SE3::VertexTupleType &vars) {
    VertexR3Jet *vpt3 = vars.GetVertex<0>();
    VertexSE3Jet *vse3 = vars.GetVertex<1>();
    auto se3 = vse3->GetEstimateJet();
    auto pt3 = vpt3->GetEstimateJet();
    Vector3j pc;
    pc = se3.Block(0, 0, 3, 3) * pt3 + se3.Block(0, 3, 3, 1);
    EvalVectorTypeJet result;
    result << fx_ * pc[0] / pc[2] + cx_, fy_ * pc[1] / pc[2] + cy_;
    return result;
}

NLLS_REGISTER_FACTOR(FactorProjectPt3SE3)

} // namespace nlls