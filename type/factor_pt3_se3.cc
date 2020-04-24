#include "type/factor_pt3_se3.h"

namespace hitnlls {

FactorPt3SE3::EvalVectorTypeJet FactorPt3SE3::operator()(FactorPt3SE3::VertexTupleType &vars) {
    VertexPt3Jet *vpt3 = vars.GetVertex<0>();
    VertexSE3Jet *vse3 = vars.GetVertex<1>();
    auto se3 = vse3->GetJetEstimate();
    auto pt3 = vpt3->GetJetEstimate();
    JetVector3d pc;
    pc = se3.Block(0, 0, 3, 3) * pt3 + se3.Block(0, 3, 3, 1);
    EvalVectorTypeJet result;
    result = camera_->Project(pc);
    return result;
}

HITNLLS_REGISTER_FACTOR(FactorPt3SE3)

} // namespace hitnlls