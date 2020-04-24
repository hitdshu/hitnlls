#include "type/factor_pt2_se2.h"

namespace hitnlls {

FactorPt2SE2::EvalVectorTypeJet FactorPt2SE2::operator()(FactorPt2SE2::VertexTupleType &vars) {
    VertexPt2Jet *vpt2 = vars.GetVertex<0>();
    VertexSE2Jet *vse2 = vars.GetVertex<1>();
    auto se2 = vse2->GetJetEstimate();
    auto pt2 = vpt2->GetJetEstimate();
    EvalVectorTypeJet result;
    result = se2.Block(0, 0, 2, 2) * pt2 + se2.Block(0, 2, 2, 1);
    return result;
}

HITNLLS_REGISTER_FACTOR(FactorPt2SE2)

} // namespace hitnlls