#include "factor_curve.h"

namespace nlls {

typename FactorCurveExp::EvalVectorTypeJet FactorCurveExp::operator()(typename FactorCurveExp::VertexTupleType &vars) {
    EvalVectorTypeJet result;
    auto var = vars.GetVertex<0>();
    auto est = var->GetEstimateJet();
    result = exp(est[0] + Jetf(x_) * est[1]);
    return result;
}

NLLS_REGISTER_FACTOR(FactorCurveExp)

} // namespace nlls