#include "../geometry/se2.h"

#include "vertex_se2.h"

namespace nlls {

void VertexSE2::Reset() {
    estimate_.SetIdentity();
}

void VertexSE2::ApplyPerturbation(const Vector3f &pert) {
    SE2fTangent left_update(pert[0], pert[1], pert[2]);
    estimate_ = left_update.Exp().ToTransform() * estimate_;
}

void VertexSE2Jet::Reset() {
    estimate_.SetIdentity();
    SetEstimate(estimate_);
}

void VertexSE2Jet::ApplyPerturbation(const Vector3j &pert) {
    SE2Tangent<Jetf> left_update(pert[0], pert[1], pert[2]);
    estimate_jet_ = left_update.Exp().ToTransform() * estimate_jet_;
}

NLLS_REGISTER_VERTEX(VertexSE2)
NLLS_REGISTER_VERTEX(VertexSE2Jet)

} // namespace nlls