#include "../geometry/se3.h"

#include "vertex_se3.h"

namespace nlls {

void VertexSE3::Reset() {
    estimate_.SetIdentity();
}

void VertexSE3::ApplyPerturbation(const Vector6f &pert) {
    Vector3f t = pert.Block(0, 0, 3, 1);
    Vector3f r = pert.Block(3, 0, 3, 1);
    SE3fTangent left_update(t, r);
    estimate_ = left_update.Exp().ToTransform() * estimate_;
}

void VertexSE3Jet::Reset() {
    estimate_.SetIdentity();
    SetEstimate(estimate_);
}

void VertexSE3Jet::ApplyPerturbation(const Vector6j &pert) {
    Vector3j t = pert.Block(0, 0, 3, 1);
    Vector3j r = pert.Block(3, 0, 3, 1);
    SE3Tangent<Jetf> left_update(t, r);
    estimate_jet_ = left_update.Exp().ToTransform() * estimate_jet_;
}

NLLS_REGISTER_VERTEX(VertexSE3)
NLLS_REGISTER_VERTEX(VertexSE3Jet)

} // namespace nlls