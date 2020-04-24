#include "type/vertex_se3.h"

namespace hitnlls {

void VertexSE3::SetZero() {
    estimate_.SetIdentity();
}

void VertexSE3::ApplyPerturbation(const matrix::Vector6d &pert) {
    matrix::Vector3d t = pert.Block(0, 0, 3, 1);
    matrix::Vector3d r = pert.Block(3, 0, 3, 1);
    geometry::SE3dTangent left_update(t, r);
    estimate_ = left_update.Exp().ToTransform() * estimate_;
}

HITNLLS_REGISTER_VERTEX(VertexSE3)
HITNLLS_REGISTER_VERTEX(VertexSE3Jet)

} // namespace hitnlls