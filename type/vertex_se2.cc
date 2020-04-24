#include "type/vertex_se2.h"

namespace hitnlls {

void VertexSE2::SetZero() {
    estimate_.SetIdentity();
}

void VertexSE2::ApplyPerturbation(const matrix::Vector3d &pert) {
    geometry::SE2dTangent left_update(pert[0], pert[1], pert[2]);
    estimate_ = left_update.Exp().ToTransform() * estimate_;
}

HITNLLS_REGISTER_VERTEX(VertexSE2)
HITNLLS_REGISTER_VERTEX(VertexSE2Jet)

} // namespace hitnlls