#include "type/vertex_pt3.h"

namespace hitnlls {

void VertexPt3::SetZero() {
    estimate_.SetZero();
}

void VertexPt3::ApplyPerturbation(const matrix::Vector3d &pert) {
    estimate_ += pert;
}

HITNLLS_REGISTER_VERTEX(VertexPt3)
HITNLLS_REGISTER_VERTEX(VertexPt3Jet)

} // namespace hitnlls