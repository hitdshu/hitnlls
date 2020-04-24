#include "type/vertex_pt2.h"

namespace hitnlls {

void VertexPt2::SetZero() {
    estimate_.SetZero();
}

void VertexPt2::ApplyPerturbation(const matrix::Vector2d &pert) {
    estimate_ += pert;
}

HITNLLS_REGISTER_VERTEX(VertexPt2)
HITNLLS_REGISTER_VERTEX(VertexPt2Jet)

} // namespace hitnlls