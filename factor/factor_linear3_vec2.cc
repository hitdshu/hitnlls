#include "factor/factor_linear3_vec2.h"
#include "node/node_vec2.h"

namespace hitnlls {
namespace factor {

::hitnlls::matrix::Matrix<float, 3, 1> FactorLinear3Vec2::Evaluate() {
    ::hitnlls::node::NodeVec2 *node_v2 = dynamic_cast<::hitnlls::node::NodeVec2 *>(GetNode(0));

    ::hitnlls::matrix::Vector2f v2 = node_v2->GetEstimate();
    return matA_ * v2;
}

::hitnlls::matrix::Matrixxf FactorLinear3Vec2::Jacobian(int nidx) {
    return matA_;
}

HITNLLS_REGISTER_FACTOR(FactorLinear3Vec2);

} // namespace factor
} // namespace hinlls