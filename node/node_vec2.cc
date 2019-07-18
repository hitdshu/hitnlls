#include "node/node_vec2.h"

namespace hitnlls {
namespace node {

void NodeVec2::UpdateInternal(const ::hitnlls::matrix::Vector2f &inc) {
    est_ += inc;
}

HITNLLS_REGISTER_NODE(NodeVec2);

} // namespace node
} // namespace hitnlls
