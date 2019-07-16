#include "node/node_point2f.h"

namespace hitnlls {
namespace node {

void NodePoint2f::UpdateInternal(const ::hitnlls::matrix::Vector2f &inc) {
    est_ += inc;
}

HITNLLS_REGISTER_NODE(NodePoint2f);

} // namespace node
} // namespace hitnlls
