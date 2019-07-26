#include "node/node_point3f.h"

namespace hitnlls {
namespace node {

void NodePoint3f::UpdateInternal(const ::hitnlls::matrix::Vector3f &inc) {
    est_ += inc;
}

HITNLLS_REGISTER_NODE(NodePoint3f);

} // namespace node
} // namespace hitnlls
