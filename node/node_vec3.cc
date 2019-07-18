#include "node/node_vec3.h"

namespace hitnlls {
namespace node {

void NodeVec3::UpdateInternal(const ::hitnlls::matrix::Vector3f &inc) {
    est_ += inc;
}

HITNLLS_REGISTER_NODE(NodeVec3);

} // namespace node
} // namespace hitnlls
