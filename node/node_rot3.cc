#include "node/node_rot3.h"

namespace hitnlls {
namespace node {

void NodeRot3::UpdateInternal(const ::hitnlls::matrix::Vector3f &inc) {
    est_ += inc;
}

HITNLLS_REGISTER_NODE(NodeRot3);

} // namespace node
} // namespace hitnlls
