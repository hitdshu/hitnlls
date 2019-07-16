#include "node/node_rot2.h"

namespace hitnlls {
namespace node {

void NodeRot2::UpdateInternal(const ::hitnlls::matrix::Vector1f &inc) {
    est_ += inc;
}

HITNLLS_REGISTER_NODE(NodeRot2);

} // namespace node
} // namespace hitnlls
