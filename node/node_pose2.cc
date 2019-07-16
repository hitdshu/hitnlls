#include "node/node_pose2.h"

namespace hitnlls {
namespace node {

void NodePose2::UpdateInternal(const ::hitnlls::matrix::Vector3f &inc) {
    est_ += inc;
}

HITNLLS_REGISTER_NODE(NodePose2);

} // namespace node
} // namespace hitnlls
