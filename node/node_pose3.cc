#include "node/node_pose3.h"

namespace hitnlls {
namespace node {

void NodePose3::UpdateInternal(const ::hitnlls::matrix::Matrix<float, 6, 1> &inc) {
    est_ += inc;
}

HITNLLS_REGISTER_NODE(NodePose3);

} // namespace node
} // namespace hitnlls
