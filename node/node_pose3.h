#pragma once

#include "node/se3.h"
#include "node/base_node_impl.h"

namespace hitnlls {
namespace node {

class NodePose3 : public BaseNodeImpl<SE3, 6> {
public:
    virtual void UpdateInternal(const ::hitnlls::matrix::Matrix<float, 6, 1> &inc) override;

};

} // namespace node
} // namespace hitnlls