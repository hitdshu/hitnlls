#pragma once

#include "node/se2.h"
#include "node/base_node_impl.h"

namespace hitnlls {
namespace node {

class NodePose2 : public BaseNodeImpl<SE2, 3> {
public:
    virtual void UpdateInternal(const ::hitnlls::matrix::Vector3f &inc) override;

};

} // namespace node
} // namespace hitnlls