#pragma once

#include "node/base_node_impl.h"
#include "node/so3.h"

namespace hitnlls {
namespace node {

class NodeRot3 : public BaseNodeImpl<SO3, 3> {
public:
    virtual void UpdateInternal(const ::hitnlls::matrix::Vector3f &inc) override;

};

} // namespace node
} // namespace hitnlls