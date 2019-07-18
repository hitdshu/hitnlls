#pragma once

#include "node/base_node_impl.h"

namespace hitnlls {
namespace node {

class NodeVec3 : public BaseNodeImpl<::hitnlls::matrix::Vector3f, 3> {
public:
    virtual void UpdateInternal(const ::hitnlls::matrix::Vector3f &inc) override;

};

} // namespace node
} // namespace hitnlls