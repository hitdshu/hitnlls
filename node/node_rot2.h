#pragma once

#include "node/base_node_impl.h"
#include "node/so2.h"

namespace hitnlls {
namespace node {

class NodeRot2 : public BaseNodeImpl<SO2, 1> {
public:
    virtual void UpdateInternal(const ::hitnlls::matrix::Vector1f &inc) override;

};

} // namespace node
} // namespace hitnlls