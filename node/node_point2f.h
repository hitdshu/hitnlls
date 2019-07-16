#pragma once

#include "node/base_node_impl.h"

namespace hitnlls {
namespace node {

class NodePoint2f : public BaseNodeImpl<::hitnlls::matrix::Vector2f, 2> {
public:
    virtual void UpdateInternal(const ::hitnlls::matrix::Vector2f &inc) override;

};

} // namespace node
} // namespace hitnlls