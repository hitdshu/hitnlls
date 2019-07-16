#pragma once

#include <iostream>
#include "node/base_node.h"
#include "matrix/matrix.h"

namespace hitnlls {
namespace node {

template<typename ValueType, int ndim>
class BaseNodeImpl : public BaseNode {
public:
    BaseNodeImpl() : BaseNode(ndim) {}
    virtual ~BaseNodeImpl() = default;

    virtual void SetEstimate(const ValueType &est) final { est_ = est; }
    virtual ValueType GetEstimate() const final { return est_; }
    virtual void SaveToCache() override final { cache_ = est_; }
    virtual void RestoreCache() override final { est_ = cache_; }

    virtual void UpdateInternal(const ::hitnlls::matrix::Matrix<float, ndim, 1> &inc) = 0;

    virtual void UpdateImp(const ::hitnlls::matrix::Matrixxf &inc) override final {
        UpdateInternal(inc);
    }

protected:
    ValueType est_;
    ValueType cache_;
};

} // namespace node
} // namespace hitnlls