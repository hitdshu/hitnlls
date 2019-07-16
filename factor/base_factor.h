#pragma once

#include "matrix/matrixx.h"
#include "node/base_node.h"
#include "common/register.h"
#include <map>

namespace hitnlls {
namespace factor {

class BaseFactor {
public:
    BaseFactor(int nnodes = -1, int ndim = -1) { nnodes_ = nnodes; ndim_ = ndim; }
    virtual ~BaseFactor() = default;

    virtual void SetNode(int nidx, ::hitnlls::node::BaseNode *node) final { nodes_[nidx] = node; }
    virtual ::hitnlls::node::BaseNode * GetNode(int nidx) final { return nodes_[nidx]; }
    virtual int GetDim() const final { return ndim_; }
    virtual int GetNnodes() const final { return nnodes_; }
    virtual int GetNodeId(int nidx) final { return nodes_[nidx]->GetId(); }

    virtual float ComputeError() = 0;
    virtual ::hitnlls::matrix::Matrixxf GetJacobTInfoRes(int nidx) = 0;
    virtual ::hitnlls::matrix::Matrixxf GetJacobTInfoJacob(int nidx1, int nidx2) = 0;

    BaseFactor(const BaseFactor &) = delete;
    BaseFactor &operator=(const BaseFactor &) = delete;

protected:
    ::std::map<int, ::hitnlls::node::BaseNode *> nodes_;
    int nnodes_;
    int ndim_;
};

HITNLLS_REGISTER_REGISTERER(BaseFactor);
#define HITNLLS_REGISTER_FACTOR(name) \
    HITNLLS_REGISTER_CLASS(BaseFactor, name)

} // namespace factor
} // namespace hitnlls