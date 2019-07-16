#pragma once

#include "node/base_node.h"
#include "factor/base_factor.h"
#include "matrix/matrixs.h"
#include "matrix/vecxs.h"
#include <map>
#include <vector>

namespace hitnlls {
namespace graph {

class BaseGraph {
public:
    BaseGraph() = default;
    virtual ~BaseGraph() = default;

    virtual void AddNode(::hitnlls::node::BaseNode *node) final { nodes_[node->GetId()] = node; }
    virtual void AddFactor(::hitnlls::factor::BaseFactor *factor) final { factors_.push_back(factor); }

    virtual void BuildProblem(::hitnlls::matrix::Matrixs<::hitnlls::matrix::Matrixxf> &matA, ::hitnlls::matrix::Vecxs<::hitnlls::matrix::Matrixxf> &vecb) = 0;
    virtual bool UpdateInc(const ::hitnlls::matrix::Vecxs<::hitnlls::matrix::Matrixxf> &inc) = 0;
    virtual float ComputeGraphError() = 0;

    BaseGraph(const BaseGraph &) = delete;
    BaseGraph &operator=(const BaseGraph &) = delete;

protected:
    ::std::map<int, ::hitnlls::node::BaseNode *> nodes_;
    ::std::vector<::hitnlls::factor::BaseFactor *> factors_;
};

HITNLLS_REGISTER_REGISTERER(BaseGraph);
#define HITNLLS_REGISTER_GRAPH(name) \
    HITNLLS_REGISTER_CLASS(BaseGraph, name)

} // namespace graph
} // namespace hitnlls