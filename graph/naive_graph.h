#pragma once

#include "graph/base_graph_impl.h"

namespace hitnlls {
namespace graph {

class NaiveGraph : public BaseGraphImpl {
public:
    virtual void SymbolicAnalysis() override;
};

} // namespace graph
} // namespace hitnlls