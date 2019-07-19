#pragma once

#include "graph/naive_graph.h"

namespace hitnlls {
namespace graph {

class MinheapGraph : public NaiveGraph {
public:
    virtual void SymbolicAnalysis() override;
};

} // namespace graph
} // namespace hitnlls