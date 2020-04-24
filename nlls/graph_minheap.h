#pragma once

#include "nlls/graph_simple.h"

namespace hitnlls {

class GraphMinheap : public GraphSimple {
public:
    virtual void SymbolicAnalysis() override;
};

} // namespace hitnlls