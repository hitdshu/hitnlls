#pragma once

#include "nlls/graph.h"

namespace hitnlls {

class GraphSimple : public Graph {
public:
    virtual void SymbolicAnalysis() override;
};

} // namespace hitnlls