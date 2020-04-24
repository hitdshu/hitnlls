#pragma once

#include "nlls/solver.h"

namespace hitnlls {

class SolverGn : public SolverBase {
public:
    SolverGn() : SolverBase() {}

    virtual void Optimize() override;
};

} // namespace hitnlls