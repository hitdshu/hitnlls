#pragma once

#include "nlls/solver.h"

namespace hitnlls {

class SolverLm : public SolverBase {
public:
    SolverLm() : SolverBase() {}

    virtual void Optimize() override;
};

} // namespace hitnlls