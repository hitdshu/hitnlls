#pragma once

#include "solver/base_solver.h"

namespace hitnlls {
namespace solver {

class LmSolver : public BaseSolver {
public:
    LmSolver() : BaseSolver() {}

    virtual void Optimize() override;
};

} // namespace solver
} // namespace hitnlls