#pragma once

#include "solver/base_solver.h"

namespace hitnlls {
namespace solver {

class GnSolver : public BaseSolver {
public:
    GnSolver() : BaseSolver() {}

    virtual void Optimize() override;
};

} // namespace solver
} // namespace hitnlls