#pragma once

#include "solver/base_solver.h"

namespace hitnlls {
namespace solver {

class PcgSolver : public BaseSolver {
public:
    PcgSolver() : BaseSolver() {}

    virtual void Optimize() override;
};

} // namespace solver
} // namespace hitnlls