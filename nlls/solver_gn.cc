#include <limits>

#include "nlls/solver_gn.h"

namespace hitnlls {

void SolverGn::Optimize() {
    double error_bef = ::std::numeric_limits<double>::max();
    double error_aft = ::std::numeric_limits<double>::max();
    for (int sidx = 0; sidx < steps_; ++sidx) {
        SparseBlockMatrix mat_a;
        MatrixSparseArray vecb;
        error_bef = graph_->ComputeGraphChi();
        graph_->BuildProblem(mat_a, vecb);
        MatrixSparseArray inc = mat_a.SolveWithLm(vecb);
        if (!graph_->UpdateInc(inc)) {
            return;
        }
        error_aft = graph_->ComputeGraphChi();
        if (early_stop_) {
            if (std::abs(error_aft - error_bef) / error_bef < error_reduction_) {
                return;
            }
        }
    }
}

HITNLLS_REGISTER_SOLVER(SolverGn)

} // namespace hitnlls