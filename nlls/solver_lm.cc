#include <limits>
#include <algorithm>

#include "nlls/solver_lm.h"

namespace hitnlls {

void SolverLm::Optimize() {
    const int lm_iter_num = 6;
    double error_bef = std::numeric_limits<double>::max();
    double error_aft = std::numeric_limits<double>::max();
    double lambda = 0;
    for (int sidx = 0; sidx < steps_; ++sidx) {
        SparseBlockMatrix mat_a;
        MatrixSparseArray vecb;
        error_bef = graph_->ComputeGraphChi();
        graph_->BuildProblem(mat_a, vecb);
        if (0 == sidx) {
            for (int ridx = 0; ridx < mat_a.Rows(); ++ridx) {
                lambda = std::max(lambda, mat_a(ridx, ridx)(0, 0));
            }
            lambda *= 1e-4;
        }
        mat_a.SetLmLambda(lambda);
        MatrixSparseArray inc = mat_a.SolveWithLm(vecb);
        int lm_idx = 0;
        while (!graph_->UpdateInc(inc) && lm_idx < lm_iter_num) {
            lambda *= 10;
            mat_a.SetLmLambda(lambda);
            inc = mat_a.SolveWithLm(vecb);
            ++lm_idx;
        }
        if (lm_iter_num == lm_idx) {
            return;
        }
        lambda /= 10;
        error_aft = graph_->ComputeGraphChi();
        if (early_stop_) {
            if (std::abs(error_aft - error_bef) / error_bef < error_reduction_) {
                return;
            }
        }
    }
}

HITNLLS_REGISTER_SOLVER(SolverLm)

} // namespace hitnlls