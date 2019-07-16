#include "solver/gn_solver.h"
#include <limits>

namespace hitnlls {
namespace solver {

void GnSolver::Optimize() {
    float error_bef = std::numeric_limits<float>::max();
    float error_aft = std::numeric_limits<float>::max();
    for (int sidx = 0; sidx < steps_; ++sidx) {
        ::hitnlls::matrix::Matrixs<::hitnlls::matrix::Matrixxf> matA;
        ::hitnlls::matrix::Vecxs<::hitnlls::matrix::Matrixxf> vecb;
        error_bef = graph_->ComputeGraphError();
        graph_->BuildProblem(matA, vecb);

        ::hitnlls::matrix::Vecxs<::hitnlls::matrix::Matrixxf> inc = matA.SolveWithlm(vecb);
        if (!graph_->UpdateInc(inc)) {
            return;
        }
        error_aft = graph_->ComputeGraphError();
        if (early_stop_) {
            if (fabs(error_aft - error_bef) / error_bef < error_reduction_) {
                return;
            }
        }
    }
}

HITNLLS_REGISTER_SOLVER(GnSolver);

} // namespace solver
} // namespace hitnlls