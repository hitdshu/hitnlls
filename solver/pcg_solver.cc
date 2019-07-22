#include "solver/pcg_solver.h"
#include <limits>

namespace hitnlls {
namespace solver {

void PcgSolver::Optimize() {
    float error_bef = std::numeric_limits<float>::max();
    float error_aft = std::numeric_limits<float>::max();
    for (int sidx = 0; sidx < steps_; ++sidx) {
        ::hitnlls::matrix::Matrixs<::hitnlls::matrix::Matrixxf> mat_a;
        ::hitnlls::matrix::Vecxs<::hitnlls::matrix::Matrixxf> vecb;
        error_bef = graph_->ComputeGraphError();
        graph_->BuildProblem(mat_a, vecb);
        ::hitnlls::matrix::Matrixs<::hitnlls::matrix::Matrixxf> m_inv = mat_a.GetPreconditioner();
        ::hitnlls::matrix::Vecxs<::hitnlls::matrix::Matrixxf> r = vecb;
        ::hitnlls::matrix::Vecxs<::hitnlls::matrix::Matrixxf> z = m_inv * r;
        ::hitnlls::matrix::Vecxs<::hitnlls::matrix::Matrixxf> p = z;
        ::hitnlls::matrix::Vecxs<::hitnlls::matrix::Matrixxf> w = mat_a * p;
        ::hitnlls::matrix::Matrixxf alpha = (r * z) / (p * w);
        ::hitnlls::matrix::Vecxs<::hitnlls::matrix::Matrixxf> x = p * alpha;
        ::hitnlls::matrix::Vecxs<::hitnlls::matrix::Matrixxf> rn = r - w * alpha;
        double init_res = r.Norm();
        double comp_res = ::std::numeric_limits<double>::max();
        int iter_idx = 0;

        while (iter_idx < vecb.Length() && comp_res > init_res * 1e-4) {
            ::hitnlls::matrix::Vecxs<::hitnlls::matrix::Matrixxf> zn = m_inv * r;
            ::hitnlls::matrix::Matrixxf beta = (rn * zn) / (r * z);
            p = zn + p * beta;
            w = mat_a * p;
            alpha = (rn * zn) / (p * w);
            x = x + p * alpha;
            r = rn;
            rn = rn - w * alpha;
            z = zn;
            comp_res = r.Norm();
            ++iter_idx;
        }

        error_aft = graph_->ComputeGraphError();
        if (early_stop_) {
            if (fabs(error_aft - error_bef) / error_bef < error_reduction_) {
                return;
            }
        }
    }
}

HITNLLS_REGISTER_SOLVER(PcgSolver);

} // namespace solver
} // namespace hitnlls