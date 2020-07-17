#include "solver.h"
#include "problem.h"
#include "dense_qr.h"
#include "dense_cholesky.h"
#include "dense_schur.h"
#include "sparse_cholesky.h"

namespace nlls {
namespace internal {

float Solver::ComputeChi() const {
    float chi = 0;
    for (const auto &f : problem_->factors_) {
        f->Compute(false);
        chi += f->GetChi();
    }
    return chi;
}

std::unordered_set<FactorBase *> &Solver::Factors() {
    return problem_->factors_;
}

const std::unordered_set<FactorBase *> &Solver::Factors() const {
    return problem_->factors_;
}

std::unordered_set<VertexBase *> &Solver::Vertices() {
    return problem_->vertices_;
}

const std::unordered_set<VertexBase *> &Solver::Vertices() const {
    return problem_->vertices_;
}

std::unordered_set<VertexBase *> &Solver::SchurVertices() {
    return problem_->schur_vertices_;
}

const std::unordered_set<VertexBase *> &Solver::SchurVertices() const {
    return problem_->schur_vertices_;
}

void Solver::Solve() {
    BuildStructure();
    const SolverOption &option = problem_->option_;
    float err_pre = ComputeChi();
    float err_cur = 0;
    {
        IterStat stat;
        stat.iter = 0;
        stat.chi = err_pre;
        stat.inc = 0;
        problem_->stat_.push_back(stat);
    }
    for (int iter = 0; iter < option.num_iters; ++iter) {
        float inc = SolveOneStep();
        err_cur = ComputeChi();
        IterStat stat;
        stat.iter = iter + 1;
        stat.chi = err_cur;
        stat.inc = inc;
        problem_->stat_.push_back(stat);
        if (((err_pre- err_cur) / err_pre) < option.term_ratio || err_cur < option.term_error) {
            break;
        }
        err_pre = err_cur;
    }
}

Solver *Solver::CreateSolver(Problem *p, Strategy s, Method m) {
    Solver *solver = nullptr;
    if (m == DENSE_QR) {
        solver = new DenseQr(s);
    } else if (m == DENSE_CHOLESKY) {
        solver = new DenseCholesky(s);
    } else if (m == DENSE_SCHUR) {
        solver = new DenseSchur(s);
    } else if (m == SPARSE_CHOLESKY) {
        solver = new SparseCholesky(s);
    }
    solver->SetProblem(p);
    return solver;
}

void Solver::PushVertices() {
    for (auto iter : problem_->vertices_) {
        iter->Push();
    }
}

void Solver::PopVertices() {
    for (auto iter : problem_->vertices_) {
        iter->Pop();
    }
}

float Solver::SolveOneStep() {
    const Strategy sttg = problem_->option_.sttg;
    if (sttg == GAUSS_NEWTON) {
        return SolveOneStepGaussNewton();
    } else if (sttg == LEVEN_MARQ) {
        return SolveOneStepLevenMarq();
    } else {
        return SolveOneStepTrustRegion();
    }
}

float Solver::SolveOneStepGaussNewton() {
    float err_bef = this->ComputeChi();
    FillStructure();
    VectorXf inc = this->ComputeGaussNewtonPoint();
    PushVertices();
    this->Update(inc);
    float inc_norm = inc.Norm();
    float err_cur = this->ComputeChi();
    if (err_cur > err_bef) {
        PopVertices();
    }
    return inc_norm;
}

float Solver::SolveOneStepLevenMarq() {
    const SolverOption &option = problem_->option_;
    FillStructure();
    if (lmp_.u < 0) {
        lmp_.u = InitLevenMarqParam();
    }
    float err_bef = this->ComputeChi();
    float inc_norm = 0;
    for (int i = 0; i < option.num_inner_iters; ++i) {
        float l_off = 0;
        if (i != 0) {
            l_off = lmp_.u / lmp_.v * 2;
        }
        VectorXf inc = this->ComputeGaussNewtonPoint(lmp_.u - l_off);
        inc_norm = inc.Norm();
        PushVertices();
        this->Update(inc);
        float err_cur = this->ComputeChi();
        if (err_cur > err_bef) {
            PopVertices();
            lmp_.u *= lmp_.v;
            lmp_.v *= 2;
        } else {
            lmp_.u *= 1.0 / 3.0;
            lmp_.v = 2;
            break;
        }
    }
    return inc_norm;
}

float Solver::SolveOneStepTrustRegion() {
    const SolverOption &option = problem_->option_;
    FillStructure();
    float err_bef = this->ComputeChi();
    VectorXf step_gn = ComputeGaussNewtonPoint();
    VectorXf step_cc = ComputeCauchyPoint();
    float gn_len = step_gn.Norm();
    float cc_len = step_cc.Norm();
    if (trp_.r < 0) {
        trp_.r = step_gn.Norm() / 2.0;
    }
    VectorXf update_step;
    for (int i = 0; i < option.num_inner_iters; ++i) {
        if (gn_len < trp_.r) {
            update_step = step_gn;
        } else if (cc_len < trp_.r) {
            float int_coeff = (gn_len - trp_.r) / (gn_len - cc_len);
            update_step = int_coeff * step_cc + (1 - int_coeff) * step_gn;
        } else {
            update_step = step_cc / cc_len * trp_.r;
        }
        PushVertices();
        Update(update_step);
        float err_cur = this->ComputeChi();
        if (err_cur > err_bef) {
            PopVertices();
            trp_.r /= trp_.v;
        } else {
            float ratio = SurrogateReduction(update_step) / (err_bef - err_cur);
            if (ratio < trp_.l) {
                trp_.r /= trp_.v;
            } else if (ratio > trp_.u) {
                trp_.r *= trp_.v;
            }
            break;
        }
    }
    return update_step.Norm();
}

} // namespace internal
} // namespace nlls