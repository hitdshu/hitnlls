#include <sstream>
#include <iomanip>

#include "problem.h"

namespace nlls {

Problem::~Problem() {
    std::unordered_set<VertexBase *> all_vertices;
    std::unordered_set<KernelBase *> all_kernels;
    for (auto iter : factors_) {
        for (int i = 0; i < iter->NumVertices(); ++i) {
            VertexBase *tmp = iter->GetVertexAt(i);
            all_vertices.insert(tmp);
        }
        if (iter->GetKernel()) {
            all_kernels.insert(iter->GetKernel());
        }
        delete iter;
    }
    for (auto iter : all_vertices) {
        delete iter;
    }
    for (auto iter : all_kernels) {
        delete iter;
    }
}

float Problem::Solve() {
    this->BuildStructure();
    internal::Solver *solver = internal::Solver::CreateSolver(this, option_.sttg, option_.method);
    solver->Solve();
    float chi = solver->ComputeChi();
    delete solver;
    return chi;
}

void Problem::BuildStructure() {
    for (auto iter = factors_.begin(); iter != factors_.end(); ++iter) {
        if ((*iter)->GetStatus() == FactorBase::NonActive) {
            delete *iter;
            iter = factors_.erase(iter);
        }
    }
    for (auto iter = factors_.begin(); iter != factors_.end(); ++iter) {
        for (int i = 0; i < (*iter)->NumVertices(); ++i) {
            VertexBase *tmp = (*iter)->GetVertexAt(i);
            if (tmp->GetStatus() == VertexBase::Active) {
                vertices_.insert(tmp);
            }
        }
    }
}

std::string Problem::BriefReport() const {
    const IterStat &ss = GetInitStat();
    const IterStat &es = GetFinalStat();
    std::stringstream strs;
    strs << "-------------------------------" << "Solver Report" << "-------------------------------" << std::endl;
    strs << "Iter:" << std::setw(4) << ss.iter << ", chi:" << std::setw(14) << ss.chi << ", inc: " << std::setw(14) << ss.inc << std::endl;
    strs << "Iter:" << std::setw(4) << es.iter << ", chi:" << std::setw(14) << es.chi << ", inc: " << std::setw(14) << ss.inc << std::endl;
    strs << "-------------------------------" << "Solver Report" << "-------------------------------";
    return strs.str();
}

std::string Problem::Report() const {
    std::stringstream strs;
    strs << "-------------------------------" << "Solver Report" << "-------------------------------" << std::endl;
    for (size_t i = 0; i < stat_.size(); ++i) {
        strs << "Iter:" << std::setw(4) << stat_[i].iter << ", chi:" << std::setw(14) << stat_[i].chi << ", inc: " << std::setw(14) << stat_[i].inc << std::endl;
    }
    strs << "-------------------------------" << "Solver Report" << "-------------------------------";
    return strs.str();
}

} // namespace nlls