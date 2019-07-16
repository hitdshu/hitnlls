#pragma once

#include "graph/base_graph.h"
#include "matrix/matrixx.h"
#include "matrix/matrixsym.h"

namespace hitnlls {
namespace graph {

class BaseGraphImpl : public BaseGraph {
public:
    virtual void SymbolicAnalysis() = 0;

    virtual float ComputeGraphError() override final {
        float error = 0;
        for (size_t idx = 0; idx < factors_.size(); ++idx) {
            error += factors_[idx]->ComputeError();
        }
        return error;
    }
    virtual void BuildProblem(::hitnlls::matrix::Matrixs<::hitnlls::matrix::Matrixxf> &matA, ::hitnlls::matrix::Vecxs<::hitnlls::matrix::Matrixxf> &vecb) override final {
        SymbolicAnalysis();
        matA = ::hitnlls::matrix::Matrixs<::hitnlls::matrix::Matrixxf>(ordering_.Rows(), ordering_.Rows());
        vecb = ::hitnlls::matrix::Vecxs<::hitnlls::matrix::Matrixxf>(ordering_.Rows());
        std::map<int, int> id2ordering;
        for (int idx = 0; idx < ordering_.Rows(); ++idx) {
            id2ordering[ordering_[idx]] = idx;
        }
        for (size_t idx = 0; idx < factors_.size(); ++idx) {
            for (int nidx1 = 0; nidx1 < factors_[idx]->GetNnodes(); ++nidx1) {
                int nid1 = factors_[idx]->GetNodeId(nidx1);
                if (0 == id2ordering.count(nid1))
                    continue;
                for (int nidx2 = 0; nidx2 < factors_[idx]->GetNnodes(); ++nidx2) {
                    int nid2 = factors_[idx]->GetNodeId(nidx2);
                    if (0 == id2ordering.count(nid2))
                        continue;
                    matA(id2ordering[nid1], id2ordering[nid2]) += factors_[idx]->GetJacobTInfoJacob(nidx1, nidx2);
                }
                vecb[id2ordering[nid1]] += factors_[idx]->GetJacobTInfoRes(nidx1);
            }
        }
    }
    virtual bool UpdateInc(const ::hitnlls::matrix::Vecxs<::hitnlls::matrix::Matrixxf> &inc) override final {
        float error_bef_update = 0;
        for (size_t idx = 0; idx < factors_.size(); ++idx) {
            error_bef_update += factors_[idx]->ComputeError();
        }
        for (size_t idx = 0; idx < ordering_.Rows(); ++idx) {
            nodes_[ordering_[idx]]->SaveToCache();
            nodes_[ordering_[idx]]->UpdatePlus(inc[idx]);
        }
        float error_aft_update = 0;
        for (size_t idx = 0; idx < factors_.size(); ++idx) {
            error_aft_update += factors_[idx]->ComputeError();
        }
        if (error_aft_update < error_bef_update) {
            return true;
        } else {
            for (size_t idx = 0; idx < ordering_.Rows(); ++idx) {
                nodes_[ordering_[idx]]->RestoreCache();
            }
            return false;
        }
        return false;
    }

protected:
    ::hitnlls::matrix::Matrixxi ordering_;
};

} // namespace graph
} // namespace hitnlls