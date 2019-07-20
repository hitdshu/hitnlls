#pragma once

#include "graph/base_graph.h"
#include "matrix/matrixx.h"
#include <set>

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
    virtual void MarginalizationAnalysis() final {
        ::std::set<int> marg_node_ids;
        ::std::set<int> nomg_node_ids;
        ::std::set<int> factor_node_ids;
        for (auto iter = nodes_.begin(); iter != nodes_.end(); ++iter) {
            if (!(iter->second)->GetFixed()) {
                if ((iter->second)->GetMarginalized()) {
                    marg_node_ids.insert(iter->first);
                } else {
                    nomg_node_ids.insert(iter->first);
                }
            }
        }
        for (size_t idx = 0; idx < factors_.size(); ++idx) {
            for (int nidx = 0; nidx < factors_[idx]->GetNnodes(); ++nidx) {
                int nid = factors_[idx]->GetNodeId(nidx);
                factor_node_ids.insert(nid);
            }
        }
        ::std::set<int> effe_marg_node_ids;
        ::std::set<int> effe_nomg_node_ids;
        for (auto iter = marg_node_ids.begin(); iter != marg_node_ids.end(); ++iter) {
            if (1 == factor_node_ids.count(*iter)) {
                effe_marg_node_ids.insert(*iter);
            }
        }
        for (auto iter = nomg_node_ids.begin(); iter != nomg_node_ids.end(); ++iter) {
            if (1 == factor_node_ids.count(*iter)) {
                effe_nomg_node_ids.insert(*iter);
            }
        }
        ordering_ = ::hitnlls::matrix::Matrixxi(effe_marg_node_ids.size() + effe_nomg_node_ids.size(), 1);
        int ridx = 0;
        for (auto iter = effe_marg_node_ids.begin(); iter != effe_marg_node_ids.end(); ++iter) {
            ordering_[ridx++] = *iter;
        }
        for (auto iter = effe_nomg_node_ids.begin(); iter != effe_nomg_node_ids.end(); ++iter) {
            ordering_[ridx++] = *iter;
        }
    }
    virtual void BuildProblem(::hitnlls::matrix::Matrixs<::hitnlls::matrix::Matrixxf> &matA, ::hitnlls::matrix::Vecxs<::hitnlls::matrix::Matrixxf> &vecb) override final {
        if (!HasMarginalization()) {
            SymbolicAnalysis();
        } else {
            MarginalizationAnalysis();
        }
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