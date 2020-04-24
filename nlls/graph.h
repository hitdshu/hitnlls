#pragma once

#include <map>
#include <vector>
#include <set>

#include "type/factor.h"
#include "type/vertex.h"
#include "common/register.h"
#include "nlls/sparse_array.h"
#include "nlls/sparse_matrix.h"

namespace hitnlls {

class GraphBase {
public:
    GraphBase() = default;
    virtual ~GraphBase() = default;

    void AddVertex(VertexBase *v) { id2v_[v->GetId()] = v; }
    void AddFactor(FactorBase *f) { factors_.push_back(f); }
    bool HasMarginalization() const {
        for (auto iter = id2v_.begin(); iter != id2v_.end(); ++iter) {
            if (iter->second->GetMarginalized()) {
                return true;
            }
        }
        return false;
    }

    virtual void BuildProblem(SparseBlockMatrix &A, MatrixSparseArray &b) = 0;
    virtual bool UpdateInc(const MatrixSparseArray &inc) = 0;
    virtual double ComputeGraphChi() const = 0;

    GraphBase(const GraphBase &) = delete;
    GraphBase &operator=(const GraphBase &) = delete;
    GraphBase &operator=(GraphBase &&) = delete;

protected:
    std::map<int, VertexBase *> id2v_;
    std::vector<FactorBase *> factors_;
};

HITNLLS_REGISTER_REGISTER(GraphBase)
#define HITNLLS_REGISTER_GRAPH(name) \
    HITNLLS_REGISTER_CLASS(GraphBase, name)

class Graph : public GraphBase {
public:
    void MarginalizationAnalysis() {
        std::set<int> marg_node_ids;
        std::set<int> nomg_node_ids;
        std::set<int> factor_node_ids;
        for (auto iter = id2v_.begin(); iter != id2v_.end(); ++iter) {
            if ((iter->second)->GetStatus() == VertexBase::Active) {
                if ((iter->second)->GetMarginalized()) {
                    marg_node_ids.insert(iter->first);
                } else {
                    nomg_node_ids.insert(iter->first);
                }
            }
        }
        for (size_t idx = 0; idx < factors_.size(); ++idx) {
            for (int nidx = 0; nidx < factors_[idx]->NumVertices(); ++nidx) {
                int nid = factors_[idx]->GetVertexIdAt(nidx);
                factor_node_ids.insert(nid);
            }
        }
        std::set<int> effe_marg_node_ids;
        std::set<int> effe_nomg_node_ids;
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
        ordering_ = SparseMatrixi(effe_marg_node_ids.size() + effe_nomg_node_ids.size(), 1);
        int ridx = 0;
        for (auto iter = effe_marg_node_ids.begin(); iter != effe_marg_node_ids.end(); ++iter) {
            ordering_(ridx++, 0) = *iter;
        }
        for (auto iter = effe_nomg_node_ids.begin(); iter != effe_nomg_node_ids.end(); ++iter) {
            ordering_(ridx++, 0) = *iter;
        }
    }
    virtual void BuildProblem(SparseBlockMatrix &A, MatrixSparseArray &b) override final {
        if (!HasMarginalization()) {
            SymbolicAnalysis();
        } else {
            MarginalizationAnalysis();
        }
        A = SparseBlockMatrix(ordering_.Rows(), ordering_.Cols());
        b = MatrixSparseArray(ordering_.Rows());
        std::map<int, int> id2ordering;
        for (int i = 0; i < ordering_.Rows(); ++i) {
            id2ordering[ordering_(i, 0)] = i;
        }
        for (int i = 0; i < factors_.size(); ++i) {
            for (int j = 0; j < factors_[i]->NumVertices(); ++j) {
                int vid1 = factors_[i]->GetVertexIdAt(j);
                if (!id2ordering.count(vid1)) {
                    continue;
                }
                for (int k = 0; k < factors_[i]->NumVertices(); ++k) {
                    int vid2 = factors_[i]->GetVertexIdAt(k);
                    if (!id2ordering.count(vid2)) {
                        continue;
                    }
                    if (A.CountIndex(id2ordering[vid1], id2ordering[vid2])) {
                        A(id2ordering[vid1], id2ordering[vid2]) += factors_[i]->JacobTInfoJacob(j, k);
                    } else {
                        A(id2ordering[vid1], id2ordering[vid2]) = factors_[i]->JacobTInfoJacob(j, k);
                    }
                }
                if (b.CountIndex(id2ordering[j])) {
                    b[id2ordering[j]] += factors_[i]->JacobTInfoRes(j);
                } else {
                    b[id2ordering[j]] = factors_[i]->JacobTInfoRes(j);
                }
            }
        }
    }
    virtual bool UpdateInc(const MatrixSparseArray &inc) override final {
        double chi_before_inc = 0;
        for (int i = 0; i < factors_.size(); ++i) {
            chi_before_inc += factors_[i]->GetChi();
        }
        for (int i = 0; i < ordering_.Rows(); ++i) {
            id2v_[ordering_(i, 0)]->Push();
            id2v_[ordering_(i, 0)]->ApplyPerturbation(inc[i].Data());
        }
        double chi_after_inc = 0;
        for (int i = 0; i < factors_.size(); ++i) {
            chi_after_inc += factors_[i]->GetChi();
        }
        if (chi_after_inc < chi_before_inc) {
            return true;
        }
        for (int i = 0; i < ordering_.Rows(); ++i) {
            id2v_[ordering_(i, 0)]->Pop();
        }
        return false;
    }
    virtual double ComputeGraphChi() const override final {
        double chi = 0;
        for (int i = 0; i < factors_.size(); ++i) {
            chi += factors_[i]->GetChi();
        }
        return chi;
    }

    virtual void SymbolicAnalysis() = 0;

protected:
    SparseMatrixi ordering_;
};

} // namespace hitnlls