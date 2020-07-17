#include <map>
#include "sparse_cholesky.h"

namespace nlls {
namespace internal {

void SparseCholesky::BuildStructure() {
    orders_.clear();
    v2i_.clear();
    vc2vr2c_.clear();
    ComputeOrdering();
    int dim_offset = 0;
    for (auto v : orders_) {
        v2i_[v] = dim_offset;
        dim_offset += v->Dof();
    }
    b_.Resize(dim_offset);
}

void SparseCholesky::FillStructure() {
    vc2vr2c_.clear();
    b_.SetZero();
    std::unordered_set<FactorBase *> &fs = Factors();
    for (auto &factor : fs) {
        factor->Compute();
        MatrixXf total_jacob = factor->Jacobian();
        VectorXf total_resid = factor->Residual();
        MatrixXf jacob_t_jacob = total_jacob.Transpose() * total_jacob;
        VectorXf jacob_t_resid = total_jacob.Transpose() * total_resid;
        for (int i = 0; i < factor->NumVertices(); ++i) {
            VertexBase *vtmp1 = factor->GetVertexAt(i);
            if (vtmp1->GetStatus() != VertexBase::Active) {
                continue;
            }
            int v_dof1 = vtmp1->Dof();
            int v_off1 = v2i_[vtmp1];
            int p_off1 = factor->GetVertexDimOffset(i);
            for (int j = 0; j < factor->NumVertices(); ++j) {
                VertexBase *vtmp2 = factor->GetVertexAt(j);
                if (vtmp2->GetStatus() != VertexBase::Active) {
                    continue;
                }
                int v_dof2 = vtmp2->Dof();
                int v_off2 = v2i_[vtmp2];
                if (v_off1 >= v_off2) {
                    int p_off2 = factor->GetVertexDimOffset(j);
                    if (vc2vr2c_[vtmp2][vtmp1].Size() == 0) {
                        vc2vr2c_[vtmp2][vtmp1] = jacob_t_jacob.Block(p_off1, p_off2, v_dof1, v_dof2);
                    } else {
                        vc2vr2c_[vtmp2][vtmp1] += jacob_t_jacob.Block(p_off1, p_off2, v_dof1, v_dof2);
                    }
                }
            }
            b_.Block(v_off1, 0, v_dof1, 1) += jacob_t_resid.Block(p_off1, 0, v_dof1, 1);
        }
    }
}

float SparseCholesky::InitLevenMarqParam() {
    using std::max;
    float lm = 0;
    for (auto v : orders_) {
        MatrixXf &c = vc2vr2c_[v][v];
        for (int i = 0; i < c.Rows(); ++i) {
            lm = max(lm, c(i, i));
        }
    }
    lm *= 1e-3;
    return lm;
}

VectorXf SparseCholesky::ComputeGaussNewtonPoint(float l) {
    if (l > 0.0) {
        for (auto v : orders_) {
            MatrixXf &c = vc2vr2c_[v][v];
            for (int i = 0; i < c.Rows(); ++i) {
                c(i, i) += l;
            }
        }
    }
    std::unordered_map<VertexBase *, std::unordered_map<VertexBase *, MatrixXf>> choll = vc2vr2c_;
    for (auto &v : orders_) {
        int v_dof = v->Dof();
        int v_off = v2i_[v];
        MatrixXf &c = choll[v][v];
        LLT<MatrixXf> cllt(c);
        MatrixXf cl = cllt.L();
        MatrixXf clti = cl.Transpose().Inverse();
        MatrixXf c_inv = c.Inverse();
        std::unordered_map<VertexBase *, MatrixXf> &rscm = choll[v];
        for (auto &rsc1 : rscm) {
            VertexBase *vtmp1 = rsc1.first;
            MatrixXf &c1 = rsc1.second;
            int v_off1 = v2i_[vtmp1];
            for (auto &rsc2 : rscm) {
                if (rsc2.first == v) {
                    continue;
                }
                VertexBase *vtmp2 = rsc2.first;
                int v_off2 = v2i_[vtmp2];
                if (v_off1 >= v_off2) {
                    MatrixXf &c2 = rsc2.second;
                    if (choll[vtmp2][vtmp1].Size() == 0) {
                        choll[vtmp2][vtmp1] = -c1 * c_inv * c2;
                    } else {
                        choll[vtmp2][vtmp1] -= c1 * c_inv * c2;
                    }
                }
            }
        }
        rscm[v] = cl;
        for (auto &rsc : rscm) {
            if (rsc.first != v) {
                rsc.second *= clti;
            }
        }
    }
    VectorXf result(b_);
    for (auto v : orders_) {
        int v_off = v2i_[v];
        int v_dof = v->Dof();
        std::unordered_map<VertexBase *, MatrixXf> &rslm = choll[v];
        MatrixXf &diag_l = rslm[v];
        MatrixXf diag_linv = diag_l.Inverse();
        VectorXf tmp_result = diag_linv * result.Block(v_off, 0, v_dof, 1);
        result.Block(v_off, 0, v_dof, 1) = tmp_result;
        for (auto &iter : rslm) {
            if (iter.first != v) {
                int tmp_off = v2i_[iter.first];
                int tmp_dof = iter.first->Dof();
                result.Block(tmp_off, 0, tmp_dof, 1) -= iter.second * tmp_result;
            }
        }
    }
    for (int i = orders_.size() - 1; i >= 0; --i) {
        auto v = orders_[i];
        int v_off = v2i_[v];
        int v_dof = v->Dof();
        std::unordered_map<VertexBase *, MatrixXf> &rslm = choll[v];
        for (auto &iter : rslm) {
            if (iter.first != v) {
                int tmp_off = v2i_[iter.first];
                int tmp_dof = iter.first->Dof();
                result.Block(v_off, 0, v_dof, 1) -= iter.second.Transpose() * result.Block(tmp_off, 0, tmp_dof, 1);
            }
        }
        result.Block(v_off, 0, v_dof, 1) = rslm[v].Transpose().Inverse() * result.Block(v_off, 0, v_dof, 1);
    }
    return result;
}

VectorXf SparseCholesky::ComputeCauchyPoint() {
    VectorXf g(-b_);
    float s = - pow(g.Norm(), 2) / WeightedNormSquared(g);
    return s * g;
}

void SparseCholesky::Update(const VectorXf &inc) {
    for (auto v : orders_) {
        v->Update(inc.Data() + v2i_[v]);
    }
}

float SparseCholesky::SurrogateReduction(const VectorXf &inc) {
    float cost_reduce = b_.Dot(inc) - 0.5 * WeightedNormSquared(inc);
    return cost_reduce;
}

void SparseCholesky::ComputeOrdering() {
    auto vs = Vertices();
    auto schur_vs = SchurVertices();
    if (schur_vs.size() > 0) {
        for (auto v : schur_vs) {
            orders_.push_back(v);
        }
        for (auto v : vs) {
            if (!schur_vs.count(v)) {
                orders_.push_back(v);
            }
        }
    } else {
        std::unordered_map<VertexBase *, std::unordered_map<VertexBase *, int>> v2v2c;
        std::unordered_set<FactorBase *> &fs = Factors();
        for (auto &factor : fs) {
            for (int i = 0; i < factor->NumVertices(); ++i) {
                VertexBase *vi = factor->GetVertexAt(i);
                for (int j = 0; j < factor->NumVertices(); ++j) {
                    VertexBase *vj = factor->GetVertexAt(j);
                    v2v2c[vi][vj] = 1;
                    v2v2c[vj][vi] = 1;
                }
            }
        }
        std::multimap<int, VertexBase *> c2v;
        for (auto tmp : v2v2c) {
            if (vs.count(tmp.first) && tmp.first->GetStatus() == VertexBase::Active) {
                c2v.insert(std::make_pair(int(tmp.second.size()), tmp.first));
            }
        }
        for (auto item : c2v) {
            orders_.push_back(item.second);
        }
    }
}

float SparseCholesky::WeightedNormSquared(const VectorXf &inc) {
    float wns = 0;
    for (auto &iter : vc2vr2c_) {
        VertexBase *sv = iter.first;
        int sv_dof = sv->Dof();
        int sv_off = v2i_[sv];
        std::unordered_map<VertexBase *, MatrixXf> &rscm = vc2vr2c_[sv];
        for (auto &rsc : rscm) {
            VertexBase *vtmp = rsc.first;
            MatrixXf &c = rsc.second;
            int v_dof = vtmp->Dof();
            int v_off = v2i_[vtmp];
            float tmp_norm = (inc.Block(v_off, 0, v_dof, 1).Transpose() * c * inc.Block(sv_off, 0, sv_dof, 1))(0, 0);
            if (vtmp == sv) {
                wns += tmp_norm;
            } else {
                wns += tmp_norm * 2;
            }
        }
    }
    return wns;
}

} // namespace internal
} // namespace nlls