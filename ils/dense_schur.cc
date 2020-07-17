#include "dense_schur.h"

namespace nlls {
namespace internal {

void DenseSchur::BuildStructure() {
    vr2i_.clear();
    vs2c_.clear();
    vs2vr2c_.clear();
    std::unordered_set<VertexBase *> &vs = Vertices();
    std::unordered_set<VertexBase *> &vss = SchurVertices();
    int dim_offset_reduce = 0;
    int dim_offset_schur = 0;
    for (auto &iter : vs) {
        if (!vss.count(iter)) {
            vr2i_[iter] = dim_offset_reduce;
            dim_offset_reduce += iter->Dof();
        } else {
            vs2i_[iter] = dim_offset_schur;
            dim_offset_schur += iter->Dof();
        }
    }
    for (auto iter : vss) {
        vs2c_[iter] = MatrixXd(iter->Dof(), iter->Dof());
    }
    Ar_.Resize(dim_offset_reduce, dim_offset_reduce);
    br_.Resize(dim_offset_reduce);
    bs_.Resize(dim_offset_schur);
}

void DenseSchur::FillStructure() {
    for (auto &iter : vs2c_) {
        iter.second.SetZero();
    }
    vs2vr2c_.clear();
    Ar_.SetZero();
    br_.SetZero();
    bs_.SetZero();
    std::unordered_set<FactorBase *> &fs = Factors();
    std::unordered_set<VertexBase *> &vss = SchurVertices();
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
            int p_off1 = factor->GetVertexDimOffset(i);
            for (int j = 0; j < factor->NumVertices(); ++j) {
                VertexBase *vtmp2 = factor->GetVertexAt(j);
                if (vtmp2->GetStatus() != VertexBase::Active) {
                    continue;
                }
                int v_dof2 = vtmp2->Dof();
                int p_off2 = factor->GetVertexDimOffset(j);
                if (vss.count(vtmp1) && vss.count(vtmp2)) {
                    vs2c_[vtmp1] += jacob_t_jacob.Block(p_off1, p_off2, v_dof1, v_dof2);
                } else if ((!vss.count(vtmp1)) && (!vss.count(vtmp2))) {
                    int v_off1 = vr2i_[vtmp1];
                    int v_off2 = vr2i_[vtmp2];
                    if (v_off1 >= v_off2) {
                        Ar_.Block(v_off1, v_off2, v_dof1, v_dof2) += jacob_t_jacob.Block(p_off1, p_off2, v_dof1, v_dof2);
                    }
                } else if (!vss.count(vtmp1) && vss.count(vtmp2)) {
                    if (vs2vr2c_[vtmp2][vtmp1].Size() == 0) {
                        vs2vr2c_[vtmp2][vtmp1] = jacob_t_jacob.Block(p_off1, p_off2, v_dof1, v_dof2);
                    } else {
                        vs2vr2c_[vtmp2][vtmp1] += jacob_t_jacob.Block(p_off1, p_off2, v_dof1, v_dof2);
                    }
                }
            }
            if (vss.count(vtmp1)) {
                int v_off1 = vs2i_[vtmp1];
                bs_.Block(v_off1, 0, v_dof1, 1) += jacob_t_resid.Block(p_off1, 0, v_dof1, 1);
            } else {
                int v_off1 = vr2i_[vtmp1];
                br_.Block(v_off1, 0, v_dof1, 1) += jacob_t_resid.Block(p_off1, 0, v_dof1, 1);
            }
        }
    }
}

float DenseSchur::InitLevenMarqParam() {
    using std::max;
    float lm = 0;
    for (int i = 0; i < Ar_.Rows(); ++i) {
        lm = max(lm, Ar_(i, i));
    }
    for (auto &iter : vs2c_) {
        for (int i = 0; i < iter.second.Rows(); ++i) {
            lm = max(lm, iter.second(i, i));
        }
    }
    lm *= 1e-3;
    return lm;
}

VectorXf DenseSchur::ComputeGaussNewtonPoint(float l) {
    if (l > 0.0) {
       for (int i = 0; i < Ar_.Rows(); ++i) {
            Ar_(i, i) += l;
        }
        for (auto &iter : vs2c_) {
            for (int i = 0; i < iter.second.Rows(); ++i) {
                iter.second(i, i) += l;
            }
        } 
    }
    MatrixXf Ar = Ar_;
    VectorXf br = br_;
    for (auto &iter : vs2c_) {
        VertexBase *sv = iter.first;
        int sv_dof = sv->Dof();
        int sv_off = vs2i_[sv];
        MatrixXf &sc = iter.second;
        MatrixXf sc_inv = sc.Inverse();
        std::unordered_map<VertexBase *, MatrixXf> &rscm = vs2vr2c_[sv];
        for (auto &rsc1 : rscm) {
            VertexBase *vtmp1 = rsc1.first;
            MatrixXf &c1 = rsc1.second;
            int v_dof1 = vtmp1->Dof();
            int v_off1 = vr2i_[vtmp1];
            for (auto &rsc2 : rscm) {
                int v_off2 = vr2i_[vtmp1];
                if (v_off1 >= v_off2) {
                    VertexBase *vtmp2 = rsc2.first;
                    MatrixXf &c2 = rsc2.second;
                    int v_dof2 = vtmp2->Dof();
                    Ar.Block(v_off1, v_off2, v_dof1, v_dof2) -= c1 * sc_inv * c2;
                }
            }
            br.Block(v_off1, 0, v_dof1, 1) -= c1 * sc_inv * bs_.Block(sv_off, 0, sv_dof, 1);
        }
    }
    LLT<MatrixXf> llt_r(Ar);
    VectorXf reduce_inc = llt_r.Solve(br);
    VectorXf update_bs = bs_;
    for (auto &iter : vs2c_) {
        VertexBase *sv = iter.first;
        int sv_dof = sv->Dof();
        int sv_off = vs2i_[sv];
        std::unordered_map<VertexBase *, MatrixXf> &rscm = vs2vr2c_[sv];
        for (auto &rsc : rscm) {
            VertexBase *vtmp = rsc.first;
            MatrixXf &c = rsc.second;
            int v_dof = vtmp->Dof();
            int v_off = vr2i_[vtmp];
            update_bs.Block(sv_off, 0, sv_dof, 1) -= c.Transpose() * reduce_inc.Block(v_off, 0, v_dof, 1);
        }
    }
    for (auto &iter : vs2c_) {
        VertexBase *sv = iter.first;
        MatrixXf &sc = iter.second;
        int sv_dof = sv->Dof();
        int sv_off = vs2i_[sv];
        update_bs.Block(sv_off, 0, sv_dof, 1) = sc.Inverse() * update_bs.Block(sv_off, 0, sv_dof, 1);
    }
    VectorXf result(reduce_inc.Size() + update_bs.Size());
    result.Block(0, 0, update_bs.Size(), 1) = update_bs;
    result.Block(update_bs.Size(), 0, reduce_inc.Size(), 1) = reduce_inc;
    return result;
}

VectorXf DenseSchur::ComputeCauchyPoint() {
    VectorXf g(bs_.Size() + br_.Size());
    g.Block(0, 0, bs_.Size(), 1) = -bs_;
    g.Block(bs_.Size(), 0, br_.Size(), 1) = -br_;
    float s = - pow(g.Norm(), 2) / WeightedNormSquared(g);
    return s * g;
}

void DenseSchur::Update(const VectorXf &inc) {
    std::unordered_set<VertexBase *> &vs = Vertices();
    std::unordered_set<VertexBase *> &vss = SchurVertices();
    int reduce_offset = bs_.Size();
    for (auto v : vs) {
        if (vss.count(v)) {
            v->Update(inc.Data() + vs2i_[v]);
        } else {
            v->Update(inc.Data() + vr2i_[v] + reduce_offset);
        }
    }
}

float DenseSchur::SurrogateReduction(const VectorXf &inc) {
    VectorXf b(br_.Size() + bs_.Size());
    b.Block(0, 0, bs_.Size(), 1) = bs_;
    b.Block(bs_.Size(), 0, br_.Size(), 1) = br_;
    float cost_reduce = b.Dot(inc) - 0.5 * WeightedNormSquared(inc);
    return cost_reduce;
}

float DenseSchur::WeightedNormSquared(const VectorXf &inc) {
    VectorXf inc_bs = inc.Block(0, 0, bs_.Size(), 1);
    VectorXf inc_br = inc.Block(bs_.Size(), 0, br_.Size(), 1);
    float wns = (inc_br.Transpose() * Ar_ * inc_br).Norm() * 2;
    for (int i = 0; i < Ar_.Rows(); ++i) {
        wns -= inc_br[i] * Ar_(i, i) * inc_br[i];
    }
    for (auto &iter : vs2c_) {
        VertexBase *sv = iter.first;
        MatrixXf &sc = iter.second;
        int sv_dof = sv->Dof();
        int sv_off = vs2i_[sv];
        VectorXf tmp = inc_bs.Block(sv_off, 0, sv_dof, 1);
        wns += (tmp.Transpose() * sc * tmp).Norm();
    }
    for (auto &iter : vs2c_) {
        VertexBase *sv = iter.first;
        int sv_dof = sv->Dof();
        int sv_off = vs2i_[sv];
        std::unordered_map<VertexBase *, MatrixXf> &rscm = vs2vr2c_[sv];
        for (auto &rsc : rscm) {
            VertexBase *vtmp = rsc.first;
            MatrixXf &c = rsc.second;
            int v_dof = vtmp->Dof();
            int v_off = vr2i_[vtmp];
            wns += (inc_br.Block(v_off, 0, v_dof, 1).Transpose() * c * inc_bs.Block(sv_off, 0, sv_dof, 1))(0, 0) * 2;
        }
    }
    return wns;
}

} // namespace internal
} // namespace nlls