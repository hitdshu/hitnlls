#include <algorithm>

#include "dense_cholesky.h"

namespace nlls {
namespace internal {
void DenseCholesky::BuildStructure() {
    v2i_.clear();
    std::unordered_set<VertexBase *> &vs = Vertices();
    int dim_offset = 0;
    for (auto &iter : vs) {
        v2i_[iter] = dim_offset;
        dim_offset += iter->Dof();
    }
    A_.Resize(dim_offset, dim_offset);
    b_.Resize(dim_offset);
}

void DenseCholesky::FillStructure() {
    A_.SetZero();
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
            int v_off1 = v2i_[vtmp1];
            int v_dof1 = vtmp1->Dof();
            int p_off1 = factor->GetVertexDimOffset(i);
            for (int j = 0; j < factor->NumVertices(); ++j) {
                VertexBase *vtmp2 = factor->GetVertexAt(j);
                if (vtmp2->GetStatus() != VertexBase::Active) {
                    continue;
                }
                int v_off2 = v2i_[vtmp2];
                int v_dof2 = vtmp1->Dof();
                int p_off2 = factor->GetVertexDimOffset(j);
                if (v_off1 >= v_off2) {
                    A_.Block(v_off1, v_off2, v_dof1, v_dof2) += jacob_t_jacob.Block(p_off1, p_off2, v_dof1, v_dof2);
                }
            }
            b_.Block(v_off1, 0, v_dof1, 1) += jacob_t_resid.Block(p_off1, 0, v_dof1, 1);
        }
    }
     for (int i = 0; i < A_.Rows(); ++i) {
        for (int j = i + 1; j < A_.Cols(); ++j) {
            A_(i, j) = A_(j, i);
        }
    }
}

float DenseCholesky::InitLevenMarqParam() {
    using std::max;
    float lm = 0;
    for (int i = 0; i < A_.Rows(); ++i) {
        lm = max(lm, A_(i, i));
    }
    lm *= 1e-3;
    return lm;
}

VectorXf DenseCholesky::ComputeGaussNewtonPoint(float l) {
    if (l > 0) {
        for (int i = 0; i < A_.Rows(); ++i) {
            A_(i, i) += l;
        }
    }
    LLT<MatrixXf> llt(A_);
    return llt.Solve(b_);
}

VectorXf DenseCholesky::ComputeCauchyPoint() {
    VectorXf g = - b_;
    float s = - pow(g.Norm(), 2) / (g.Transpose() * A_ * g).Norm();
    return s * g;
}

void DenseCholesky::Update(const VectorXf &inc) {
    std::unordered_set<VertexBase *> &vs = Vertices();
    for (auto v : vs) {
        v->Update(inc.Data() + v2i_[v]);
    }
}

float DenseCholesky::SurrogateReduction(const VectorXf &inc) {
    float cost_reduce = b_.Dot(inc) - 0.5 * (inc.Transpose() * A_ * inc).Norm();
    return cost_reduce;
}
} // namespace internal
} // namespace nlls