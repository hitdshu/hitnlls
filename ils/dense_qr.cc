#include <algorithm>

#include "dense_qr.h"

namespace nlls {
namespace internal {
void DenseQr::BuildStructure() {
    v2i_.clear();
    f2i_.clear();
    std::unordered_set<FactorBase *> &fs = Factors();
    int fac_offset = 0;
    for (auto &iter : fs) {
        f2i_[iter] = fac_offset;
        fac_offset += iter->NumObservations();
    }
    std::unordered_set<VertexBase *> &vs = Vertices();
    int dim_offset = 0;
    for (auto &iter : vs) {
        v2i_[iter] = dim_offset;
        dim_offset += iter->Dof();
    }
    if (s_ == LEVEN_MARQ) {
        fac_offset += dim_offset;
    }
    A_.Resize(fac_offset, dim_offset);
    b_.Resize(fac_offset);
}

void DenseQr::FillStructure() {
    A_.SetZero();
    b_.SetZero();
    std::unordered_set<FactorBase *> &fs = Factors();
    for (auto &factor : fs) {
        int fac_off = f2i_[factor];
        int fac_nobs = factor->NumObservations();
        factor->Compute();
        for (int i = 0; i < factor->NumVertices(); ++i) {
            VertexBase *vtmp = factor->GetVertexAt(i);
            if (vtmp->GetStatus() != VertexBase::Active) {
                continue;
            }
            int v_off = v2i_[vtmp];
            int v_dof = vtmp->Dof();
            A_.Block(fac_off, v_off, fac_nobs, v_dof) += factor->Jacobian(i);
        }
        b_.Block(fac_off, 0, fac_nobs, 1) += factor->Residual();
    }
}

float DenseQr::InitLevenMarqParam() {
    using std::max;
    float lm = 0;
    for (int j = 0; j < A_.Cols(); ++j) {
        float tmp = 0;
        for (int i = 0; i < A_.Rows(); ++i) {
            tmp += pow(A_(i, j), 2);
        }
        lm = max(lm, tmp);
    }
    lm *= 1e-3;
    return lm;
}

VectorXf DenseQr::ComputeGaussNewtonPoint(float l) {
    if (l > 0) {
        int fac_off = A_.Rows();
        int dim_off = A_.Cols();
        A_.Block(fac_off - dim_off, 0, dim_off, dim_off) = MatrixXd::Identity(dim_off, dim_off) * l;
    }
    QR<MatrixXf> qr(A_);
    return qr.Solve(b_);
}

VectorXf DenseQr::ComputeCauchyPoint() {
    VectorXf g = - A_.Transpose() * b_;
    float s = - pow(g.Norm(), 2) / pow((A_ * g).Norm(), 2);
    return s * g;
}

void DenseQr::Update(const VectorXf &inc) {
    std::unordered_set<VertexBase *> &vs = Vertices();
    for (auto v : vs) {
        v->Update(inc.Data() + v2i_[v]);
    }
}

float DenseQr::SurrogateReduction(const VectorXf &inc) {
    float cost_bef = 0.5 * pow(b_.Norm(), 2);
    float cost_aft = 0.5 * pow((A_ * inc - b_).Norm(), 2);
    float cost_reduce = cost_bef - cost_aft;
    return cost_reduce;
}

} // namespace internal
} // namespace nlls