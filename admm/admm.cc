#include "admm.h"

namespace nlls {
namespace internal {

void Admm::Solve() {
    p_->Build();
    MatrixXf P = p_->GetP();
    VectorXf q = p_->Getq();
    MatrixXf A = p_->GetA();
    VectorXf b = p_->Getb();
    MatrixXf C = p_->GetC();
    VectorXf d = p_->Getd();
    int xdim = p_->GetXDim();
    int edim = p_->GetEDim();
    int sdim = p_->GetSDim();
    std::vector<Variable *> cones = p_->GetCones();
    VectorXf x(xdim);
    x.SetZero();
    VectorXf s(sdim);
    s.SetZero();
    VectorXf s_hat(sdim);
    s_hat.SetZero();
    VectorXf u(sdim);
    u.SetZero();
    float pho = 2;
    const float step = 4;
    const float ratio = 10;
    for (int idx = 0; idx < iter_num_; ++idx) {
        VectorXf res(xdim + sdim + edim + sdim);
        res.Block(0, 0, xdim, 1) = q;
        res.Block(xdim, 0, sdim, 1) = pho * (u - s_hat);
        res.Block(xdim + sdim, 0, edim, 1) = -b;
        res.Block(xdim + sdim + edim, 0, sdim, 1) = -d;
        res = -res;
        MatrixXf coeff(xdim + sdim + edim + sdim, xdim + sdim + edim + sdim);
        coeff.SetZero();
        coeff.Block(0, 0, xdim, xdim) = P;
        coeff.Block(0, xdim + sdim, xdim, edim) = A.Transpose();
        coeff.Block(0, xdim + sdim + edim, xdim, sdim) = C.Transpose();
        coeff.Block(xdim, xdim, sdim, sdim) = pho * MatrixXf::Identity(sdim, sdim);
        coeff.Block(xdim, xdim + sdim + edim, sdim, sdim) = -MatrixXf::Identity(sdim, sdim);
        coeff.Block(xdim + sdim, 0, edim, xdim) = A;
        coeff.Block(xdim + sdim + edim, 0, sdim, xdim) = C;
        coeff.Block(xdim + sdim + edim, xdim, sdim, sdim) = -MatrixXf::Identity(sdim, sdim);
        VectorXf sol = coeff.Inverse() * res;
        x = sol.Block(0, 0, xdim, 1);
        s = sol.Block(xdim, 0, sdim, 1);
        VectorXf sho = s_hat;
        s_hat = s + u;
        int shidx = 0;
        for (auto &cone : cones) {
            int cone_dim = cone->GetDim();
            cone->SetVector(s_hat.Block(shidx, 0, cone_dim, 1));
            cone->Project();
            s_hat.Block(shidx, 0, cone_dim, 1) = cone->GetVector();
        }
        u += s - s_hat;
        float pr = (s - s_hat).Norm();
        float dr = (s_hat - sho).Norm() * pho;
        if (pr > ratio * dr) {
            pho *= step;
            u = u / step;
        } else if (ratio * pr < dr) {
            pho /= step;
            u *= step;
        }
    }
    p_->SetEstimates(x);
}

} // namespace internal
} // namespace nlls