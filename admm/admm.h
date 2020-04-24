#pragma once

#include "admm/problem.h"

namespace hitcadmm {

class Admm {
public:
    typedef std::shared_ptr<Admm> Ptr;

    explicit Admm(int iter_num = 20) { iter_num_ = iter_num; }

    void SetProblem(const Problem::Ptr &p) { p_ = p; }
    void Solve();

private:
    Problem::Ptr p_;
    int iter_num_;
};

void Admm::Solve() {
    p_->Build();

    hitnlls::matrix::MatrixXd P = p_->GetP();
    hitnlls::matrix::VectorXd q = p_->Getq();
    hitnlls::matrix::MatrixXd A = p_->GetA();
    hitnlls::matrix::VectorXd b = p_->Getb();
    hitnlls::matrix::MatrixXd C = p_->GetC();
    hitnlls::matrix::VectorXd d = p_->Getd();
    int xdim = p_->GetXDim();
    int edim = p_->GetEDim();
    int sdim = p_->GetSDim();
    std::vector<Variable::Ptr> cones = p_->GetCones();
    hitnlls::matrix::VectorXd x(xdim);
    x.SetZero();
    hitnlls::matrix::VectorXd s(sdim);
    s.SetZero();
    hitnlls::matrix::VectorXd s_hat(sdim);
    s_hat.SetZero();
    hitnlls::matrix::VectorXd u(sdim);
    u.SetZero();

    double pho = 2;
    for (int idx = 0; idx < iter_num_; ++idx) {
        hitnlls::matrix::VectorXd res(xdim + sdim + edim + sdim);
        res.Block(0, 0, xdim, 1) = q;
        res.Block(xdim, 0, sdim, 1) = pho * (u - s_hat);
        res.Block(xdim + sdim, 0, edim, 1) = -b;
        res.Block(xdim + sdim + edim, 0, sdim, 1) = -d;
        res = -res;
        hitnlls::matrix::MatrixXd coeff(xdim + sdim + edim + sdim, xdim + sdim + edim + sdim);
        coeff.SetZero();
        coeff.Block(0, 0, xdim, xdim) = P;
        coeff.Block(0, xdim + sdim, xdim, edim) = A.Transpose();
        coeff.Block(0, xdim + sdim + edim, xdim, sdim) = C.Transpose();
        coeff.Block(xdim, xdim, sdim, sdim) = pho * hitnlls::matrix::MatrixXd::Identity(sdim, sdim);
        coeff.Block(xdim, xdim + sdim + edim, sdim, sdim) = -hitnlls::matrix::MatrixXd::Identity(sdim, sdim);
        coeff.Block(xdim + sdim, 0, edim, xdim) = A;
        coeff.Block(xdim + sdim + edim, 0, sdim, xdim) = C;
        coeff.Block(xdim + sdim + edim, xdim, sdim, sdim) = -hitnlls::matrix::MatrixXd::Identity(sdim, sdim);
        hitnlls::matrix::VectorXd sol = coeff.Inverse() * res;
        x = sol.Block(0, 0, xdim, 1);
        s = sol.Block(xdim, 0, sdim, 1);
        hitnlls::matrix::VectorXd sho = s_hat;
        s_hat = s + u;
        int shidx = 0;
        for (auto &cone : cones) {
            int cone_dim = cone->GetDim();
            cone->SetVector(s_hat.Block(shidx, 0, cone_dim, 1));
            cone->Project();
            s_hat.Block(shidx, 0, cone_dim, 1) = cone->GetVector();
        }
        u += s - s_hat;
        double pr = (s - s_hat).Norm();
        double dr = (s_hat - sho).Norm() * pho;
        if (pr > 10 * dr) {
            pho *= 4;
            u = u / 4;
        } else if (10 * pr < dr) {
            pho /= 4;
            u *= 4;
        }
    }

    p_->SetEstimates(x);
}

} // namespace hitcadmm