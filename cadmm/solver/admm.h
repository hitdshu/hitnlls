#pragma once

#include "solver/problem.h"

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

    Eigen::MatrixXd P = p_->GetP();
    Eigen::VectorXd q = p_->Getq();
    Eigen::MatrixXd A = p_->GetA();
    Eigen::VectorXd b = p_->Getb();
    Eigen::MatrixXd C = p_->GetC();
    Eigen::VectorXd d = p_->Getd();
    int xdim = p_->GetXDim();
    int edim = p_->GetEDim();
    int sdim = p_->GetSDim();
    std::vector<Vertex::Ptr> cones = p_->GetCones();
    Eigen::VectorXd x(xdim);
    x.setZero();
    Eigen::VectorXd s(sdim);
    s.setZero();
    Eigen::VectorXd s_hat(sdim);
    s_hat.setZero();
    Eigen::VectorXd u(sdim);
    u.setZero();

    double pho = 2;
    for (int idx = 0; idx < iter_num_; ++idx) {
        Eigen::VectorXd res(xdim + sdim + edim + sdim);
        res.block(0, 0, xdim, 1) = q;
        res.block(xdim, 0, sdim, 1) = pho * (u - s_hat);
        res.block(xdim + sdim, 0, edim, 1) = -b;
        res.block(xdim + sdim + edim, 0, sdim, 1) = -d;
        res = -res;
        Eigen::MatrixXd coeff(xdim + sdim + edim + sdim, xdim + sdim + edim + sdim);
        coeff.setZero();
        coeff.block(0, 0, xdim, xdim) = P;
        coeff.block(0, xdim + sdim, xdim, edim) = A.transpose();
        coeff.block(0, xdim + sdim + edim, xdim, sdim) = C.transpose();
        coeff.block(xdim, xdim, sdim, sdim) = pho * Eigen::MatrixXd::Identity(sdim, sdim);
        coeff.block(xdim, xdim + sdim + edim, sdim, sdim) = -Eigen::MatrixXd::Identity(sdim, sdim);
        coeff.block(xdim + sdim, 0, edim, xdim) = A;
        coeff.block(xdim + sdim + edim, 0, sdim, xdim) = C;
        coeff.block(xdim + sdim + edim, xdim, sdim, sdim) = -Eigen::MatrixXd::Identity(sdim, sdim);
        Eigen::VectorXd sol = coeff.inverse() * res;
        x = sol.block(0, 0, xdim, 1);
        s = sol.block(xdim, 0, sdim, 1);
        Eigen::VectorXd sho = s_hat;
        s_hat = s + u;
        int shidx = 0;
        for (auto &cone : cones) {
            int cone_dim = cone->GetDim();
            cone->SetVector(s_hat.block(shidx, 0, cone_dim, 1));
            cone->Project();
            s_hat.block(shidx, 0, cone_dim, 1) = cone->GetVector();
        }
        u += s - s_hat;
        double pr = (s - s_hat).norm();
        double dr = (s_hat - sho).norm() * pho;
        if (pr > 10 * dr) {
            pho *= 4;
            u /= 4;
        } else if (10 * pr < dr) {
            pho /= 4;
            u *= 4;
        }
    }

    p_->SetEstimates(x);
}

} // namespace hitcadmm