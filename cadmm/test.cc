#include <Eigen/Dense>
#include <iostream>

#include "vertex/nnc.h"
#include "vertex/soc.h"
#include "vertex/sdc.h"
#include "vertex/var.h"
#include "solver/problem.h"
#include "solver/admm.h"
#include "factor/equality.h"
#include "factor/loss.h"
#include "factor/inequality.h"

using namespace hitcadmm;

int main(int argc, char **argv) {
    Eigen::Matrix2d P;
    P << 1, 0.5, 
        0.5, 1;
    Eigen::Vector2d q;
    q << 1, 1;
    double pho = 200;
    Eigen::Matrix<double, 1, 2> A;
    A << 1, -0.5;
    double b = -5;
    Eigen::Matrix<double, 1, 2> C;
    C << -1, 1;
    double d = 2;

    Problem::Ptr problem(new Problem());
    Var<2>::Ptr var(new Var<2>());
    Loss::Ptr loss(new Loss(var, P, q));
    Eigen::Matrix<double, 1, 1> b_v;
    b_v << b;
    Equality::Ptr eq(new Equality(var, A, b_v));
    Eigen::Matrix<double, 1, 1> d_v;
    d_v << d;
    Nnc<1>::Ptr iec(new Nnc<1>());
    Inequality::Ptr ie(new Inequality(var, C, d_v, iec));
    problem->AddLoss(loss);
    problem->AddEquality(eq);
    problem->AddInequality(ie);
    Admm::Ptr admm(new Admm());
    admm->SetProblem(problem);
    admm->Solve();

    std::cout << "Admm solved var: " << var->GetVector().transpose() << std::endl;
    
    return 0;
}
