#pragma once

#include <vector>

#include "loss.h"
#include "equality.h"
#include "inequality.h"

namespace nlls {
namespace internal {
class Admm;
} // namespace internal
class Question {
public:
    ~Question();
    
    void AddLoss(Loss *loss) { losses_.push_back(loss); }
    void AddEquality(Equality *equa) { equas_.push_back(equa); }
    void AddInequality(Inequality *ine) { ines_.push_back(ine); }

    void Solve();

private:
    friend class internal::Admm;

    void Build();
    void SetEstimates(const VectorXf &est);

    MatrixXf GetP() const { return P_; }
    VectorXf Getq() const { return q_; }
    MatrixXf GetA() const { return A_; }
    VectorXf Getb() const { return b_; }
    MatrixXf GetC() const { return C_; }
    VectorXf Getd() const { return d_; }
    int GetXDim() const { return xdim_; }
    int GetEDim() const { return edim_; }
    int GetSDim() const { return sdim_; }
    std::vector<Variable *> GetCones() const { return cones_; }

    std::vector<Loss *> losses_;
    std::vector<Equality *> equas_;
    std::vector<Inequality *> ines_;
    std::vector<Variable *> xs_;
    MatrixXf P_;
    VectorXf q_;
    MatrixXf A_;
    VectorXf b_;
    MatrixXf C_;
    VectorXf d_;
    int xdim_;
    int edim_;
    int sdim_;
    std::vector<Variable *> cones_;
};

} // namespace nlls