#pragma once

#include "variable.h"

namespace nlls {

class Loss {
public:
    explicit Loss(Variable *v, const MatrixXf &P, const VectorXf &q) { v_ = v; P_ = P; q_ = q; }
    
    Variable *GetV() const { return v_; }
    MatrixXf GetP() const { return P_; }
    VectorXf Getq() const { return q_; }

private:
    Variable *v_;
    MatrixXf P_;
    VectorXf q_;
};

} // namespace nlls