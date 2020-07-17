#pragma once

#include "variable.h"

namespace nlls {

class Equality {
public:
    explicit Equality(Variable *v, const MatrixXf &A, const VectorXf &b) { v_ = v; A_ = A; b_ = b; }
    
    Variable *GetV() const { return v_; }
    MatrixXf GetA() const { return A_; }
    VectorXf Getb() const { return b_; }

private:
    Variable *v_;
    MatrixXf A_;
    VectorXf b_;
};

} // namespace nlls