#pragma once

#include "variable.h"

namespace nlls {

class Inequality {
public:
    explicit Inequality(Variable *v, const MatrixXf &C, const VectorXf &d, Variable *c) { v_ = v; C_ = C; d_ = d; c_ = c; }
    
    Variable *GetV() const { return v_; }
    MatrixXf GetC() const { return C_; }
    VectorXf Getd() const { return d_; }
    
    Variable *GetCone() const { return c_; }

private:
    Variable *v_;
    Variable *c_;
    MatrixXf C_;
    VectorXf d_;
};

} // namespace nlls