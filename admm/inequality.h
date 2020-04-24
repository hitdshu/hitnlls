#pragma once

#include <memory>

#include "matrix/dense.h"
#include "admm/variable.h"

namespace hitcadmm {

class Inequality {
public:
    typedef std::shared_ptr<Inequality> Ptr;

    explicit Inequality(const Variable::Ptr &v, const hitnlls::matrix::MatrixXd &C, const hitnlls::matrix::VectorXd &d, const Variable::Ptr &c) { v_ = v; C_ = C; d_ = d; c_ = c; }
    
    Variable::Ptr GetV() const { return v_; }
    hitnlls::matrix::MatrixXd GetC() const { return C_; }
    hitnlls::matrix::VectorXd Getd() const { return d_; }
    
    Variable::Ptr GetCone() const { return c_; }

private:
    Variable::Ptr v_;
    Variable::Ptr c_;
    hitnlls::matrix::MatrixXd C_;
    hitnlls::matrix::VectorXd d_;
};

} // namespace hitcadmm