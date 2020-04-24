#pragma once

#include <memory>

#include "matrix/dense.h"
#include "admm/variable.h"

namespace hitcadmm {

class Equality {
public:
    typedef std::shared_ptr<Equality> Ptr;

    explicit Equality(const Variable::Ptr &v, const hitnlls::matrix::MatrixXd &A, const hitnlls::matrix::VectorXd &b) { v_ = v; A_ = A; b_ = b; }
    
    Variable::Ptr GetV() const { return v_; }
    hitnlls::matrix::MatrixXd GetA() const { return A_; }
    hitnlls::matrix::VectorXd Getb() const { return b_; }

private:
    Variable::Ptr v_;
    hitnlls::matrix::MatrixXd A_;
    hitnlls::matrix::VectorXd b_;
};

} // namespace hitcadmm