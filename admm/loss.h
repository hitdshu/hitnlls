#pragma once

#include <memory>

#include "matrix/dense.h"
#include "admm/variable.h"

namespace hitcadmm {

class Loss {
public:
    typedef std::shared_ptr<Loss> Ptr;

    explicit Loss(const Variable::Ptr &v, const hitnlls::matrix::MatrixXd &P, const hitnlls::matrix::VectorXd &q) { v_ = v; P_ = P; q_ = q; }
    
    Variable::Ptr GetV() const { return v_; }
    hitnlls::matrix::MatrixXd GetP() const { return P_; }
    hitnlls::matrix::VectorXd Getq() const { return q_; }

private:
    Variable::Ptr v_;
    hitnlls::matrix::MatrixXd P_;
    hitnlls::matrix::VectorXd q_;
};

} // namespace hitcadmm