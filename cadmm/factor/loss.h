#pragma once

#include <memory>
#include <Eigen/Dense>
#include "vertex/vertex.h"

namespace hitcadmm {

class Loss {
public:
    typedef std::shared_ptr<Loss> Ptr;

    explicit Loss(const Vertex::Ptr &v, const Eigen::MatrixXd &P, const Eigen::VectorXd &q) { v_ = v; P_ = P; q_ = q; }
    
    Vertex::Ptr GetV() const { return v_; }
    Eigen::MatrixXd GetP() const { return P_; }
    Eigen::VectorXd Getq() const { return q_; }

private:
    Vertex::Ptr v_;
    Eigen::MatrixXd P_;
    Eigen::VectorXd q_;
};

} // namespace hitcadmm