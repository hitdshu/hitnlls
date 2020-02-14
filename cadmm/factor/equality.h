#pragma once

#include <memory>
#include <Eigen/Dense>
#include "vertex/vertex.h"

namespace hitcadmm {

class Equality {
public:
    typedef std::shared_ptr<Equality> Ptr;

    explicit Equality(const Vertex::Ptr &v, const Eigen::MatrixXd &A, const Eigen::VectorXd &b) { v_ = v; A_ = A; b_ = b; }
    
    Vertex::Ptr GetV() const { return v_; }
    Eigen::MatrixXd GetA() const { return A_; }
    Eigen::VectorXd Getb() const { return b_; }

private:
    Vertex::Ptr v_;
    Eigen::MatrixXd A_;
    Eigen::VectorXd b_;
};

} // namespace hitcadmm