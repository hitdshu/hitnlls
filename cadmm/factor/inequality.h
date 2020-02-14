#pragma once

#include <memory>
#include <Eigen/Dense>
#include "vertex/vertex.h"

namespace hitcadmm {

class Inequality {
public:
    typedef std::shared_ptr<Inequality> Ptr;

    explicit Inequality(const Vertex::Ptr &v, const Eigen::MatrixXd &C, const Eigen::VectorXd &d, const Vertex::Ptr &c) { v_ = v; C_ = C; d_ = d; c_ = c; }
    
    Vertex::Ptr GetV() const { return v_; }
    Eigen::MatrixXd GetC() const { return C_; }
    Eigen::VectorXd Getd() const { return d_; }
    
    Vertex::Ptr GetCone() const { return c_; }

private:
    Vertex::Ptr v_;
    Vertex::Ptr c_;
    Eigen::MatrixXd C_;
    Eigen::VectorXd d_;
};

} // namespace hitcadmm