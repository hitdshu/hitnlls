#pragma once

#include <Eigen/Dense>
#include <memory>

namespace hitcadmm {

class Vertex {
public:
    typedef std::shared_ptr<Vertex> Ptr;

    explicit Vertex(int ndim) { ndim_ = ndim; }

    virtual void Project() = 0;
    virtual void SetVector(const Eigen::VectorXd &v) = 0;
    virtual Eigen::VectorXd GetVector() const = 0;

    virtual int GetDim() final { return ndim_; }

private:
    int ndim_;
};

} // namespace hitcadmm
