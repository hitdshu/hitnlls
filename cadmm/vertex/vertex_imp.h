#pragma once

#include "vertex/vertex.h"

namespace hitcadmm {

template <int ndim, typename ValueType>
class VertexImp : public Vertex {
public:
    explicit VertexImp() : Vertex(ndim) {}

    virtual void Project() = 0;
    virtual void SetVector(const Eigen::VectorXd &v) = 0;
    virtual Eigen::VectorXd GetVector() const = 0;

    virtual void SetValue(const ValueType &val) final { val_ = val; }
    virtual ValueType GetValue() const final { return val_; }

protected:
    bool CheckDim(const Eigen::VectorXd &v) const {
        if (v.size() < ndim) {
            return false;
        } else {
            return true;
        }
    }

private:
    ValueType val_;
};

} // namespace hitcadmm