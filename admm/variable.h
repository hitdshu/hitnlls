#pragma once

#include <memory>

#include "matrix/dense.h"

namespace hitcadmm {

class Variable {
public:
    typedef std::shared_ptr<Variable> Ptr;

    explicit Variable(int ndim) { ndim_ = ndim; }

    virtual void Project() = 0;
    virtual void SetVector(const hitnlls::matrix::VectorXd &v) = 0;
    virtual hitnlls::matrix::VectorXd GetVector() const = 0;

    int GetDim() { return ndim_; }

private:
    int ndim_;
};

template <int ndim, typename ValueType>
class VariableImp : public Variable {
public:
    explicit VariableImp() : Variable(ndim) {}

    virtual void Project() = 0;
    virtual void SetVector(const hitnlls::matrix::VectorXd &v) = 0;
    virtual hitnlls::matrix::VectorXd GetVector() const = 0;

    void SetValue(const ValueType &val) { val_ = val; }
    ValueType GetValue() const { return val_; }

protected:
    bool CheckDim(const hitnlls::matrix::VectorXd &v) const {
        if (v.Size() < ndim) {
            return false;
        } else {
            return true;
        }
    }

private:
    ValueType val_;
};

} // namespace hitcadmm