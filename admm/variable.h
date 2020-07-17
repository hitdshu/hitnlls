#pragma once

#include "../matrix/dense.h"

namespace nlls {

class Variable {
public:
    explicit Variable(int ndim) { ndim_ = ndim; }
    virtual void Project() = 0;
    virtual void SetVector(const VectorXf &v) = 0;
    virtual VectorXf GetVector() const = 0;
    int GetDim() const { return ndim_; }
private:
    int ndim_;
};

template <int N, typename ValueType>
class VariableImp : public Variable {
public:
    explicit VariableImp() : Variable(N) {}
    virtual void Project() = 0;
    virtual void SetVector(const VectorXf &v) = 0;
    virtual VectorXf GetVector() const = 0;
    void SetValue(const ValueType &val) { val_ = val; }
    const ValueType &GetValue() const { return val_; }
protected:
    bool CheckDim(const VectorXf &v) const { if (v.Size() < N) { return false; } else { return true; } }
private:
    ValueType val_;
};

template <typename ValueType>
class VariableImpX : public Variable {
public:
    explicit VariableImpX(int n) : Variable(n) {}
    virtual void Project() = 0;
    virtual void SetVector(const VectorXf &v) = 0;
    virtual VectorXf GetVector() const = 0;
    void SetValue(const ValueType &val) { val_ = val; }
    const ValueType &GetValue() const { return val_; }
protected:
    bool CheckDim(const VectorXf &v) const { if (v.Size() < Variable::GetDim()) { return false; } else { return true; } }
private:
    ValueType val_;
};

} // namespace nlls