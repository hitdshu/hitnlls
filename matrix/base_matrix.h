#pragma once

#include <string>

namespace hitnlls {
namespace matrix {

struct MatrixSize {
    int rows;
    int cols;
};

template <typename ValueType>
class BaseMatrix {
public:
    BaseMatrix() = default;
    virtual ~BaseMatrix() = default;

    virtual inline ValueType &operator()(int ridx, int cidx = 0) = 0;
    virtual inline ValueType operator()(int ridx, int cidx = 0) const = 0;
    virtual inline ValueType &operator[](int didx) = 0;
    virtual inline ValueType operator[](int didx) const = 0;

    virtual inline MatrixSize Size() const = 0;
    virtual inline int Rows() const = 0;
    virtual inline int Cols() const = 0;
    virtual void Print(const ::std::string &name = "") const = 0;
};

} // namespace matrix
} // namespace hitnlls