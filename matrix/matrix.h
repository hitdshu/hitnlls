#pragma once

#include "matrix/matrixx.h"

namespace hitnlls {
namespace matrix {

template <typename ValueType, int nrows, int ncols>
class Matrix : public Matrixx<ValueType> {
public:
    Matrix() : Matrixx<ValueType>(nrows, ncols) {}
    Matrix(const Matrix &mat) : Matrixx<ValueType>(mat) {}
    Matrix(const Matrixx<ValueType> &matx) { this->operator=(matx); }
    Matrix(const ValueType &val) { Matrixx<ValueType>::Init(nrows, ncols); this->operator=(val); }

    inline ValueType &operator()(int ridx, int cidx = 0) override;
    inline ValueType operator()(int ridx, int cidx = 0) const override;
    inline ValueType &operator[](int didx) override;
    inline ValueType operator[](int didx) const override;
    Matrix operator+(const Matrix &mat) const;
    Matrix operator+(const ValueType &val) const;
    Matrixx<ValueType> operator+(const Matrixx<ValueType> &matx) const;
    Matrix operator-(const Matrix &mat) const;
    Matrix operator-() const;
    Matrix operator-(const ValueType &val) const;
    Matrixx<ValueType> operator-(const Matrixx<ValueType> &matx) const;
    template <int sncols>
    Matrix<ValueType, ncols, sncols> operator*(const Matrix<ValueType, ncols, sncols> &mat) const;
    Matrixx<ValueType> operator*(const Matrixx<ValueType> &matx) const;
    Matrix &operator=(const ValueType &val) { Matrixx<ValueType>::operator=(val); return *this; }
    Matrix &operator=(const Matrixx<ValueType> &matx) { assert(matx.Rows() == nrows && matx.Cols() == ncols); Matrixx<ValueType>::operator=(matx); return *this; }
    Matrix operator*(const ValueType &val) const;
    Matrix &operator+=(const Matrix &mat);
    Matrix &operator+=(const ValueType &val);
    Matrix &operator-=(const Matrix &mat);
    Matrix &operator-=(const ValueType &val);
    Matrix &operator+=(const Matrixx<ValueType> &matx);
    Matrix &operator-=(const Matrixx<ValueType> &matx);
    Matrix &operator*=(const Matrixx<ValueType> &matx);
    Matrix &operator/=(const Matrixx<ValueType> &matx);

    static Matrix Identity();

    friend ::std::ostream &operator<<(::std::ostream &out, const Matrix<ValueType, nrows, ncols> &mat) { mat.Print(); return out; }
};

template <typename ValueType, int nrows, int ncols>
ValueType &Matrix<ValueType, nrows, ncols>::operator()(int ridx, int cidx) {
    return Matrixx<ValueType>::operator()(ridx, cidx);
}

template <typename ValueType, int nrows, int ncols>
ValueType Matrix<ValueType, nrows, ncols>::operator()(int ridx, int cidx) const {
    return Matrixx<ValueType>::operator()(ridx, cidx);
}

template <typename ValueType, int nrows, int ncols>
ValueType &Matrix<ValueType, nrows, ncols>::operator[](int didx) {
    return Matrixx<ValueType>::operator[](didx);
}

template <typename ValueType, int nrows, int ncols>
ValueType Matrix<ValueType, nrows, ncols>::operator[](int didx) const {
    return Matrixx<ValueType>::operator[](didx);
}

template <typename ValueType, int nrows, int ncols>
Matrix<ValueType, nrows, ncols> Matrix<ValueType, nrows, ncols>::operator+(const Matrix<ValueType, nrows, ncols> &mat) const {
    Matrix<ValueType, nrows, ncols> result;
    for (int ridx = 0; ridx < nrows; ++ridx) {
        for (int cidx = 0; cidx < ncols; ++cidx) {
            result(ridx, cidx) = this->operator()(ridx, cidx) + mat(ridx, cidx);
        }
    }
    return result;
}

template <typename ValueType, int nrows, int ncols>
Matrix<ValueType, nrows, ncols> Matrix<ValueType, nrows, ncols>::operator+(const ValueType &val) const {
    Matrix<ValueType, nrows, ncols> result;
    for (int ridx = 0; ridx < nrows; ++ridx) {
        for (int cidx = 0; cidx < ncols; ++cidx) {
            result(ridx, cidx) = this->operator()(ridx, cidx) + val;
        }
    }
    return result;
}

template <typename ValueType, int nrows, int ncols>
Matrixx<ValueType> Matrix<ValueType, nrows, ncols>::operator+(const Matrixx<ValueType> &matx) const {
    return Matrixx<ValueType>::operator+(matx);
}

template <typename ValueType, int nrows, int ncols>
Matrix<ValueType, nrows, ncols> Matrix<ValueType, nrows, ncols>::operator-(const Matrix<ValueType, nrows, ncols> &mat) const {
    Matrix<ValueType, nrows, ncols> result;
    for (int ridx = 0; ridx < nrows; ++ridx) {
        for (int cidx = 0; cidx < ncols; ++cidx) {
            result(ridx, cidx) = this->operator()(ridx, cidx) - mat(ridx, cidx);
        }
    }
    return result;
}

template <typename ValueType, int nrows, int ncols>
Matrix<ValueType, nrows, ncols> Matrix<ValueType, nrows, ncols>::operator-() const {
    Matrix<ValueType, nrows, ncols> result;
    for (int ridx = 0; ridx < nrows; ++ridx) {
        for (int cidx = 0; cidx < ncols; ++cidx) {
            result(ridx, cidx) = - this->operator()(ridx, cidx);
        }
    }
    return result;
}

template <typename ValueType, int nrows, int ncols>
Matrix<ValueType, nrows, ncols> Matrix<ValueType, nrows, ncols>::operator-(const ValueType &val) const {
    Matrix<ValueType, nrows, ncols> result;
    for (int ridx = 0; ridx < nrows; ++ridx) {
        for (int cidx = 0; cidx < ncols; ++cidx) {
            result(ridx, cidx) = this->operator()(ridx, cidx) - val;
        }
    }
    return result;
}

template <typename ValueType, int nrows, int ncols>
Matrixx<ValueType> Matrix<ValueType, nrows, ncols>::operator-(const Matrixx<ValueType> &matx) const {
    return Matrixx<ValueType>::operator-(matx);
}

template <typename ValueType, int nrows, int ncols>
template <int sncols>
Matrix<ValueType, ncols, sncols> Matrix<ValueType, nrows, ncols>::operator*(const Matrix<ValueType, ncols, sncols> &mat) const {
    Matrix<ValueType, nrows, sncols> result;
    for (int ridx = 0; ridx < nrows; ++ridx) {
        for (int cidx = 0; cidx < sncols; ++cidx) {
            result(ridx, cidx) = 0;
            for (int iidx = 0; iidx < ncols; ++iidx) {
                result(ridx, cidx) += this->operator()(ridx, iidx) * mat(iidx, cidx);
            }
        }
    }
    return result;
}

template <typename ValueType, int nrows, int ncols>
Matrix<ValueType, nrows, ncols> Matrix<ValueType, nrows, ncols>::operator*(const ValueType &val) const {
    Matrix<ValueType, nrows, ncols> result;
    for (int ridx = 0; ridx < nrows; ++ridx) {
        for (int cidx = 0; cidx < ncols; ++cidx) {
            result(ridx, cidx) = this->operator()(ridx, cidx) * val;
        }
    }
    return result;
}

template <typename ValueType, int nrows, int ncols>
Matrixx<ValueType> Matrix<ValueType, nrows, ncols>::operator*(const Matrixx<ValueType> &matx) const {
    return Matrixx<ValueType>::operator*(matx);
}

template <typename ValueType, int nrows, int ncols>
Matrix<ValueType, nrows, ncols> &Matrix<ValueType, nrows, ncols>::operator+=(const Matrix<ValueType, nrows, ncols> &mat) {
    for (int ridx = 0; ridx < nrows; ++ridx) {
        for (int cidx = 0; cidx < ncols; ++cidx) {
            this->operator()(ridx, cidx) += mat(ridx, cidx);
        }
    }
    return *this;
}

template <typename ValueType, int nrows, int ncols>
Matrix<ValueType, nrows, ncols> &Matrix<ValueType, nrows, ncols>::operator+=(const ValueType &val) {
    for (int ridx = 0; ridx < nrows; ++ridx) {
        for (int cidx = 0; cidx < ncols; ++cidx) {
            this->operator()(ridx, cidx) += val;
        }
    }
    return *this;
}

template <typename ValueType, int nrows, int ncols>
Matrix<ValueType, nrows, ncols> &Matrix<ValueType, nrows, ncols>::operator-=(const Matrix<ValueType, nrows, ncols> &mat) {
    for (int ridx = 0; ridx < nrows; ++ridx) {
        for (int cidx = 0; cidx < ncols; ++cidx) {
            this->operator()(ridx, cidx) -= mat(ridx, cidx);
        }
    }
    return *this;
}

template <typename ValueType, int nrows, int ncols>
Matrix<ValueType, nrows, ncols> &Matrix<ValueType, nrows, ncols>::operator-=(const ValueType &val) {
    for (int ridx = 0; ridx < nrows; ++ridx) {
        for (int cidx = 0; cidx < ncols; ++cidx) {
            this->operator()(ridx, cidx) -= val;
        }
    }
    return *this;
}

template <typename ValueType, int nrows, int ncols>
Matrix<ValueType, nrows, ncols> &Matrix<ValueType, nrows, ncols>::operator+=(const Matrixx<ValueType> &matx) {
    Matrixx<ValueType>::operator+=(matx);
    return *this;
}

template <typename ValueType, int nrows, int ncols>
Matrix<ValueType, nrows, ncols> &Matrix<ValueType, nrows, ncols>::operator-=(const Matrixx<ValueType> &matx) {
    Matrixx<ValueType>::operator+=(matx);
    return *this;
}

template <typename ValueType, int nrows, int ncols>
Matrix<ValueType, nrows, ncols> &Matrix<ValueType, nrows, ncols>::operator*=(const Matrixx<ValueType> &matx) {
    Matrixx<ValueType>::operator+=(matx);
    return *this;
}

template <typename ValueType, int nrows, int ncols>
Matrix<ValueType, nrows, ncols> &Matrix<ValueType, nrows, ncols>::operator/=(const Matrixx<ValueType> &matx) {
    Matrixx<ValueType>::operator+=(matx);
    return *this;
}

template <typename ValueType, int nrows, int ncols>
Matrix<ValueType, nrows, ncols> Matrix<ValueType, nrows, ncols>::Identity() {
    Matrix mat;
    for (int idx = 0; idx < nrows && idx < ncols; ++idx) {
        mat(idx, idx) = 1;
    }
    return mat;
}

using Matrix33f = Matrix<float, 3, 3>;
using Matrix33d = Matrix<double, 3, 3>;
using Matrix33i = Matrix<int, 3, 3>;
using Matrix22f = Matrix<float, 2, 2>;
using Matrix22d = Matrix<double, 2, 2>;
using Matrix22i = Matrix<int, 2, 2>;
using Matrix44f = Matrix<float, 4, 4>;
using Matrix44d = Matrix<double, 4, 4>;
using Matrix44i = Matrix<int, 4, 4>;

using Vector3f = Matrix<float, 3, 1>;
using Vector3d = Matrix<double, 3, 1>;
using Vector3i = Matrix<int, 3, 1>;
using Vector2f = Matrix<float, 2, 1>;
using Vector2d = Matrix<double, 2, 1>;
using Vector2i = Matrix<int, 2, 1>;
using Vector4f = Matrix<float, 4, 1>;
using Vector4d = Matrix<double, 4, 1>;
using Vector4i = Matrix<int, 4, 1>;
using Vector1f = Matrix<float, 1, 1>;
using Vector1d = Matrix<double, 1, 1>;
using Vector1i = Matrix<int, 1, 1>;

} // namespace matrix
} // namespace hitnlls