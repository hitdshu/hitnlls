#pragma once

#include "matrix/base_matrix.h"
#include <memory>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <cassert>

namespace hitnlls {
namespace matrix {

template <typename ValueType>
class Matrixx : public BaseMatrix<ValueType> {
public:
    Matrixx();
    Matrixx(int nrows, int ncols = 1);
    Matrixx(const Matrixx &matx);

    inline ValueType &operator()(int ridx, int cidx = 0) override;
    inline ValueType operator()(int ridx, int cidx = 0) const override;
    inline ValueType &operator[](int didx) override;
    inline ValueType operator[](int didx) const override;
    Matrixx<ValueType> operator+(const Matrixx<ValueType> &matx) const;
    Matrixx<ValueType> operator+(const ValueType &val) const;
    Matrixx<ValueType> operator-(const Matrixx<ValueType> &matx) const;
    Matrixx<ValueType> operator-() const;
    Matrixx<ValueType> operator-(const ValueType &val) const;
    Matrixx<ValueType> operator*(const Matrixx<ValueType> &matx) const;
    Matrixx<ValueType> operator*(const ValueType &val) const;
    Matrixx<ValueType> operator/(const Matrixx<ValueType> &matx) const;
    Matrixx<ValueType> &operator+=(const Matrixx<ValueType> &matx);
    Matrixx<ValueType> &operator+=(const ValueType &val);
    Matrixx<ValueType> &operator-=(const Matrixx<ValueType> &matx);
    Matrixx<ValueType> &operator-=(const ValueType &val);
    Matrixx<ValueType> &operator*=(const Matrixx<ValueType> &matx);
    Matrixx<ValueType> &operator/=(const Matrixx<ValueType> &matx);
    Matrixx<ValueType> &operator=(const Matrixx<ValueType> &matx);
    Matrixx<ValueType> &operator=(const ValueType &val);

    ValueType Norm() const;
    Matrixx<ValueType> Transpose() const;
    Matrixx<ValueType> Block(int rsidx, int csidx, int rows, int cols) const;
    Matrixx<ValueType> &SetBlock(int rsidx, int csidx, int rows, int cols, const Matrixx<ValueType> &matx);
    Matrixx<ValueType> SolveAxb(const Matrixx<ValueType> &b) const;
    Matrixx<ValueType> Inverse() const;
    Matrixx<ValueType> DiagonalInverse() const;
    ValueType MaxDiagonalValue() const;

    void SetLowtri();
    bool CholeskyLLT(Matrixx<ValueType> &chol) const;
    Matrixx<ValueType> CholeskyLLT() const;
    bool LUPDecomp(Matrixx<ValueType> &l, Matrixx<ValueType> &u, Matrixx<ValueType> &p) const;
    Matrixx<ValueType> LUPSolve(const Matrixx<ValueType> &l, const Matrixx<ValueType> &u, const Matrixx<ValueType> &p, const Matrixx<ValueType> &b) const;

    inline MatrixSize Size() const override;
    inline int Rows() const override;
    inline int Cols() const override;
    void Print(const ::std::string &name = "") const override;

    static Matrixx<ValueType> Identity(int nrows);

    friend Matrixx<ValueType> operator+(const ValueType &val, const Matrixx<ValueType> &matx) { return matx.operator+(val); }
    friend Matrixx<ValueType> operator-(const ValueType &val, const Matrixx<ValueType> &matx) { return matx.operator-(val); }
    friend Matrixx<ValueType> operator*(const ValueType &val, const Matrixx<ValueType> &matx) { return matx.operator*(val); }
    friend ::std::ostream &operator<<(::std::ostream &out, const Matrixx<ValueType> &matx) { matx.Print(); return out; }

protected:
    void Init(int nrows, int ncols);

private:
    ::std::unique_ptr<ValueType[]> data_;
    int nrows_;
    int ncols_;
};

template <typename ValueType>
Matrixx<ValueType>::Matrixx() {
    nrows_ = -1;
    ncols_ = -1;
    data_.reset();
}

template <typename ValueType>
Matrixx<ValueType>::Matrixx(int nrows, int ncols) {
    assert(nrows > 0);
    assert(ncols > 0);
    Init(nrows, ncols);
}

template <typename ValueType>
Matrixx<ValueType>::Matrixx(const Matrixx<ValueType> &matx) {
    Init(matx.nrows_, matx.ncols_);
    memcpy(data_.get(), matx.data_.get(), nrows_ * ncols_ * sizeof(ValueType));
}

template <typename ValueType>
Matrixx<ValueType> &Matrixx<ValueType>::operator=(const Matrixx<ValueType> &matx) {
    if (&matx == this) {
        return *this;
    }
    Init(matx.nrows_, matx.ncols_);
    memcpy(data_.get(), matx.data_.get(), nrows_ * ncols_ * sizeof(ValueType));
    return *this;
}

template <typename ValueType>
Matrixx<ValueType> &Matrixx<ValueType>::operator=(const ValueType &val) {
    if (-1 == nrows_ || -1 == ncols_) {
        return *this;
    }
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = 0; cidx < ncols_; ++cidx) {
            this->operator()(ridx, cidx) = val;
        }
    }
    return *this;
}

template <typename ValueType>
ValueType &Matrixx<ValueType>::operator()(int ridx, int cidx) {
    assert(ridx < nrows_ && 0 <= ridx);
    assert(cidx < ncols_ && 0 <= cidx);
    int didx = cidx + ridx * ncols_;
    return data_.get()[didx];
}

template <typename ValueType>
ValueType Matrixx<ValueType>::operator()(int ridx, int cidx) const {
    assert(ridx < nrows_ && 0 <= ridx);
    assert(cidx < ncols_ && 0 <= cidx);
    int didx = cidx + ridx * ncols_;
    return data_.get()[didx];
}

template <typename ValueType>
ValueType &Matrixx<ValueType>::operator[](int didx) {
    assert(didx < nrows_ * ncols_ && 0 <= didx);
    return data_.get()[didx];
}

template <typename ValueType>
ValueType Matrixx<ValueType>::operator[](int didx) const {
    assert(didx < nrows_ * ncols_ && 0 <= didx);
    return data_.get()[didx];
}

template <typename ValueType>
Matrixx<ValueType> Matrixx<ValueType>::operator+(const Matrixx<ValueType> &matx) const {
    if (nrows_ == -1 || ncols_ == -1) {
        return matx;
    }
    assert(nrows_ == matx.nrows_ && ncols_ == matx.ncols_);
    Matrixx<ValueType> result(nrows_, ncols_);
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = 0; cidx < ncols_; ++cidx) {
            result(ridx, cidx) = this->operator()(ridx, cidx) + matx(ridx, cidx);
        }
    }
    return result;
}

template <typename ValueType>
Matrixx<ValueType> Matrixx<ValueType>::operator+(const ValueType &val) const {
    Matrixx<ValueType> result(nrows_, ncols_);
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = 0; cidx < ncols_; ++cidx) {
            result(ridx, cidx) = this->operator()(ridx, cidx) + val;
        }
    }
    return result;
}

template <typename ValueType>
Matrixx<ValueType> Matrixx<ValueType>::operator-(const Matrixx<ValueType> &matx) const {
    if (nrows_ == -1 || ncols_ == -1) {
        return - matx;
    }
    assert(nrows_ == matx.nrows_ && ncols_ == matx.ncols_);
    Matrixx<ValueType> result(nrows_, ncols_);
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = 0; cidx < ncols_; ++cidx) {
            result(ridx, cidx) = this->operator()(ridx, cidx) - matx(ridx, cidx);
        }
    }
    return result;
}

template <typename ValueType>
Matrixx<ValueType> Matrixx<ValueType>::operator-() const {
    Matrixx<ValueType> result(nrows_, ncols_);
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = 0; cidx < ncols_; ++cidx) {
            result(ridx, cidx) = - this->operator()(ridx, cidx);
        }
    }
    return result;
}

template <typename ValueType>
Matrixx<ValueType> Matrixx<ValueType>::operator-(const ValueType &val) const {
    Matrixx<ValueType> result(nrows_, ncols_);
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = 0; cidx < ncols_; ++cidx) {
            result(ridx, cidx) = this->operator()(ridx, cidx) - val;
        }
    }
    return result;
}

template <typename ValueType>
Matrixx<ValueType> Matrixx<ValueType>::operator*(const Matrixx<ValueType> &matx) const {
    if (ncols_ == matx.nrows_) {
        Matrixx<ValueType> result(nrows_, matx.Cols());
        for (int ridx = 0; ridx < nrows_; ++ridx) {
            for (int cidx = 0; cidx < matx.Cols(); ++cidx) {
                result(ridx, cidx) = 0;
                for (int iidx = 0; iidx < ncols_; ++iidx) {
                    result(ridx, cidx) += this->operator()(ridx, iidx) * matx(iidx, cidx);
                }
            }
        }
        return result;
    } else if (nrows_ == matx.nrows_ && ncols_ == matx.ncols_) {
        ValueType sum = 0;
        for (int ridx = 0; ridx < nrows_; ++ridx) {
            for (int cidx = 0; cidx < ncols_; ++cidx) {
                sum += this->operator()(ridx, cidx) * matx(ridx, cidx);
            }
        }
        Matrixx<ValueType> result(1, 1);
        result(0, 0) = sum;
        return result;
    } else {
        return Matrixx<ValueType>();
    }
}

template <typename ValueType>
Matrixx<ValueType> Matrixx<ValueType>::operator/(const Matrixx<ValueType> &matx) const {
    assert(matx.nrows_ == ncols_ && matx.nrows_ == matx.ncols_);
    return this->operator*(matx.Inverse());
}

template <typename ValueType>
Matrixx<ValueType> Matrixx<ValueType>::operator*(const ValueType &val) const {
    Matrixx<ValueType> result(nrows_, ncols_);
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = 0; cidx < ncols_; ++cidx) {
            result(ridx, cidx) = this->operator()(ridx, cidx) * val;
        }
    }
    return result;
}

template <typename ValueType>
Matrixx<ValueType> &Matrixx<ValueType>::operator+=(const Matrixx<ValueType> &matx) {
    if (nrows_ == -1 || ncols_ == -1) {
        this->operator=(matx);
        return *this;
    }
    assert(nrows_ == matx.nrows_ && ncols_ == matx.ncols_);
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = 0; cidx < ncols_; ++cidx) {
            this->operator()(ridx, cidx) += matx(ridx, cidx);
        }
    }
    return *this;
}

template <typename ValueType>
Matrixx<ValueType> &Matrixx<ValueType>::operator+=(const ValueType &val) {
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = 0; cidx < ncols_; ++cidx) {
            this->operator()(ridx, cidx) += val;
        }
    }
    return *this;
}

template <typename ValueType>
Matrixx<ValueType> &Matrixx<ValueType>::operator-=(const Matrixx<ValueType> &matx) {
    if (nrows_ == -1 || ncols_ == -1) {
        this->operator=(- matx);
        return *this;
    }
    assert(nrows_ == matx.nrows_ && ncols_ == matx.ncols_);
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = 0; cidx < ncols_; ++cidx) {
            this->operator()(ridx, cidx) -= matx(ridx, cidx);
        }
    }
    return *this;
}

template <typename ValueType>
Matrixx<ValueType> &Matrixx<ValueType>::operator*=(const Matrixx<ValueType> &matx) {
    assert(matx.nrows_ == ncols_ && matx.nrows_ == matx.ncols_);
    Matrixx<ValueType> tmp(nrows_, ncols_);
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = 0; cidx < ncols_; ++cidx) {
            tmp(ridx, cidx) = 0;
            for (int iidx = 0; iidx < ncols_; ++iidx) {
                tmp(ridx, cidx) += this->operator()(ridx, iidx) * matx(iidx, cidx);
            }
        }
    }
    this->operator=(tmp);
    return *this;
}

template <typename ValueType>
Matrixx<ValueType> &Matrixx<ValueType>::operator/=(const Matrixx<ValueType> &matx) {
    assert(matx.nrows_ == matx.ncols_ && matx.nrows_ == matx.ncols_);
    this->operator*=(matx.Inverse());
    return *this;
}

template <typename ValueType>
Matrixx<ValueType> &Matrixx<ValueType>::operator-=(const ValueType &val) {
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = 0; cidx < ncols_; ++cidx) {
            this->operator()(ridx, cidx) -= val;
        }
    }
    return *this;
}

template <typename ValueType>
ValueType Matrixx<ValueType>::Norm() const {
    ValueType norm = 0;
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = 0; cidx < ncols_; ++cidx) {
            norm += this->operator()(ridx, cidx) * this->operator()(ridx, cidx);
        }
    }
    return sqrt(norm);
}

template <typename ValueType>
Matrixx<ValueType> Matrixx<ValueType>::Transpose() const {
    Matrixx<ValueType> result(ncols_, nrows_);
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = 0; cidx < ncols_; ++cidx) {
            result(cidx, ridx) = this->operator()(ridx, cidx);
        }
    }
    return result;
}

template <typename ValueType>
Matrixx<ValueType> Matrixx<ValueType>::Block(int rsidx, int csidx, int rows, int cols) const {
    assert(rsidx + rows <= nrows_ && csidx + cols <= ncols_);
    Matrixx<ValueType> result(rows, cols);
    for (int ridx = 0; ridx < rows; ++ridx) {
        for (int cidx = 0; cidx < cols; ++cidx) {
            result(cidx, ridx) = this->operator()(ridx + rsidx, cidx + csidx);
        }
    }
    return result;
}

template <typename ValueType>
Matrixx<ValueType> &Matrixx<ValueType>::SetBlock(int rsidx, int csidx, int rows, int cols, const Matrixx<ValueType> &matx) {
    assert(matx.Rows() >= rows && matx.Cols() >= cols);
    assert(rsidx + rows <= nrows_ && csidx + cols <= ncols_);
    for (int ridx = 0; ridx < rows; ++ridx) {
        for (int cidx = 0; cidx < cols; ++cidx) {
            this->operator()(ridx + rsidx, cidx + csidx) = matx(ridx, cidx);
        }
    }
    return *this;
}

template <typename ValueType>
ValueType Matrixx<ValueType>::MaxDiagonalValue() const {
    ValueType result = 0;
    for (int idx = 0; idx < ::std::min(nrows_, ncols_); ++idx) {
        if (this->operator()(idx, idx) > result) {
            result = this->operator()(idx, idx);
        }
    }
    return result;
}

template <typename ValueType>
Matrixx<ValueType> Matrixx<ValueType>::SolveAxb(const Matrixx<ValueType> &b) const {
    assert(nrows_ == ncols_ && nrows_ == b.nrows_ && b.ncols_ == 1);
    Matrixx<ValueType> l;
    Matrixx<ValueType> u;
    Matrixx<ValueType> p;
    if (LUPDecomp(l, u, p)) {
        return LUPSolve(l, u, p, b);
    }
    return Matrixx<ValueType>();
}

template <typename ValueType>
Matrixx<ValueType> Matrixx<ValueType>::Inverse() const {
    assert(nrows_ == ncols_);
    Matrixx<ValueType> result(nrows_, nrows_);
    Matrixx<ValueType> l;
    Matrixx<ValueType> u;
    Matrixx<ValueType> p;
    if (LUPDecomp(l, u, p)) {
        for (int cidx = 0; cidx < nrows_; ++cidx) {
            Matrixx<ValueType> b(nrows_);
            b[cidx] = 1;
            Matrixx<ValueType> vtmp = LUPSolve(l, u, p, b);
            for (int ridx = 0; ridx < nrows_; ++ridx) {
                result(ridx, cidx) = vtmp[ridx];
            }
        }
    }
    return result;
}

template <typename ValueType>
Matrixx<ValueType> Matrixx<ValueType>::DiagonalInverse() const {
    assert(nrows_ == ncols_);
    Matrixx<ValueType> result(nrows_, nrows_);
    for (int idx = 0; idx < nrows_; ++idx) {
        result(idx, idx) = 1 / this->operator()(idx, idx);
    }
    return result;
}

template <typename ValueType>
void Matrixx<ValueType>::SetLowtri() {
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = ridx + 1; cidx < nrows_; ++cidx) {
             this->operator()(ridx, cidx) = 0;
        }
    }
}

template <typename ValueType>
bool Matrixx<ValueType>::CholeskyLLT(Matrixx<ValueType> &chol) const {
    assert(nrows_ == ncols_);
    chol = Matrixx<ValueType>(*this);
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        ValueType dii = chol(ridx, ridx);
        if (dii < 0) {
            ::std::cout << "Matrix of cholesky is not positive definite!" << ::std::endl;
            return false;
        }
        ValueType diis = sqrt(dii);
        for (int sridx = ridx; sridx < nrows_; ++sridx) {
            chol(sridx, ridx) /= diis;
        }
        for (int sridx = ridx + 1; sridx < nrows_; ++sridx) {
            for (int scidx = ridx + 1; scidx <= sridx; ++scidx) {
                chol(sridx, scidx) -= chol(sridx, ridx) * chol(scidx, ridx);
            }
        }
    }
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = ridx + 1; cidx < nrows_; ++cidx) {
            chol(ridx, cidx) = 0;
        }
    }
    return true;
}

template <typename ValueType>
Matrixx<ValueType> Matrixx<ValueType>::CholeskyLLT() const {
    assert(nrows_ == ncols_);
    Matrixx<ValueType> chol = Matrixx<ValueType>(*this);
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        ValueType dii = chol(ridx, ridx);
        ValueType diis = sqrt(dii);
        for (int sridx = ridx; sridx < nrows_; ++sridx) {
            chol(sridx, ridx) /= diis;
        }
        for (int sridx = ridx + 1; sridx < nrows_; ++sridx) {
            for (int scidx = ridx + 1; scidx <= sridx; ++scidx) {
                chol(sridx, scidx) -= chol(sridx, ridx) * chol(scidx, ridx);
            }
        }
    }
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = ridx + 1; cidx < nrows_; ++cidx) {
            chol(ridx, cidx) = 0;
        }
    }
    return chol;
}

template <typename ValueType>
bool Matrixx<ValueType>::LUPDecomp(Matrixx<ValueType> &l, Matrixx<ValueType> &u, Matrixx<ValueType> &p) const {
    assert(nrows_ == ncols_);
    l = Matrixx<ValueType>(nrows_, nrows_);
    u = Matrixx<ValueType>(nrows_, nrows_);
    p = Matrixx<ValueType>(nrows_);
    for (int idx = 0; idx < nrows_; ++idx) {
        p[idx] = idx;
    }
    Matrixx<ValueType> mtmp(*this);
    for (int rcidx = 0; rcidx < nrows_; ++rcidx) {
        int pr = -1;
        for (int ridx = rcidx; ridx < nrows_; ++ridx) {
            if (fabs(mtmp(ridx, rcidx)) > 1e-5) {
                pr = ridx;
                break;
            }
        }
        std::swap(p[rcidx], p[pr]);
        for (int cidx = 0; cidx < nrows_; ++cidx) {
            std::swap(mtmp(rcidx, cidx), mtmp(pr, cidx));
        }
        ValueType dii = mtmp(rcidx, rcidx);
        if (fabs(dii) < 1e-5) {
            ::std::cout << "Matrix of LUDecomp is not positive in pivot!" << ::std::endl;
            return false;
        }
        l(rcidx, rcidx) = 1;
        for (int cidx = rcidx; cidx < nrows_; ++cidx) {
            u(rcidx, cidx) = mtmp(rcidx, cidx);
            l(cidx, rcidx) = mtmp(cidx, rcidx) / dii;
        }
        for (int sridx = rcidx + 1; sridx < nrows_; ++sridx) {
            for (int scidx = rcidx + 1; scidx < nrows_; ++scidx) {
                mtmp(sridx, scidx) -= l(sridx, rcidx) * u(rcidx, scidx);
            }
        }
    }
    return true;
}

template <typename ValueType>
Matrixx<ValueType> Matrixx<ValueType>::LUPSolve(const Matrixx<ValueType> &l, const Matrixx<ValueType> &u, const Matrixx<ValueType> &p, const Matrixx<ValueType> &b) const {
    Matrixx<ValueType> x(l.nrows_, 1);
    Matrixx<ValueType> y(l.nrows_, 1);
    for (int ridx = 0; ridx < l.nrows_; ++ridx) {
        y[ridx] = b[p[ridx]];
        for (int cidx = 0; cidx < ridx; ++cidx) {
            y[ridx] = y[ridx] - l(ridx, cidx) * y[cidx];
        }
    }
    for (int ridx = nrows_ - 1; ridx >= 0; --ridx) {
        x[ridx] = y[ridx];
        for (int cidx = nrows_ - 1; cidx > ridx; --cidx) {
            x[ridx] = x[ridx] - u(ridx, cidx) * x[cidx];
        }
        x[ridx] /= u(ridx, ridx);
    }
    return x;
}

template <typename ValueType>
MatrixSize Matrixx<ValueType>::Size() const {
    MatrixSize size;
    size.rows = nrows_;
    size.cols = ncols_;
    return size;
}

template <typename ValueType>
int Matrixx<ValueType>::Rows() const {
    return nrows_;
}

template <typename ValueType>
int Matrixx<ValueType>::Cols() const {
    return ncols_;
}

template <typename ValueType>
void Matrixx<ValueType>::Print(const ::std::string &name) const {
    if (nrows_ == -1 || ncols_ == -1) {
        ::std::cout << "--Matrix " << name << " is empty!" << ::std::endl;
        return;
    }
    ::std::cout << "--Matrix " << name << " size: rows " << nrows_ << ", cols " << ncols_ << ::std::endl;
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        for (int cidx = 0; cidx < ncols_; ++cidx) {
            ::std::cout << ::std::setw(12) << this->operator()(ridx, cidx) << " ";
        }
        ::std::cout << ::std::endl;
    }
}

template <typename ValueType>
Matrixx<ValueType> Matrixx<ValueType>::Identity(int nrows) {
    Matrixx<ValueType> result(nrows, nrows);
    for (int idx = 0; idx < nrows; ++idx) {
        result(idx, idx) = 1;
    }
    return result;
}

template <typename ValueType>
void Matrixx<ValueType>::Init(int nrows, int ncols) {
    nrows_ = nrows;
    ncols_ = ncols;
    data_.reset(new ValueType[nrows * ncols]);
    memset(data_.get(), 0, nrows * ncols * sizeof(ValueType));
}

using Matrixxf = Matrixx<float>;
using Matrixxd = Matrixx<double>;
using Matrixxi = Matrixx<int>;

} // namespace matrix
} // namespace hitnlls