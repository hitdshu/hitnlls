#pragma once

#include "matrix/vecxs.h"
#include "matrix/mtraits.h"
#include <vector>

namespace hitnlls {
namespace matrix {

template <typename ValueType>
class Matrixs {
public:
    Matrixs() { nrows_ = -1; ncols_ = -1; lm_lambda_ = 0.0; }
    Matrixs(int nrows, int ncols) { Init(nrows, ncols); lm_lambda_ = 0.0; }
    Matrixs(int nrows) { Init(nrows, nrows); lm_lambda_ = 0.0; }
    Matrixs(const Matrixs &smat) { nrows_ = smat.nrows_; ncols_ = smat.ncols_; rvecs_ = smat.rvecs_; }

    inline ValueType &operator()(int ridx, int cidx);
    inline ValueType operator()(int ridx, int cidx) const;
    inline Vecxs<ValueType> &operator[](int ridx) { return rvecs_[ridx]; }
    inline Vecxs<ValueType> operator[](int ridx) const { return rvecs_[ridx]; }

    inline int Rows() const { return nrows_; }
    inline int Cols() const { return ncols_; }

    inline void Delete(int ridx, int cidx) { rvecs_[ridx].Delete(cidx); }
    inline bool HasValue(int ridx, int cidx) const { return rvecs_[ridx].HasValue(cidx); }
    inline int NumElementsInRow(int ridx, int cidx = 0) const { return rvecs_[ridx].GetIndicesget(cidx).size(); }
    inline ::std::vector<int> ElementsInRow(int ridx, int cidx = 0) const { return rvecs_[ridx].GetIndicesget(cidx); }
    inline void SetLambdalm(double lm) { lm_lambda_ = lm; }

    Matrixs<ValueType> operator+(const Matrixs<ValueType> &smat) const;
    Matrixs<ValueType> operator-(const Matrixs<ValueType> &smat) const;
    Matrixs<ValueType> operator*(const Matrixs<ValueType> &smat) const;

    Matrixs<ValueType> Transpose() const;
    Matrixs<ValueType> CholeskyLLT() const;
    Vecxs<ValueType> SolveWithlm(const Vecxs<ValueType> &b) const;

    void Print(const ::std::string &name = "") const;

    friend ::std::ostream &operator<<(::std::ostream &out, const Matrixs<ValueType> &smat) { smat.Print(); return out; }

protected:
    void Init(int nrows, int ncols);
    Vecxs<ValueType> ForwardSubst(const Matrixs<ValueType> &choll, const Vecxs<ValueType> &b) const;
    Vecxs<ValueType> BackwardSubst(const Matrixs<ValueType> &choll, const Vecxs<ValueType> &y) const;

private:
    int nrows_;
    int ncols_;
    std::vector<Vecxs<ValueType> > rvecs_;
    double lm_lambda_;
};

template <typename ValueType>
ValueType &Matrixs<ValueType>::operator()(int ridx, int cidx) {
    return rvecs_[ridx](cidx);
}

template <typename ValueType>
ValueType Matrixs<ValueType>::operator()(int ridx, int cidx) const {
    return rvecs_[ridx](cidx);
}

template <typename ValueType>
Matrixs<ValueType> Matrixs<ValueType>::operator+(const Matrixs<ValueType> &smat) const {
    assert(nrows_ == smat.nrows_ && ncols_ == smat.ncols_);
    Matrixs<ValueType> result(*this);
    for (int idx = 0; idx < nrows_; ++idx) {
        result[idx] += smat[idx];
    }
    return result;
}

template <typename ValueType>
Matrixs<ValueType> Matrixs<ValueType>::operator-(const Matrixs<ValueType> &smat) const {
    assert(nrows_ == smat.nrows_ && ncols_ == smat.ncols_);
    Matrixs<ValueType> result(*this);
    for (int idx = 0; idx < nrows_; ++idx) {
        result[idx] -= smat[idx];
    }
    return result;
}

template <typename ValueType>
Matrixs<ValueType> Matrixs<ValueType>::operator*(const Matrixs<ValueType> &smat) const {
    assert(ncols_ == smat.nrows_);
    Matrixs<ValueType> result(nrows_, smat.ncols_);
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        ::std::vector<int> crcols = this->operator[](ridx).GetIndices();
        for (size_t ccidx = 0; ccidx < crcols.size(); ++ccidx) {
            int cidx = crcols[ccidx];
            result[ridx] += smat[cidx] * this->operator()(ridx, cidx);
        }
    }
    return result;
}

template <typename ValueType>
void Matrixs<ValueType>::Init(int nrows, int ncols) {
    nrows_ = nrows;
    ncols_ = ncols;
    for (int idx = 0; idx < nrows_; ++idx) {
        rvecs_.push_back(Vecxs<ValueType>(ncols_));
    }
}

template <typename ValueType>
Matrixs<ValueType> Matrixs<ValueType>::Transpose() const {
    Matrixs<ValueType> result(ncols_, nrows_);
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        ::std::vector<int> cidxs = this->operator[](ridx).GetIndices();
        for (size_t cidxidx = 0; cidxidx < cidxs.size(); ++cidxidx) {
            int cidx = cidxs[cidxidx];
            result(cidx, ridx) = TransposeTraits<ValueType>::mtranspose(this->operator()(ridx, cidx));
        }
    }
    return result;
}

template <typename ValueType>
Matrixs<ValueType> Matrixs<ValueType>::CholeskyLLT() const {
    assert(nrows_ == ncols_);
    Matrixs<ValueType> chol = Matrixs<ValueType>(*this);
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        ValueType &dii = chol(ridx, ridx);
        dii += lm_lambda_ * IdentityTraits<ValueType>::midentity(dii);
        ValueType diis = MsqrtTraits<ValueType>::msqrt(dii);
        ValueType diistransinv = InverseTraits<ValueType>::minv(TransposeTraits<ValueType>::mtranspose(diis));
        ::std::vector<int> sridxs = chol[ridx].GetIndicesget(ridx);
        for (size_t sridxidx = 0; sridxidx < sridxs.size(); ++sridxidx) {
            int sridx = sridxs[sridxidx];
            chol(sridx, ridx) *= diistransinv;
        }
        sridxs = chol[ridx].GetIndicesget(ridx + 1);
        for (size_t sridxidx = 0; sridxidx < sridxs.size(); ++sridxidx) {
            int sridx = sridxs[sridxidx];
            for (size_t scidxidx = 0; scidxidx < sridxs.size(); ++scidxidx) {
                int scidx = sridxs[scidxidx];
                chol(sridx, scidx) -= chol(sridx, ridx) * TransposeTraits<ValueType>::mtranspose(chol(scidx, ridx));
            }
        }
    }
    for (int ridx = 0; ridx < nrows_; ++ridx) {
        ::std::vector<int> sridxs = chol[ridx].GetIndicesget(ridx + 1);
        for (size_t cidxidx = 0; cidxidx < sridxs.size(); ++cidxidx) {
            int cidx = sridxs[cidxidx];
            chol.Delete(ridx, cidx);
        }
        SetLowtriTraits<ValueType>::msetlt(chol(ridx, ridx));
    }
    return chol;
}

template <typename ValueType>
Vecxs<ValueType> Matrixs<ValueType>::SolveWithlm(const Vecxs<ValueType> &b) const {
    assert(ncols_ == b.Length() && nrows_ == ncols_);
    Matrixs<ValueType> choll = CholeskyLLT();
    Vecxs<ValueType> y = ForwardSubst(choll, b);
    return BackwardSubst(choll, y);
}

template <typename ValueType>
Vecxs<ValueType> Matrixs<ValueType>::ForwardSubst(const Matrixs<ValueType> &choll, const Vecxs<ValueType> &b) const {
    Vecxs<ValueType> y(nrows_);
    for (int idx = 0; idx < nrows_; ++idx) {
        ValueType tmp = b[idx];
        ::std::vector<int> sindices = choll[idx].GetIndices();
        for (size_t sidxidx = 0; sidxidx < sindices.size(); ++sidxidx) {
            int sidx = sindices[sidxidx];
            if (sidx >= idx)
                continue;
            tmp -= choll(idx, sidx) * y[sidx];
        }
        y[idx] = InverseTraits<ValueType>::minv(choll(idx, idx)) * tmp;
    }
    return y;
}

template <typename ValueType>
Vecxs<ValueType> Matrixs<ValueType>::BackwardSubst(const Matrixs<ValueType> &choll, const Vecxs<ValueType> &y) const {
    Matrixs<ValueType> chollt = choll.Transpose();
    Vecxs<ValueType> x(nrows_);
    for (int idx = nrows_ - 1; idx >= 0; --idx) {
        ValueType tmp = y[idx];
        ::std::vector<int> sindices = chollt[idx].GetIndices();
        for (size_t sidxidx = 0; sidxidx < sindices.size(); ++sidxidx) {
            int sidx = sindices[sidxidx];
            if (sidx <= idx)
                continue;
            tmp -= chollt(idx, sidx) * x[sidx];
        }
        x[idx] = InverseTraits<ValueType>::minv(chollt(idx, idx)) * tmp;
    }
    return x;
}

template <typename ValueType>
void Matrixs<ValueType>::Print(const ::std::string &name) const {
    ::std::cout << "********Sparse matrix " << name << " " << nrows_ << " * " << ncols_ << std::endl;
    for (int idx = 0; idx < nrows_; ++idx) {
        std::cout << "****Row " << idx << " " << rvecs_[idx];
    }
    ::std::cout << "------------------------------------" << ::std::endl;
}

using Matrixsf = Matrixs<float>;
using Matrixsd = Matrixs<double>;
using Matrixsi = Matrixs<int>;

using Matrixsxf = Matrixs<Matrixxf>;
using Matrixsxd = Matrixs<Matrixxd>;
using Matrixsxi = Matrixs<Matrixxi>;

} // namespace matrix
} // namespace hitnlls