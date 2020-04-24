#pragma once

#include "matrix/dense.h"
#include "nlls/trait.h"
#include "nlls/sparse_array.h"

namespace hitnlls {

template <typename T>
class SparseMatrix {
public:
    SparseMatrix() { nrows_ = 0; ncols_ = 0; lm_lambda_ = 0.0; }
    SparseMatrix(int nrows, int ncols) { Init(nrows, ncols); lm_lambda_ = 0.0; }
    SparseMatrix(int nrows) { Init(nrows, nrows); lm_lambda_ = 0.0; }

    inline int Rows() const { return nrows_; }
    inline int Cols() const { return ncols_; }
    inline bool CountIndex(int r, int c) const { return rows_[r].CountIndex(c); }
    inline int NumElementsAtRow(int r) const { return rows_[r].NonzeroLen(); }
    inline std::vector<int> IndicesAtRow(int r) const { return rows_[r].GetIndices(); }
    SparseMatrix Tranpose() const {
        SparseMatrix tm(ncols_, nrows_);
        for (int r = 0; r < nrows_; ++r) {
            std::vector<int> indices = rows_[r].GetIndices();
            for (auto idx : indices) {
                tm(idx, r) = TransposeTraits<T>::Transpose(this->operator()(r, idx));
            }
        }
        return tm;
    }

    inline T &operator()(int r, int c) { return rows_[r][c]; }
    inline T operator()(int r, int c) const { return rows_[r][c]; }

    inline void SetLmLambda(double lambda) { lm_lambda_ = lambda; }
    SparseArray<T> SolveWithLm(const SparseArray<T> &b) const {
        SparseMatrix choll = CholeskyLLT();
        SparseArray<T> y = ForwardSubst(choll, b);
        return BackwardSubst(choll, y);
    }
    SparseMatrix CholeskyLLT() const {
        SparseMatrix chol(*this);
        for (int r = 0; r < nrows_; ++r) {
            T &dii = chol(r, r);
            dii += lm_lambda_ * IdentityTraits<T>::Identity(dii);
            T diis = SqrtTraits<T>::Sqrt(dii);
            T diistransinv = InverseTraits<T>::Inverse(TransposeTraits<T>::Transpose(diis));
            std::vector<int> idxs = chol.rows_[r].GetIndices();
            for (auto idx : idxs) {
                if (idx >= r) {
                    chol(idx, r) *= diistransinv;
                }
            }
            for (auto idx1 : idxs) {
                if (idx1 > r) {
                    for (auto idx2 : idxs) {
                        if (idx2 > r) {
                            if (CountIndex(idx1, idx2)) {
                                chol(idx1, idx2) -= chol(idx1, r) * TransposeTraits<T>::Transpose(chol(idx2, r));
                            } else {
                                chol(idx1, idx2) = -chol(idx1, r) * TransposeTraits<T>::Transpose(chol(idx2, r));
                            }
                        }
                    }
                }
            }
        }
        for (int r = 0; r < nrows_; ++r) {
            std::vector<int> idxs = chol.rows_[r].GetIndices();
            for (auto idx : idxs) {
                if (idx > r) {
                    chol.rows_[r].Delete(idx);
                }
            }
        }
        return chol;
    }
    SparseArray<T> ForwardSubst(const SparseMatrix &choll, const SparseArray<T> &b) const {
        SparseArray<T> y(nrows_);
        for (int r = 0; r < nrows_; ++r) {
            T tmp = b[r];
            std::vector<int> idxs = choll.rows_[r].GetIndices();
            for (auto idx : idxs) {
                if (idx >= r) {
                    continue;
                }
                tmp -= choll(r, idx) * y[idx];
            }
            y[r] = InverseTraits<T>::Inverse(choll(r, r)) * tmp;
        }
        return y;
    }
    SparseArray<T> BackwardSubst(const SparseMatrix &choll, const SparseArray<T> &y) const {
        SparseArray<T> x(nrows_);
        SparseMatrix chollt = choll.Tranpose();
        for (int r = nrows_ - 1; r >= 0; --r) {
            T tmp = y[r];
            std::vector<int> idxs = chollt.rows_[r].GetIndices();
            for (auto idx : idxs) {
                if (idx <= r) {
                    continue;
                }
                tmp -= chollt(r, idx) * x[idx];
            }
            x[r] = InverseTraits<T>::Inverse(chollt(r, r)) * tmp;
        }
        return x;
    }

    friend std::ostream &operator<<(std::ostream &os, const SparseMatrix &sm) {
        for (int i = 0; i < sm.Rows(); ++i) {
            for (int j = 0; j < sm.Cols(); ++j) {
                os << "(" << i << "," << j << ")" << std::setw(12) << sm(i, j) << " ";
            }
            if (i != sm.Rows()) {
                os << std::endl;
            }
        }
        return os;
    }

protected:
    void Init(int nrows, int ncols) {
        nrows_ = nrows;
        ncols_ = ncols;
        for (int i = 0; i < nrows_; ++i) {
            rows_.push_back(SparseArray<T>(ncols_));
        }
    }

private:
    int nrows_;
    int ncols_;
    double lm_lambda_;
    std::vector<SparseArray<T>> rows_;
};

using SparseMatrixf = SparseMatrix<float>;
using SparseMatrixd = SparseMatrix<double>;
using SparseMatrixi = SparseMatrix<int>;
using SparseBlockMatrix = SparseMatrix<matrix::MatrixXd>;

} // namespace hitnlls