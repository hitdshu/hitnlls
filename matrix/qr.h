#pragma once

namespace hitnlls {
namespace matrix {

template <typename E> 
class QR {
public:
    using QType = typename ExpressionTypeTraits<E>::EvalReturnType;
    using QShape = typename ExpressionTypeTraits<E>::ShapeType;
    using ElementType = typename ExpressionTypeTraits<E>::ElementType;
    using RType = Matrix<ElementType, QShape::NCOLS, QShape::NCOLS>;
    using BType = Matrix<ElementType, QShape::NROWS, 1>;
    using XType = Matrix<ElementType, QShape::NCOLS, 1>;

    HITNLLS_INLINE explicit QR(const Expression<E> &e) : e_(e.Cast()) { Resize(e.Rows(), e.Cols()); }

    HITNLLS_INLINE void Resize(int m, int n) { shape_.Resize(m, n); q_.Resize(m, n); r_.Resize(n, n); }
    void Compute() {
        int m = Rows();
        int n = Cols();
        q_ = e_;
        r_ = ElementType(0);
        for (int j = 0; j < n; ++j) {
            ElementType rjj = 0;
            for (int i = 0; i < m; ++i) {
                rjj += std::pow(q_(i, j), 2);
            }
            rjj = std::sqrt(rjj);
            r_(j, j) = rjj;
            for (int i = 0; i < m; ++i) {
                q_(i, j) /= rjj;
            }
            for (int k = j + 1; k < n; ++k) {
                ElementType rjk = 0;
                for (int i = 0; i < m; ++i) {
                    rjk += q_(i, j) * q_(i, k);
                }
                r_(j, k) = rjk;
                for (int i = 0; i < m; ++i) {
                    q_(i, k) -= rjk * q_(i, j);
                }
            }
        }
    }
    template <typename VE>
    XType Solve(const Expression<VE> &b) const {
        HITNLLS_CHECK_DIM_COMPARATIBLE(BType, VE)
        int m = Rows();
        int n = Cols();
        XType x(n);
        x.NoAlias() = q_.Transpose() * b;
        for (int i = n - 1; i >= 0; --i) {
            for (int j = n - 1; j > i; --j) {
                x[i] -= r_(i, j) * x[j];
            }
            x[i] /= r_(i, i);
        }
        return x;
    }
    HITNLLS_INLINE const QType &Q() const { return q_; }
    HITNLLS_INLINE const RType &R() const { return r_; }

    HITNLLS_INLINE int Rows() const { return shape_.Rows(); }
    HITNLLS_INLINE int Cols() const { return shape_.Cols(); }
    HITNLLS_INLINE int Size() const { return shape_.Size(); }
    HITNLLS_INLINE const QShape &Shape() const { return shape_; }

private:
    QType q_;
    RType r_;

    const E &e_;
    QShape shape_;
};

} // namespace matrix
} // namespace hitnlls