#pragma once

namespace hitnlls {
namespace matrix {

template <typename E> 
class LUP {
public:
    using EvalReturnType = typename ExpressionTypeTraits<E>::EvalReturnType;
    using ElementType = typename ExpressionTypeTraits<E>::ElementType;
    using ShapeType = typename ExpressionTypeTraits<E>::ShapeType;
    using VectorType = Matrix<ElementType, ShapeType::NROWS, 1>;
    using PType = Blob<int, ShapeType::NROWS>;

    HITNLLS_INLINE explicit LUP(const Expression<E> &e) : e_(e.Cast()) { HITNLLS_CHECK_SQUARE_MATRIX(ShapeType) Resize(e.Rows(), e.Cols()); }

    HITNLLS_INLINE void Resize(int m, int n) { shape_.Resize(m, n); lu_.Resize(m, n); inv_.Resize(m, n); p_.Resize(m); }
    void Compute() {
        lu_.NoAlias() = e_;
        for (int i = 0; i < Rows(); ++i) {
            p_[i] = i;
        }
        for (int i = 0; i < Rows(); ++i) {
            int p = i;
            ElementType pv = std::abs(lu_(i, i));
            for (int j = i + 1; j < Rows(); ++j) {
                ElementType tmp_pv = std::abs(lu_(j, i));
                if (tmp_pv > pv) {
                    p = j;
                    pv = tmp_pv;
                }
            }
            if (i != p) {
                std::swap(p_[i], p_[p]);
                for (int j = 0; j < Rows(); ++j) {
                    std::swap(lu_(i, j), lu_(p, j));
                }
            }
            for (int j = i + 1; j < Rows(); ++j) {
                lu_(j, i) /= lu_(i, i);
                for (int k = i + 1; k < Rows(); ++k) {
                    lu_(j, k) -= lu_(j, i) * lu_(i, k);
                }
            }
        }
        VectorType b(Rows(), ElementType(0));
        for (int i = 0; i < Rows(); ++i) {
            b[i] = 1;
            VectorType x = Solve(b);
            for (int j = 0; j < Rows(); ++j) {
                inv_(j, i) = x[j];
            }
            b[i] = 0;
        }
    }
    template <typename VE>
    VectorType Solve(const Expression<VE> &b) const {
        HITNLLS_CHECK_DIM_COMPARATIBLE(VectorType, VE)
        VectorType x(Rows());
        for (int i = 0; i < Rows(); ++i) {
            x[i] = b.Cast()(p_[i], 0);
            for (int j = 0; j < i; ++j) {
                x[i] -= lu_(i, j) * x[j];
            }
        }
        for (int i = Rows() - 1; i >= 0; --i) {
            for (int j = Rows() - 1; j > i; --j) {
                x[i] -= lu_(i, j) * x[j];
            }
            x[i] /= lu_(i, i);
        }
        return x;
    }
    HITNLLS_INLINE const EvalReturnType &Inverse() const { return inv_; }

    HITNLLS_INLINE int Rows() const { return shape_.Rows(); }
    HITNLLS_INLINE int Cols() const { return shape_.Cols(); }
    HITNLLS_INLINE int Size() const { return shape_.Size(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return shape_; }

private:
    EvalReturnType lu_;
    EvalReturnType inv_;
    PType p_;

    const E &e_;
    ShapeType shape_;
};

} // namespace matrix
} // namespace hitnlls