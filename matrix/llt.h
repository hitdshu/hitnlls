#pragma once

namespace hitnlls {
namespace matrix {

template <typename E> 
class LLT {
public:
    using EvalReturnType = typename ExpressionTypeTraits<E>::EvalReturnType;
    using ElementType = typename ExpressionTypeTraits<E>::ElementType;
    using ShapeType = typename ExpressionTypeTraits<E>::ShapeType;
    using VectorType = Matrix<ElementType, ShapeType::NROWS, 1>;

    HITNLLS_INLINE explicit LLT(const Expression<E> &e) : e_(e.Cast()) { HITNLLS_CHECK_SQUARE_MATRIX(ShapeType) Resize(e.Rows(), e.Cols()); }

    HITNLLS_INLINE void Resize(int m, int n) { shape_.Resize(m, n); l_.Resize(m, n); }
    void Compute() {
        int size = Rows();
        l_.NoAlias() = e_;
        for (int i = 0; i < size; ++i) {
            ElementType dii = l_(i, i);
            ElementType diis = std::sqrt(dii);
            for (int j = i; j < size; ++j) {
                l_.At(j, i) /= diis;
            }
            for (int j = i + 1; j < size; ++j) {
                for (int k = i + 1; k <= j; ++k) {
                    l_.At(j, k) -= l_.At(j, i) * l_.At(k, i);
                }
            }
        }
        for (int i = 0; i < size; ++i) {
            for (int j = i + 1; j < size; ++j) {
                l_.At(i, j) = 0;
            }
        }
    }
    template <typename VE>
    VectorType Solve(const Expression<VE> &b) const {
        HITNLLS_CHECK_DIM_COMPARATIBLE(VectorType, VE)
        int size = Rows();
        VectorType x(size);
        for (int i = 0; i < size; ++i) {
            x[i] = b.Cast()(i, 0);
            for (int j = 0; j < i; ++j) {
                x[i] -= l_(i, j) * x[j];
            }
            x[i] /= l_(i, i);
        }
        for (int i = size - 1; i >= 0; --i) {
            for (int j = size - 1; j > i; --j) {
                x[i] -= l_(j, i) * x[j];
            }
            x[i] /= l_(i, i);
        }
        return x;
    }
    HITNLLS_INLINE const EvalReturnType &L() const { return l_; }

    HITNLLS_INLINE int Rows() const { return shape_.Rows(); }
    HITNLLS_INLINE int Cols() const { return shape_.Cols(); }
    HITNLLS_INLINE int Size() const { return shape_.Size(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return shape_; }

private:
    EvalReturnType l_;

    const E &e_;
    ShapeType shape_;
};

} // namespace matrix
} // namespace hitnlls