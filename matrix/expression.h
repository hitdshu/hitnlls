#pragma once

namespace hitnlls {
namespace matrix {

template <class Derived>
class Expression {
public:
    HITNLLS_CRTP
    HITNLLS_CRTP_ENSURE(Expression)

    HITNLLS_INLINE void Resize(int m, int n) { Cast().Resize(m, n); }

    HITNLLS_INLINE int Rows() const { return Cast().Rows(); }
    HITNLLS_INLINE int Cols() const { return Cast().Cols(); }
    HITNLLS_INLINE int Size() const { return Cast().Size(); }

    HITNLLS_EXPRESSION_DEFINE_ALIASING
    HITNLLS_INLINE Derived &NoAlias() { no_aliasing_ = true; return Cast(); }

    HITNLLS_INLINE bool Aliasing() const { return Derived::ALIASING && !no_aliasing_; }
    HITNLLS_INLINE TransposeExpression<Derived> Transpose() const { return TransposeExpression<Derived>(Cast()); }
    HITNLLS_INLINE InverseExpression<Derived> Inverse() const { return InverseExpression<Derived>(Cast()); }
    HITNLLS_INLINE BlockExpression<Derived> Block(int i, int j, int p, int q, int r = 1, int c = 1) const { return BlockExpression<Derived>(Cast(), i, j, p, q, r, c); }
    HITNLLS_INLINE BlockExpression<Derived> Row(int i) const { return BlockExpression<Derived>(Cast(), i, 0, 1, Cols()); }
    HITNLLS_INLINE BlockExpression<Derived> Col(int j) const { return BlockExpression<Derived>(Cast(), 0, j, Rows(), 1); }
    HITNLLS_INLINE SqrtExpression<Derived> Sqrt() const { return SqrtExpression<Derived>(Cast()); }
    HITNLLS_INLINE NormalizedExpression<Derived> Normalized() const { return NormalizedExpression<Derived>(Cast()); }
    HITNLLS_INLINE LUP<Derived> Lup() const { return LUP<Derived>(Cast()); }
    HITNLLS_INLINE LLT<Derived> Llt() const { return LLT<Derived>(Cast()); }
    HITNLLS_INLINE QR<Derived> Qr() const { return QR<Derived>(Cast()); }
    HITNLLS_INLINE EVD<Derived> Evd() const { return EVD<Derived>(Cast()); }
    HITNLLS_INLINE SVD<Derived> Svd() const { return SVD<Derived>(Cast()); }
};

template <typename Derived> 
std::ostream &operator<<(std::ostream &os, const Expression<Derived> &m) {
    for (int i = 0; i < m.Rows(); ++i) {
        for (int j = 0; j < m.Cols(); ++j) {
            os << std::setw(12) << m.Cast().At(i, j) << " ";
        }
        if (i != m.Rows() - 1) {
            os << std::endl;
        }
    }
    return os;
}

template <typename T>
class ScalarExpression : public Expression<ScalarExpression<T>> {
public:
    using EvalReturnType = T;
    using ElementType = T;
    using ShapeType = MatrixShape<DYNAMIC, DYNAMIC>;
    static const bool ALIASING = false;

    HITNLLS_INLINE explicit ScalarExpression(const T &v) : scalar_(v) { Resize(1, 1); }

    HITNLLS_INLINE operator ElementType() { return scalar_; }

    HITNLLS_INLINE void Resize(int m, int n) { shape_.Resize(m, n); }

    HITNLLS_INLINE const ElementType &At(int i, int j) const { return scalar_; }
    HITNLLS_INLINE int Rows() const { return shape_.Rows(); }
    HITNLLS_INLINE int Cols() const { return shape_.Cols(); }
    HITNLLS_INLINE int Size() const { return shape_.Size(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return shape_; }
    HITNLLS_EXPRESSION_DEFINE_NORM
    HITNLLS_EXPRESSION_DEFINE_TRACE
    HITNLLS_INLINE const ElementType &operator()(int i, int j) const { return At(i, j); }

private:
    T scalar_;
    ShapeType shape_;
};

template <typename E>
class TransposeExpression : public Expression<TransposeExpression<E>> {
public:
    using EvalReturnType = typename TransposeTypeTraits<E>::EvalReturnType;
    using ElementType = typename TransposeTypeTraits<E>::ElementType;
    using ShapeType = typename TransposeTypeTraits<E>::ShapeType;
    static const bool ALIASING = true;

    HITNLLS_INLINE explicit TransposeExpression(const Expression<E> &e) : e_(e.Cast()) { Resize(e.Cols(), e.Rows()); }

    HITNLLS_INLINE void Resize(int m, int n) { shape_.Resize(m, n); }

    HITNLLS_INLINE ElementType At(int i, int j) const { return e_.Cast().At(j, i); }
    HITNLLS_INLINE int Rows() const { return shape_.Rows(); }
    HITNLLS_INLINE int Cols() const { return shape_.Cols(); }
    HITNLLS_INLINE int Size() const { return shape_.Size(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return shape_; }
    HITNLLS_EXPRESSION_DEFINE_NORM
    HITNLLS_EXPRESSION_DEFINE_TRACE
    HITNLLS_INLINE ElementType operator()(int i, int j) const { return At(i, j); }

private:
    const E &e_;
    ShapeType shape_;
};

template <typename E>
class InverseExpression : public Expression<InverseExpression<E>> {
public:
    using EvalReturnType = typename ExpressionTypeTraits<E>::EvalReturnType;
    using ElementType = typename ExpressionTypeTraits<E>::ElementType;
    using ShapeType = typename ExpressionTypeTraits<E>::ShapeType;
    static const bool ALIASING = false;

    HITNLLS_INLINE explicit InverseExpression(const Expression<E> &e) : lup_(e.Cast()) { HITNLLS_CHECK_SQUARE_MATRIX(ShapeType) Resize(e.Cols(), e.Rows()); lup_.Compute(); }

    HITNLLS_INLINE void Resize(int m, int n) { shape_.Resize(m, n); lup_.Resize(m, n); }

    HITNLLS_INLINE const ElementType &At(int i, int j) const { return lup_.Inverse().At(i, j); }
    HITNLLS_INLINE int Rows() const { return shape_.Rows(); }
    HITNLLS_INLINE int Cols() const { return shape_.Cols(); }
    HITNLLS_INLINE int Size() const { return shape_.Size(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return shape_; }
    HITNLLS_EXPRESSION_DEFINE_NORM
    HITNLLS_EXPRESSION_DEFINE_TRACE
    HITNLLS_INLINE const ElementType &operator()(int i, int j) const { return At(i, j); }

private:
    LUP<E> lup_;
    ShapeType shape_;
};

template <typename E>
class SqrtExpression : public Expression<SqrtExpression<E>> {
public:
    using EvalReturnType = typename ExpressionTypeTraits<E>::EvalReturnType;
    using ElementType = typename ExpressionTypeTraits<E>::ElementType;
    using ShapeType = typename ExpressionTypeTraits<E>::ShapeType;
    static const bool ALIASING = false;

    HITNLLS_INLINE explicit SqrtExpression(const Expression<E> &e) : llt_(e.Cast()) { HITNLLS_CHECK_SQUARE_MATRIX(ShapeType) Resize(e.Rows(), e.Rows()); llt_.Compute(); }

    HITNLLS_INLINE void Resize(int m, int n) { shape_.Resize(m, n); llt_.Resize(m, n); }

    HITNLLS_INLINE const ElementType &At(int i, int j) const { return llt_.L().At(i, j); }
    HITNLLS_INLINE int Rows() const { return shape_.Rows(); }
    HITNLLS_INLINE int Cols() const { return shape_.Cols(); }
    HITNLLS_INLINE int Size() const { return shape_.Size(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return shape_; }
    HITNLLS_EXPRESSION_DEFINE_NORM
    HITNLLS_EXPRESSION_DEFINE_TRACE
    HITNLLS_INLINE const ElementType &operator()(int i, int j) const { return At(i, j); }

private:
    LLT<E> llt_;
    ShapeType shape_;
};

template <typename E>
class NormalizedExpression : public Expression<NormalizedExpression<E>> {
public:
    using EvalReturnType = typename ExpressionTypeTraits<E>::EvalReturnType;
    using ElementType = typename ExpressionTypeTraits<E>::ElementType;
    using ShapeType = typename ExpressionTypeTraits<E>::ShapeType;
    static const bool ALIASING = false;

    HITNLLS_INLINE explicit NormalizedExpression(const Expression<E> &e) : e_(e.Cast()) { Resize(e.Rows(), e.Cols()); norm_ = e_.Norm(); }

    HITNLLS_INLINE void Resize(int m, int n) { shape_.Resize(m, n); }

    HITNLLS_INLINE ElementType At(int i, int j) const { return e_.At(i, j) / norm_; }
    HITNLLS_INLINE int Rows() const { return shape_.Rows(); }
    HITNLLS_INLINE int Cols() const { return shape_.Cols(); }
    HITNLLS_INLINE int Size() const { return shape_.Size(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return shape_; }
    HITNLLS_EXPRESSION_DEFINE_NORM
    HITNLLS_EXPRESSION_DEFINE_TRACE
    HITNLLS_INLINE ElementType operator()(int i, int j) const { return At(i, j); }

private:
    const E &e_;
    ShapeType shape_;
    ElementType norm_;
};

} // namespace matrix
} // namespace hitnlls