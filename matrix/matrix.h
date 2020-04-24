#pragma once 

namespace hitnlls {
namespace matrix {

template <typename T, int M = DYNAMIC, int N = DYNAMIC>
class Matrix : public Expression<Matrix<T, M, N>> {
public:
    using ThisType = Matrix<T, M, N>;
    using BlobType = typename ExpressionTypeTraits<ThisType>::BlobType;
    using ShapeType = typename ExpressionTypeTraits<ThisType>::ShapeType;
    using EvalReturnType = typename ExpressionTypeTraits<ThisType>::EvalReturnType;
    using ElementType = typename ExpressionTypeTraits<ThisType>::ElementType;
    using LayoutType = typename ExpressionTypeTraits<ThisType>::LayoutType;
    static const bool ALIASING = false;

    HITNLLS_INLINE explicit Matrix() {}
    HITNLLS_INLINE explicit Matrix(int m) : blob_(m), layout_(m, m) { HITNLLS_CHECK_VECTOR(ThisType) }
    HITNLLS_INLINE explicit Matrix(int m, int n) : blob_(m * n), layout_(m, n) { HITNLLS_CHECK_NOT_VECTOR2(ThisType) }
    HITNLLS_INLINE explicit Matrix(int m, const T &v) : blob_(m, v), layout_(m, m) { HITNLLS_CHECK_VECTOR(ThisType) }
    HITNLLS_INLINE explicit Matrix(int m, int n, const T &v) : blob_(m * n, v), layout_(m, n) {}
    HITNLLS_INLINE explicit Matrix(const T &v1, const T &v2) { HITNLLS_CHECK_VECTOR2(ThisType) blob_[0] = v1; blob_[1] = v2; }
    HITNLLS_INLINE explicit Matrix(const T &v1, const T &v2, const T &v3) { HITNLLS_CHECK_VECTOR3(ThisType) blob_[0] = v1; blob_[1] = v2; blob_[2] = v3; }
    HITNLLS_INLINE explicit Matrix(Matrix &&m) : blob_(std::move(m.blob_)), layout_(std::move(m.layout_)) {}
    HITNLLS_INLINE explicit Matrix(const Matrix &m) : blob_(m.blob_), layout_(m.layout_) {}
    template <class E>
    HITNLLS_INLINE Matrix(const Expression<E> &e) {
        HITNLLS_CHECK_DIM_COMPARATIBLE(ThisType, E)
        Assign(e, false);
    }

    HITNLLS_INLINE static Matrix Identity() { Matrix tmp; tmp.SetIdentity(); return tmp; }
    HITNLLS_INLINE static Matrix Identity(int m, int n) { static_assert(M == 0 && N == 0, "MUST BE A DYNAMIC MATRIX"); Matrix tmp(m, n); tmp.SetIdentity(); return tmp; }

    HITNLLS_INLINE void SetZero() { blob_.Fill(0); }
    HITNLLS_INLINE void SetConstant(const ElementType &v) { blob_.Fill(v); }
    HITNLLS_INLINE void SetIdentity() {
        SetZero();
        int slen = std::min<int>(Rows(), Cols());
        for (int i = 0; i < slen; ++i) {
            At(i, i) = 1;
        }
    }
    HITNLLS_INLINE void Resize(int m, int n) { layout_.Resize(m, n); blob_.Resize(m * n); }
    template <class E>
    HITNLLS_INLINE void Assign(const Expression<E> &e, bool aliasing) {
        Resize(e.Rows(), e.Cols());
        if (aliasing) {
            ThisType tmp(e);
            Evaluator::Assign(*this, e.Cast());
        } else {
            Evaluator::Assign(*this, e.Cast());
        }
    }
    HITNLLS_INLINE BlockExpression<ThisType> Block(int i, int j, int p, int q, int r = 1, int c = 1) const { return BlockExpression<ThisType>(*this, i, j, p, q, r, c); }
    HITNLLS_INLINE MutableBlockExpression<ThisType> Block(int i, int j, int p, int q, int r = 1, int c = 1) { 
        return MutableBlockExpression<ThisType>(*this, i, j, p, q, r, c);
    }
    HITNLLS_INLINE BlockExpression<ThisType> Row(int i) const { return BlockExpression<ThisType>(*this, i, 0, 1, Cols()); }
    HITNLLS_INLINE BlockExpression<ThisType> Col(int j) const { return BlockExpression<ThisType>(*this, 0, j, Rows(), 1); }
    HITNLLS_INLINE MutableBlockExpression<ThisType> Row(int i) { return MutableBlockExpression<ThisType>(*this, i, 0, 1, Cols()); }
    HITNLLS_INLINE MutableBlockExpression<ThisType> Col(int j) { return MutableBlockExpression<ThisType>(*this, 0, j, Rows(), 1); }
    HITNLLS_INLINE ThisType &Normalize() {
        ElementType norm = Norm();
        for (int i = 0; i < Rows(); ++i) {
            for (int j = 0; j < Cols(); ++j) {
                At(i, j) /= norm;
            }
        }
        return *this;
    }

    HITNLLS_INLINE ElementType *Data() { return blob_.Data(); }
    HITNLLS_INLINE const ElementType *Data() const { return blob_.Data(); }
    HITNLLS_INLINE ElementType &At(int i, int j) { return blob_[Offset(i, j)]; }
    HITNLLS_INLINE const ElementType &At(int i, int j) const { return blob_[Offset(i, j)]; }
    HITNLLS_INLINE int Rows() const { return layout_.Rows(); }
    HITNLLS_INLINE int Cols() const { return layout_.Cols(); }
    HITNLLS_INLINE int Size() const { return layout_.Size(); }
    HITNLLS_INLINE int Offset(int i, int j) const { return layout_.Offset(i, j); }
    HITNLLS_INLINE const ShapeType &Shape() const { return layout_.Shape(); }
    HITNLLS_EXPRESSION_DEFINE_NORM
    HITNLLS_EXPRESSION_DEFINE_TRACE
    HITNLLS_INLINE const ElementType &operator()(int i, int j) const { return At(i, j); }
    HITNLLS_INLINE ElementType &operator()(int i, int j) { return At(i, j); }
    HITNLLS_INLINE ElementType &operator[](int i) { HITNLLS_CHECK_VECTOR(ThisType) return blob_[i]; }
    HITNLLS_INLINE const ElementType &operator[](int i) const { HITNLLS_CHECK_VECTOR(ThisType) return blob_[i]; }
    HITNLLS_INLINE ThisType &operator=(const ThisType &m) { blob_ = m.blob_; layout_ = m.layout_; return *this; }
    HITNLLS_INLINE ThisType &operator=(ThisType &&m) { blob_ = std::move(m.blob_); layout_ = std::move(m.layout_); return *this; }
    template <class E>
    HITNLLS_INLINE ThisType &operator=(const Expression<E> &e) { 
        HITNLLS_CHECK_DIM_COMPARATIBLE(ThisType, E)
        Assign(e, e.Aliasing() && !Expression<ThisType>::no_aliasing_);
        return *this;
    }
    template <class TO, EnableIfNotExpression<TO> * = nullptr>
    HITNLLS_INLINE ThisType &operator=(const TO &v) {
        ScalarExpression<TO> e(v); 
        Assign(e, e.Aliasing() && !Expression<ThisType>::no_aliasing_);
        return *this;
    }
    template <class E> ThisType &operator+=(const Expression<E> &e) { return *this = *this + e; }
    template <class E> ThisType &operator-=(const Expression<E> &e) { return *this = *this - e; }
    template <class E> ThisType &operator*=(const Expression<E> &e) { return *this = *this * e; }
    template <class TO, EnableIfNotExpression<TO> * = nullptr> ThisType &operator+=(const TO &v) { return *this = *this + v; }
    template <class TO, EnableIfNotExpression<TO> * = nullptr> ThisType &operator-=(const TO &v) { return *this = *this - v; }
    template <class TO, EnableIfNotExpression<TO> * = nullptr> ThisType &operator*=(const TO &v) { return *this = *this * v; }
    HITNLLS_MATRIX_INPUT_OVERLOAD

private:
    BlobType blob_;
    LayoutType layout_;
};

HITNLLS_DEFINE_MATRIX(22, 2, 2)
HITNLLS_DEFINE_MATRIX(23, 2, 3)
HITNLLS_DEFINE_MATRIX(33, 3, 3)
HITNLLS_DEFINE_MATRIX(34, 3, 4)
HITNLLS_DEFINE_MATRIX(44, 4, 4)
HITNLLS_DEFINE_MATRIX(X, 0, 0)
HITNLLS_DEFINE_VECTOR(2, 2)
HITNLLS_DEFINE_VECTOR(3, 3)
HITNLLS_DEFINE_VECTOR(4, 4)
HITNLLS_DEFINE_VECTOR(5, 5)
HITNLLS_DEFINE_VECTOR(6, 6)
HITNLLS_DEFINE_VECTOR(X, 0)

} // namespace matrix
} // namespace hitnlls