#pragma once

namespace nlls {
namespace internal {
template <class M>
struct ExpressionTraits<MatrixBase<M>> : public ExpressionTraits<M> {
};
} // namespace internal
template <class Derived> 
class MatrixBase : public internal::Expression<MatrixBase<Derived>> {
public:
    using BaseType = internal::Expression<MatrixBase<Derived>>;
    NLLS_CRTP_REF
    NLLS_MATRIX_EXTRACT_TYPES(MatrixBase)
    NLLS_INLINE int Rows() const { return Cast().Rows(); }
    NLLS_INLINE int Cols() const { return Cast().Cols(); }
    NLLS_INLINE ConstElementType operator()(int r, int c) const { return Cast()(r, c); }
    NLLS_INLINE ElementType &operator()(int r, int c) { return Cast()(r, c); }
    using BaseType::Block;
    using BaseType::Row;
    using BaseType::Col;
    NLLS_INLINE internal::MutableBlockExp<Derived> Block(int i, int j, int p, int q, int rs = 1, int cs = 1) { return internal::MutableBlockExp<Derived>(Cast(), i, j, p, q, rs, cs); }
    NLLS_INLINE internal::MutableBlockExp<Derived> Row(int i) { return internal::MutableBlockExp<Derived>(Cast(), i, 0, 1, Cols()); }
    NLLS_INLINE internal::MutableBlockExp<Derived> Col(int j) { return internal::MutableBlockExp<Derived>(Cast(), 0, j, Rows(), 1); }
    NLLS_INLINE Derived &SetConstant(ConstElementType v) {
        for (int r = 0; r < Rows(); ++r) {
            for (int c = 0; c < Cols(); ++c) {
                this->operator()(r, c) = v;
            }
        }
        return Cast();
    }
    NLLS_INLINE Derived &SetZero() { return SetConstant(ElementType(0.0)); }
    Derived &SetIdentity() {
        SetZero();
        int nrows = Rows();
        int ncols = Cols();
        int n = (nrows > ncols) ? ncols : nrows;
        for (int i = 0; i < n; ++i) {
            this->operator()(i, i) = ConstElementType(1.0);
        }
        return Cast();
    }
    void TransformInPlace(const std::function<ElementType(ConstElementType)> &op) {
        for (int r = 0; r < Rows(); ++r) {
            for (int c = 0; c < Cols(); ++c) {
                this->operator()(r, c) = op(this->operator()(r, c));
            }
        }
    }
    void Normalize() {
        ElementType norm = BaseType::Norm();
        TransformInPlace([norm](ConstElementType v) { return ElementType(v / norm); });
    }
    template <class E, EnableIfType<E::A> * = nullptr>
    void Assign(const internal::Expression<E> &e) {
        Matrix<ElementType, NR, NC> tmp(Rows(), Cols());
        internal::Evaluator::Assign(tmp, e);
        internal::Evaluator::Assign(*this, tmp);
    }
    template <class E, EnableIfNotType<E::A> * = nullptr>
    void Assign(const internal::Expression<E> &e) {
        internal::Evaluator::Assign(*this, e);
    }
    NLLS_CRTP_DEC(MatrixBase)
};
namespace internal {
template <typename T, int M, int N>
struct ExpressionTraits<Matrix<T, M, N>> {
    using ConstElementType = const RemoveConstRefType<T> &;
    using ElementType = RemoveConstRefType<T>;
    static constexpr int NR = M;
    static constexpr int NC = N;
    static constexpr bool A = false;
    static constexpr int C = 1;
};
} // namespace internal
template <typename T, int M = internal::DYNAMIC, int N = internal::DYNAMIC>
class Matrix : public MatrixBase<Matrix<T, M, N>> {
public:
    NLLS_MATRIX_EXTRACT_TYPES(Matrix)
    using BaseType = MatrixBase<Matrix<T, M, N>>;
    using BlobType = typename internal::Blob<T, M * N>;
    using LayoutType = typename internal::Layout<M, N>;
    NLLS_INLINE explicit Matrix() {}
    NLLS_INLINE explicit Matrix(int m) : blob_(m), layout_(m, m) {}
    NLLS_INLINE explicit Matrix(int m, int n) : blob_(m * n), layout_(m, n) {}
    NLLS_INLINE explicit Matrix(int m, const T &v) : blob_(m, v), layout_(m, m) {}
    NLLS_INLINE explicit Matrix(int m, int n, const T &v) : blob_(m * n, v), layout_(m, n) {}
    NLLS_INLINE explicit Matrix(ConstElementType v) { NLLS_CHECK_DIM_ISVEC1D(M, N) blob_[0] = v; }
    NLLS_INLINE explicit Matrix(ConstElementType v1, ConstElementType v2) { NLLS_CHECK_DIM_ISVEC2D(M, N) blob_[0] = v1; blob_[1] = v2; }
    NLLS_INLINE explicit Matrix(ConstElementType v1, ConstElementType v2, ConstElementType v3) { NLLS_CHECK_DIM_ISVEC3D(M, N) blob_[0] = v1; blob_[1] = v2; blob_[2] = v3; }
    NLLS_INLINE explicit Matrix(const Matrix &m) : blob_(m.blob_), layout_(m.layout_) {}
    NLLS_INLINE explicit Matrix(Matrix &&m) : blob_(std::move(m.blob_)), layout_(std::move(m.layout_)) {}
    template <class Derived>
    NLLS_INLINE Matrix(const internal::Expression<Derived> &e) {
        NLLS_CHECK_DIM_COMPATIBLE(NR, Derived::NR)
        NLLS_CHECK_DIM_COMPATIBLE(NC, Derived::NC)
        Resize(e.Rows(), e.Cols());
        internal::Evaluator::Assign(*this, e);
    }
    explicit operator ElementType() const { NLLS_CHECK_DIM_ISVEC1D(M, N) return blob_[0]; }
    NLLS_INLINE int Rows() const { return layout_.Rows(); }
    NLLS_INLINE int Cols() const { return layout_.Cols(); }
    NLLS_INLINE ConstElementType operator()(int r, int c) const { return blob_[layout_.Offset(r, c)]; }
    NLLS_INLINE ElementType &operator()(int r, int c) { return blob_[layout_.Offset(r, c)]; }
    NLLS_INLINE ConstElementType operator[](int i) const { NLLS_CHECK_DIM_ISVECTOR(M, N) return blob_[i]; }
    NLLS_INLINE ElementType &operator[](int i) { NLLS_CHECK_DIM_ISVECTOR(M, N) return blob_[i]; }
    NLLS_INLINE void Resize(int m, int n) { blob_.Resize(m * n); layout_.Resize(m, n); }
    NLLS_INLINE void Resize(int m) { NLLS_CHECK_DIM_ISVECTOR(M, N) blob_.Resize(m); layout_.Resize(m, m); }
    NLLS_INLINE ElementType *Data() { return blob_.Data(); }
    NLLS_INLINE const ElementType *Data() const { return blob_.Data(); }
    template <typename TO>
    Matrix<TO, M, N> CastType() {
        Matrix<TO, M, N> result(Rows(), Cols());
        for (int r = 0; r < Rows(); ++r) {
            for (int c = 0; c < Cols(); ++c) {
                result(r, c) = TO(this->operator()(r, c));
            }
        }
        return result;
    }
    template <class Derived>
    NLLS_INLINE Matrix &operator=(const internal::Expression<Derived> &e) { 
        NLLS_CHECK_DIM_COMPATIBLE(NR, Derived::NR)
        NLLS_CHECK_DIM_COMPATIBLE(NC, Derived::NC)
        Resize(e.Rows(), e.Cols());
        BaseType::Assign(e);
        return *this;
    }
    NLLS_INLINE Matrix &operator=(const Matrix &m) { blob_ = m.blob_; layout_ = m.layout_; return *this; }
    NLLS_INLINE Matrix &operator=(Matrix &&m) { blob_ = std::move(m.blob_); layout_ = std::move(m.layout_); return *this; }
    NLLS_INLINE Matrix &operator=(const ElementType &v) { internal::ScalarExp<ElementType> e(v); BaseType::Assign(e); return *this; }
    template <class EO> Matrix &operator+=(const internal::Expression<EO> &e) { *this = *this + e; return *this; }
    template <class EO> Matrix &operator-=(const internal::Expression<EO> &e) { *this = *this - e; return *this; }
    template <class EO> Matrix &operator*=(const internal::Expression<EO> &e) { *this = *this * e; return *this; }
    NLLS_INLINE Matrix &operator+=(const ElementType &v) { *this = *this + v; return *this; }
    NLLS_INLINE Matrix &operator-=(const ElementType &v) { *this = *this - v; return *this; }
    NLLS_INLINE Matrix &operator*=(const ElementType &v) { *this = *this * v; return *this; }
    NLLS_INLINE Matrix &operator/=(const ElementType &v) { *this = *this / v; return *this; }
    static Matrix Identity() { NLLS_CHECK_DIM_ISFIXMAT(M, N) static const Matrix i(Matrix().SetIdentity()); return i; }
    static Matrix Identity(int m, int n) { NLLS_CHECK_DIM_ISDYNMAT(M, N) Matrix i(m, n); return i.SetIdentity(); }
    static Matrix Random() {
        NLLS_CHECK_DIM_ISFIXMAT(M, N)
        Matrix result;
        for (int r = 0; r < result.Rows(); ++r) {
            for (int c = 0; c < result.Cols(); ++c) {
                result(r, c) = internal::RandomNormal<T>();
            }
        }
        return result;
    }
    static Matrix Random(int m, int n) {
        NLLS_CHECK_DIM_ISDYNMAT(M, N)
        Matrix result(m, n);
        for (int r = 0; r < result.Rows(); ++r) {
            for (int c = 0; c < result.Cols(); ++c) {
                result(r, c) = internal::RandomNormal<T>();
            }
        }
        return result;
    }
private:
    class MatrixStream {
    public:
        NLLS_INLINE explicit MatrixStream(Matrix &m) : m_(m), i_(1) {}
        NLLS_INLINE MatrixStream &operator,(ConstElementType v) { m_.blob_[i_++] = v; return *this; }
        NLLS_INLINE MatrixStream &operator<<(ConstElementType v) { m_.blob_[i_++] = v; return *this; }
        NLLS_INLINE Matrix &Finished() const { return m_; }
    private:
        Matrix &m_;
        int i_;
    };
    friend class MatrixStream;
    BlobType blob_;
    LayoutType layout_;
public:
    NLLS_INLINE MatrixStream operator<<(ConstElementType v) { blob_[0] = v; return MatrixStream(*this); }
};
namespace internal {
template <class E>
struct ExpressionTraits<MutableBlockExp<E>> {
    using ConstElementType = typename E::ConstElementType;
    using ElementType = typename E::ElementType;
    static constexpr int NR = 0;
    static constexpr int NC = 0;
    static constexpr bool A = true;
    static constexpr int C = E::C + 1;
};
template <typename E>
class MutableBlockExp : public MatrixBase<MutableBlockExp<E>> {
public:
    using BaseType = MatrixBase<MutableBlockExp<E>>;
    NLLS_MATRIX_EXTRACT_TYPES(MutableBlockExp)
    NLLS_INLINE explicit MutableBlockExp(MatrixBase<E> &e, int i, int j, int p, int q, int rs = 1, int cs = 1) : e_(e.Cast()), i_(i), j_(j), p_(p), q_(q), rs_(rs), cs_(cs) {}
    NLLS_INLINE int Rows() const { return p_; }
    NLLS_INLINE int Cols() const { return q_; }
    NLLS_INLINE ConstElementType operator()(int r, int c) const { return e_(r * rs_ + i_, c * cs_ + j_); }
    NLLS_INLINE ElementType &operator()(int r, int c) { return e_(r * rs_ + i_, c * cs_ + j_); }
    template <class EO> NLLS_INLINE MutableBlockExp &operator=(const internal::Expression<EO> &e) { BaseType::Assign(e); return *this; }
    NLLS_INLINE MutableBlockExp &operator=(const ElementType &v) { internal::ScalarExp<ElementType> e(v); BaseType::Assign(e); return *this; }
    template <class EO> MutableBlockExp &operator+=(const internal::Expression<EO> &e) { *this = *this + e; return *this; }
    template <class EO> MutableBlockExp &operator-=(const internal::Expression<EO> &e) { *this = *this - e; return *this; }
    template <class EO> MutableBlockExp &operator*=(const internal::Expression<EO> &e) { *this = *this * e; return *this; }
    NLLS_INLINE MutableBlockExp &operator+=(const ElementType &v) { *this = *this + v; return *this; }
    NLLS_INLINE MutableBlockExp &operator-=(const ElementType &v) { *this = *this - v; return *this; }
    NLLS_INLINE MutableBlockExp &operator*=(const ElementType &v) { *this = *this * v; return *this; }
    NLLS_INLINE MutableBlockExp &operator/=(const ElementType &v) { *this = *this / v; return *this; }
private:
    E &e_;
    int i_;
    int j_;
    int p_;
    int q_;
    int rs_;
    int cs_;
};
} // namespace internal
NLLS_DEFINE_MATRIX(22, 2, 2)
NLLS_DEFINE_MATRIX(23, 2, 3)
NLLS_DEFINE_MATRIX(24, 2, 4)
NLLS_DEFINE_MATRIX(25, 2, 5)
NLLS_DEFINE_MATRIX(26, 2, 6)
NLLS_DEFINE_MATRIX(32, 3, 2)
NLLS_DEFINE_MATRIX(33, 3, 3)
NLLS_DEFINE_MATRIX(34, 3, 4)
NLLS_DEFINE_MATRIX(35, 3, 5)
NLLS_DEFINE_MATRIX(36, 3, 6)
NLLS_DEFINE_MATRIX(42, 4, 2)
NLLS_DEFINE_MATRIX(43, 4, 3)
NLLS_DEFINE_MATRIX(44, 4, 4)
NLLS_DEFINE_MATRIX(45, 4, 5)
NLLS_DEFINE_MATRIX(46, 4, 6)
NLLS_DEFINE_MATRIX(52, 5, 2)
NLLS_DEFINE_MATRIX(53, 5, 3)
NLLS_DEFINE_MATRIX(54, 5, 4)
NLLS_DEFINE_MATRIX(55, 5, 5)
NLLS_DEFINE_MATRIX(56, 5, 6)
NLLS_DEFINE_MATRIX(62, 6, 2)
NLLS_DEFINE_MATRIX(63, 6, 3)
NLLS_DEFINE_MATRIX(64, 6, 4)
NLLS_DEFINE_MATRIX(65, 6, 5)
NLLS_DEFINE_MATRIX(66, 6, 6)
NLLS_DEFINE_MATRIX(X, 0, 0)
NLLS_DEFINE_VECTOR(1, 1)
NLLS_DEFINE_VECTOR(2, 2)
NLLS_DEFINE_VECTOR(3, 3)
NLLS_DEFINE_VECTOR(4, 4)
NLLS_DEFINE_VECTOR(5, 5)
NLLS_DEFINE_VECTOR(6, 6)
NLLS_DEFINE_VECTOR(7, 7)
NLLS_DEFINE_VECTOR(8, 8)
NLLS_DEFINE_VECTOR(X, 0)
} // namespace nlls
