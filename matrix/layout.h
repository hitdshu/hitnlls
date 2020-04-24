#pragma once

namespace hitnlls {
namespace matrix {

template <int M, int N>
class MatrixShape {
public:
    static const int NROWS = M;
    static const int NCOLS = N;

    HITNLLS_INLINE explicit MatrixShape() {}
    HITNLLS_INLINE explicit MatrixShape(int m, int n) {}

    HITNLLS_INLINE void Resize(int m, int n) {}

    HITNLLS_INLINE int Rows() const { return M; }
    HITNLLS_INLINE int Cols() const { return N; }
    HITNLLS_INLINE int Size() const { return M * N; }
};

template <int M>
class MatrixShape<M, DYNAMIC> {
public:
    static const int NROWS = M;
    static const int NCOLS = DYNAMIC;

    HITNLLS_INLINE explicit MatrixShape() { n_ = 0; }
    HITNLLS_INLINE explicit MatrixShape(int m, int n) { n_ = n;}

    HITNLLS_INLINE void Resize(int m, int n) { n_ = n; }

    HITNLLS_INLINE int Rows() const { return M; }
    HITNLLS_INLINE int Cols() const { return n_; }
    HITNLLS_INLINE int Size() const { return M * n_; }

private:
    int n_;
};

template <int N>
class MatrixShape<DYNAMIC, N> {
public:
    static const int NROWS = DYNAMIC;
    static const int NCOLS = N;

    HITNLLS_INLINE explicit MatrixShape() { m_ = 0; }
    HITNLLS_INLINE explicit MatrixShape(int m, int n) { m_ = m;}

    HITNLLS_INLINE void Resize(int m, int n) { m_ = m; }

    HITNLLS_INLINE int Rows() const { return m_; }
    HITNLLS_INLINE int Cols() const { return N; }
    HITNLLS_INLINE int Size() const { return m_ * N; }

private:
    int m_;
};

template <>
class MatrixShape<DYNAMIC, DYNAMIC> {
public:
    static const int NROWS = DYNAMIC;
    static const int NCOLS = DYNAMIC;

    HITNLLS_INLINE explicit MatrixShape() { m_ = 0; n_ = 0; }
    HITNLLS_INLINE explicit MatrixShape(int m, int n) { m_ = m; n_ = n; }

    HITNLLS_INLINE void Resize(int m, int n) { m_ = m; n_ = n; }

    HITNLLS_INLINE int Rows() const { return m_; }
    HITNLLS_INLINE int Cols() const { return n_; }
    HITNLLS_INLINE int Size() const { return m_ * n_; }

private:
    int m_;
    int n_;
};

template <class Derived>
class LayoutBase {
public:
    using ShapeType = typename LayoutTraits<Derived>::ShapeType;
    
    HITNLLS_CRTP
    HITNLLS_CRTP_ENSURE(LayoutBase)

    HITNLLS_INLINE void Resize(int m, int n) { Cast().Resize(); }

    HITNLLS_INLINE int Rows() const { return Cast().Rows(); }
    HITNLLS_INLINE int Cols() const { return Cast().Cols(); }
    HITNLLS_INLINE int Size() const { return Cast().Size(); }
    HITNLLS_INLINE int Offset(int i, int j) const { return Cast().Offset(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return Cast().Shape(); }
};

template <int M, int N>
class LayoutCont : public LayoutBase<LayoutCont<M, N>> {
public:
    using ThisType = LayoutCont<M, N>;
    using ShapeType = typename LayoutTraits<ThisType>::ShapeType;

    HITNLLS_INLINE explicit LayoutCont() : shape_() {}
    HITNLLS_INLINE explicit LayoutCont(int m, int n) : shape_(m, n) {}

    HITNLLS_INLINE void Resize(int m, int n) { shape_.Resize(m, n); }

    HITNLLS_INLINE int Rows() const { return shape_.Rows(); }
    HITNLLS_INLINE int Cols() const { return shape_.Cols(); }
    HITNLLS_INLINE int Size() const { return shape_.Size(); }
    HITNLLS_INLINE int Offset(int i, int j) const { return i * shape_.Cols() + j; }
    HITNLLS_INLINE const ShapeType &Shape() const { return shape_; }

private:
    ShapeType shape_;
};

class LayoutView : public LayoutBase<LayoutView> {
public:
    using ThisType = LayoutView;
    using ShapeType = MatrixShape<0, 0>;

    HITNLLS_INLINE explicit LayoutView(int i, int j, int p, int q, int r = 1, int c = 1) : shape_(p, q) {
        i_ = i; j_ = j; r_ = r; c_ = c;
    }

    HITNLLS_INLINE void Resize(int m, int n) {}

    HITNLLS_INLINE int Rows() const { return shape_.Rows(); }
    HITNLLS_INLINE int Cols() const { return shape_.Cols(); }
    HITNLLS_INLINE int Size() const { return shape_.Size(); }
    HITNLLS_INLINE int Offset(int i, int j) const { return 0; }
    HITNLLS_INLINE int OffsetRow(int i, int j) const { return i_ + i * r_; }
    HITNLLS_INLINE int OffsetCol(int i, int j) const { return j_ + j * c_; }
    HITNLLS_INLINE const ShapeType &Shape() const { return shape_; }

private:
    ShapeType shape_;
    int i_;
    int j_;
    int r_;
    int c_;
};

} // namespace matrix
} // namespace hitnlls