#pragma once

namespace nlls {
namespace internal {
template <int M, int N>
class Layout {
public:
    NLLS_INLINE explicit Layout() {}
    NLLS_INLINE explicit Layout(int m, int n) {}
    NLLS_INLINE void Resize(int m, int n) {}
    NLLS_INLINE int Rows() const { return M; }
    NLLS_INLINE int Cols() const { return N; }
    NLLS_INLINE int Size() const { return M * N; }
    NLLS_INLINE int Offset(int r, int c) const { return r * Cols() + c;  }
};
template <int M>
class Layout<M, DYNAMIC> {
public:
    NLLS_INLINE explicit Layout() { n_ = 0; }
    NLLS_INLINE explicit Layout(int m, int n) { n_ = n;}
    NLLS_INLINE void Resize(int m, int n) { n_ = n; }
    NLLS_INLINE int Rows() const { return M; }
    NLLS_INLINE int Cols() const { return n_; }
    NLLS_INLINE int Size() const { return M * n_; }
    NLLS_INLINE int Offset(int r, int c) const { return r * Cols() + c;  }
private:
    int n_;
};
template <int N>
class Layout<DYNAMIC, N> {
public:
    NLLS_INLINE explicit Layout() { m_ = 0; }
    NLLS_INLINE explicit Layout(int m, int n) { m_ = m;}
    NLLS_INLINE void Resize(int m, int n) { m_ = m; }
    NLLS_INLINE int Rows() const { return m_; }
    NLLS_INLINE int Cols() const { return N; }
    NLLS_INLINE int Size() const { return m_ * N; }
    NLLS_INLINE int Offset(int r, int c) const { return r * Cols() + c;  }
private:
    int m_;
};
template <>
class Layout<DYNAMIC, DYNAMIC> {
public:
    NLLS_INLINE explicit Layout() { m_ = 0; n_ = 0; }
    NLLS_INLINE explicit Layout(int m, int n) { m_ = m; n_ = n; }
    NLLS_INLINE void Resize(int m, int n) { m_ = m; n_ = n; }
    NLLS_INLINE int Rows() const { return m_; }
    NLLS_INLINE int Cols() const { return n_; }
    NLLS_INLINE int Size() const { return m_ * n_; }
    NLLS_INLINE int Offset(int r, int c) const { return r * Cols() + c;  }
private:
    int m_;
    int n_;
};
} // namespace internal
} // namespace nlls