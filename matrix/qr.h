#pragma once

namespace nlls {
template <typename E> 
class QR {
public:
    using ElementType = typename E::ElementType;
    using ConstElementType = typename E::ConstElementType;
    static constexpr int NR = E::NR;
    static constexpr int NC = E::NC;
    using QType = Matrix<ElementType, NR, NC>;
    using RType = Matrix<ElementType, NC, NC>;
    using BType = Matrix<ElementType, NR, 1>;
    using XType = Matrix<ElementType, NC, 1>;
    NLLS_INLINE explicit QR(const internal::Expression<E> &e) : e_(e.Cast()) { Init(); }
    NLLS_INLINE explicit QR(const MatrixBase<E> &e) : e_(e.Cast()) { Init(); }
    void Compute() {
        if (ic_) {
            return;
        }
        int m = Rows();
        int n = Cols();
        q_ = e_;
        r_ = ElementType(0.0);
        for (int j = 0; j < n; ++j) {
            ElementType rjj(0.0);
            for (int i = 0; i < m; ++i) {
                rjj += nlls::pow(q_(i, j), 2);
            }
            rjj = nlls::sqrt(rjj);
            r_(j, j) = rjj;
            for (int i = 0; i < m; ++i) {
                q_(i, j) /= rjj;
            }
            for (int k = j + 1; k < n; ++k) {
                ElementType rjk(0.0);
                for (int i = 0; i < m; ++i) {
                    rjk += q_(i, j) * q_(i, k);
                }
                r_(j, k) = rjk;
                for (int i = 0; i < m; ++i) {
                    q_(i, k) -= rjk * q_(i, j);
                }
            }
        }
        ic_ = true;
    }
    template <typename VE>
    XType Solve(const internal::Expression<VE> &b) {
        NLLS_CHECK_DIM_COMPATIBLE(NR, VE::NR)
        NLLS_CHECK_DIM_COMPATIBLE(1, VE::NC)
        if (!ic_) Compute();
        int n = Cols();
        XType x(n);
        x = q_.Transpose().NoAlias() * b;
        for (int i = n - 1; i >= 0; --i) {
            for (int j = n - 1; j > i; --j) {
                x[i] -= r_(i, j) * x[j];
            }
            x[i] /= r_(i, i);
        }
        return x;
    }
    NLLS_INLINE const QType &Q() { if (!ic_) Compute(); return q_; }
    NLLS_INLINE const RType &R() { if (!ic_) Compute(); return r_; }
    NLLS_INLINE int Rows() const { return e_.Rows(); }
    NLLS_INLINE int Cols() const { return e_.Cols(); }
    NLLS_INLINE int Size() const { return e_.Size(); }
private:
    void Init() {
        int m = e_.Rows(); 
        int n = e_.Cols(); 
        q_.Resize(m, n); 
        r_.Resize(n, n); 
        ic_ = false;
    }

    bool ic_;
    const E &e_;
    QType q_;
    RType r_;
};
} // namespace nlls