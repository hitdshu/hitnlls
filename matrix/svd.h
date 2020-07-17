#pragma once

namespace nlls {
template <typename E> 
class SVD {
public:
    using ElementType = typename E::ElementType;
    using ConstElementType = typename E::ConstElementType;
    static constexpr int NR = E::NR;
    static constexpr int NC = E::NC;
    using EvalType = Matrix<ElementType, NR, NC>;
    using UType = Matrix<ElementType, NR, (NR > NC) ? NC : NR>;
    using SType = Matrix<ElementType, (NR > NC) ? NC : NR, 1>;
    using VType = Matrix<ElementType, NC, NC>;
    NLLS_INLINE explicit SVD(const internal::Expression<E> &e) : e_(e.Cast()) { Init(); }
    NLLS_INLINE explicit SVD(const MatrixBase<E> &e) : e_(e.Cast()) { Init(); }
    void Compute() {
        if (ic_) {
            return;
        }
        const int max_iter_num = 10;
        const ElementType tol(1e-6);
        const int m = Rows();
        const int n = Cols();
        const int d = std::min(m, n);
        EvalType a(e_);
        v_.SetIdentity();
        ElementType ratio(1.0);
        ElementType norm = a.Norm();
        ElementType norm_sq = norm * norm;
        for (int iter = 0; iter < max_iter_num && ratio > tol ; ++iter) {
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    ElementType cost;
                    ElementType sint;
                    ElementType aii(0.0);
                    ElementType aij(0.0);
                    ElementType ajj(0.0);
                    for (int k = 0; k < m; ++k) {
                        aii += a(k, i) * a(k, i);
                        aij += a(k, i) * a(k, j);
                        ajj += a(k, j) * a(k, j);
                    }
                    if (internal::CalcJacobRotCoeff<ElementType>(aii, aij, ajj, cost, sint)) {
                        for (int k = 0; k < m; ++k) {
                            ElementType aki = a(k, i);
                            a(k, i) = cost * a(k, i) - sint * a(k, j);
                            a(k, j) = cost * a(k, j) + sint * aki;
                        }
                        for (int k = 0; k < n; ++k) {
                            ElementType vki = v_(k, i);
                            v_(k, i) = cost * v_(k, i) - sint * v_(k, j);
                            v_(k, j) = cost * v_(k, j) + sint * vki;
                        }
                    }
                }
            }
            ElementType diag_norm = 0;
            for (int i = 0; i < n; ++i) {
                diag_norm += nlls::pow(a(i, i), 2);
            }
            ElementType off_norm = nlls::sqrt(norm_sq - diag_norm);
            ratio = off_norm / norm;
        }
        Matrix<ElementType, NC, 1> os(n);
        for (int i = 0; i < n; ++i) {
            ElementType tmp(0.0);
            for (int k = 0; k < m; ++k) {
                tmp += a(k, i) * a(k, i);
            }
            os[i] = nlls::sqrt(tmp);
        }
        std::vector<int> order(n, 0);
        for (int i = 0; i < n; ++i) {
            order[i] = i;
        }
        std::sort(order.begin(), order.end(), [&](int i, int j) { return os[i] > os[j]; });
        VType ov = v_;
        for (int i = 0; i < n; ++i) {
            int ni = order[i];
            if (i < d) {
                for (int k = 0; k < m; ++k) {
                    u_(k, i) = a(k, ni) / os[ni];
                }
                s_[i] = os[ni];
            }
            for (int k = 0; k < n; ++k) {
                v_(k, i) = ov(k, ni);
            }
        }
        ic_ = true;
    }
    NLLS_INLINE const UType &U() { if (!ic_) Compute(); return u_; }
    NLLS_INLINE const SType &S() { if (!ic_) Compute(); return s_; }
    NLLS_INLINE const VType &V() { if (!ic_) Compute(); return v_; }
    NLLS_INLINE int Rows() const { return e_.Rows(); }
    NLLS_INLINE int Cols() const { return e_.Cols(); }
    NLLS_INLINE int Size() const { return e_.Size(); }
private:
    void Init() {
        int m = e_.Rows();
        int n = e_.Cols();
        u_.Resize(m, std::min(m, n));
        s_.Resize(std::min(m, n), 1);
        v_.Resize(n, n); 
        ic_ = false;
    }

    bool ic_;
    const E &e_;
    UType u_;
    SType s_;
    VType v_;
};
} // namespace nlls