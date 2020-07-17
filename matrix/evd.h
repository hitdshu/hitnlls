#pragma once

namespace nlls {
template <typename E> 
class EVD {
public:
    using ElementType = typename E::ElementType;
    using ConstElementType = typename E::ConstElementType;
    static constexpr int NR = E::NR;
    static constexpr int NC = E::NC;
    using EvalType = Matrix<ElementType, NR, NC>;
    using VectorType = Matrix<ElementType, NR, 1>;
    NLLS_INLINE explicit EVD(const internal::Expression<E> &e) : e_(e.Cast()) { NLLS_CHECK_DIM_COMPATIBLE(NR, NC) Init(); }
    NLLS_INLINE explicit EVD(const MatrixBase<E> &e) : e_(e.Cast()) { NLLS_CHECK_DIM_COMPATIBLE(NR, NC) Init(); }
    void Compute() {
        if (ic_) {
            return;
        }
        const int max_iter_num = 10;
        const ElementType tol(1e-6);
        const int n = Rows();
        EvalType a(e_);
        u_.SetIdentity();
        ElementType ratio(1.0);
        ElementType norm = a.Norm();
        ElementType norm_sq = norm * norm;
        for (int iter = 0; iter < max_iter_num && ratio > tol ; ++iter) {
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    ElementType cost;
                    ElementType sint;
                    if (internal::CalcJacobRotCoeff<ElementType>(a(i, i), a(i, j), a(j, j), cost, sint)) {
                        for (int k = 0; k < n; ++k) {
                            ElementType aki = a(k, i);
                            a(k, i) = cost * a(k, i) - sint * a(k, j);
                            a(k, j) = cost * a(k, j) + sint * aki;
                            ElementType uki = u_(k, i);
                            u_(k, i) = cost * u_(k, i) - sint * u_(k, j);
                            u_(k, j) = cost * u_(k, j) + sint * uki;
                        }
                        for (int k = 0; k < n; ++k) {
                            ElementType aik = a(i, k);
                            a(i, k) = cost * a(i, k) - sint * a(j, k);
                            a(j, k) = cost * a(j, k) + sint * aik;
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
        for (int i = 0; i < n; ++i) {
            v_[i] = a(i, i);
        }
        a = u_;
        std::vector<int> order(n, 0);
        for (int i = 0; i < n; ++i) {
            order[i] = i;
        }
        std::sort(order.begin(), order.end(), [this](int i, int j) { return v_[i] < v_[j]; });
        VectorType vtmp = v_;
        for (int i = 0; i < n; ++i) {
            int ni = order[i];
            v_[i] = vtmp[ni];
            for (int j = 0; j < n; ++j) {
                u_(j, i) = a(j, ni);
            }
        }
        ic_ = true;
    }
    NLLS_INLINE const EvalType &U() { if (!ic_) Compute(); return u_; }
    NLLS_INLINE const VectorType &V() { if (!ic_) Compute(); return v_; }
    NLLS_INLINE int Rows() const { return e_.Rows(); }
    NLLS_INLINE int Cols() const { return e_.Cols(); }
    NLLS_INLINE int Size() const { return e_.Size(); }
private:
    void Init() {
        int m = e_.Rows();
        int n = e_.Cols();
        u_.Resize(m, n);
        v_.Resize(m, 1);
        ic_ = false;
    }

    bool ic_;
    const E &e_;
    EvalType u_;
    VectorType v_;
};
} // namespace nlls