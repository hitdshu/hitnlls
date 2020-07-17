#pragma once

namespace nlls {
template <typename E> 
class LUP {
public:
    using ElementType = typename E::ElementType;
    using ConstElementType = typename E::ConstElementType;
    static constexpr int NR = E::NR;
    static constexpr int NC = E::NC;
    using EvalType = Matrix<ElementType, NR, NC>;
    using VectorType = Matrix<ElementType, NR, 1>;
    using PType = internal::Blob<int, NR>;
    NLLS_INLINE explicit LUP(const internal::Expression<E> &e) : e_(e.Cast()) { NLLS_CHECK_DIM_COMPATIBLE(NR, NC) Init(); }
    NLLS_INLINE explicit LUP(const MatrixBase<E> &e) : e_(e.Cast()) { NLLS_CHECK_DIM_COMPATIBLE(NR, NC) Init(); }
    void Compute() {
        lu_ = e_;
        for (int i = 0; i < Rows(); ++i) {
            p_[i] = i;
        }
        for (int i = 0; i < Rows(); ++i) {
            int p = i;
            ElementType pv = nlls::abs(lu_(i, i));
            for (int j = i + 1; j < Rows(); ++j) {
                ElementType tmp_pv = nlls::abs(lu_(j, i));
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
        VectorType b(Rows(), ElementType(0.0));
        for (int i = 0; i < Rows(); ++i) {
            b[i] = ElementType(1.0);
            VectorType x = Solve(b);
            for (int j = 0; j < Rows(); ++j) {
                inv_(j, i) = x[j];
            }
            b[i] = ElementType(0.0);
        }
    }
    template <typename VE>
    VectorType Solve(const internal::Expression<VE> &b) const {
        NLLS_CHECK_DIM_COMPATIBLE(NR, VE::NR)
        NLLS_CHECK_DIM_COMPATIBLE(1, VE::NC)
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
    NLLS_INLINE const EvalType &Inverse() const { return inv_; }
    NLLS_INLINE int Rows() const { return e_.Rows(); }
    NLLS_INLINE int Cols() const { return e_.Cols(); }
    NLLS_INLINE int Size() const { return e_.Size(); }
private:
    void Init() {
        int m = e_.Rows();
        int n = e_.Cols();
        lu_.Resize(m, n); 
        inv_.Resize(m, n); 
        p_.Resize(m);
        Compute();
    }

    const E &e_;
    PType p_;
    EvalType lu_;
    EvalType inv_;
};
} // namespace nlls