#pragma once

namespace nlls {
template <typename E> 
class LLT {
public:
    using ElementType = typename E::ElementType;
    using ConstElementType = typename E::ConstElementType;
    static constexpr int NR = E::NR;
    static constexpr int NC = E::NC;
    using EvalType = Matrix<ElementType, NR, NC>;
    using VectorType = Matrix<ElementType, NR, 1>;
    NLLS_INLINE explicit LLT(const internal::Expression<E> &e) : e_(e.Cast()) { NLLS_CHECK_DIM_COMPATIBLE(NR, NC) Init(); }
    NLLS_INLINE explicit LLT(const MatrixBase<E> &e) : e_(e.Cast()) { NLLS_CHECK_DIM_COMPATIBLE(NR, NC) Init(); }
    void Compute() {
        if (ic_) {
            return;
        }
        int size = Rows();
        l_ = e_;
        for (int i = 0; i < size; ++i) {
            ElementType dii = l_(i, i);
            ElementType diis = nlls::sqrt(dii);
            for (int j = i; j < size; ++j) {
                l_(j, i) /= diis;
            }
            for (int j = i + 1; j < size; ++j) {
                for (int k = i + 1; k <= j; ++k) {
                    l_(j, k) -= l_(j, i) * l_(k, i);
                }
            }
        }
        for (int i = 0; i < size; ++i) {
            for (int j = i + 1; j < size; ++j) {
                l_(i, j) = 0;
            }
        }
        ic_ = true;
    }
    template <typename VE>
    VectorType Solve(const internal::Expression<VE> &b) {
        NLLS_CHECK_DIM_COMPATIBLE(NR, VE::NR)
        NLLS_CHECK_DIM_COMPATIBLE(NC, VE::NC)
        if (!ic_) Compute();
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
    NLLS_INLINE const EvalType &L() { if (!ic_) Compute(); return l_; }
    NLLS_INLINE int Rows() const { return e_.Rows(); }
    NLLS_INLINE int Cols() const { return e_.Cols(); }
    NLLS_INLINE int Size() const { return e_.Size(); }
private:
    void Init() {
        l_.Resize(e_.Rows(), e_.Cols());
        ic_ = false;
    }

    bool ic_;
    const E &e_;
    EvalType l_;
};
} // namespace nlls