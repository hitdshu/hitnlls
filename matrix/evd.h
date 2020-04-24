#pragma once

namespace hitnlls {
namespace matrix {

template <typename E> 
class EVD {
public:
    using EvalReturnType = typename ExpressionTypeTraits<E>::EvalReturnType;
    using ElementType = typename ExpressionTypeTraits<E>::ElementType;
    using ShapeType = typename ExpressionTypeTraits<E>::ShapeType;
    using VectorType = Matrix<ElementType, ShapeType::NROWS, 1>;

    HITNLLS_INLINE explicit EVD(const Expression<E> &e) : e_(e.Cast()) { HITNLLS_CHECK_SQUARE_MATRIX(ShapeType) Resize(e.Rows(), e.Rows()); }

    HITNLLS_INLINE void Resize(int m, int n) { shape_.Resize(m, n); u_.Resize(m, n); v_.Resize(m, 1); }
    void Compute() {
        const int max_iter_num = 10;
        const ElementType tol = 1e-7;
        const int n = Rows();
        EvalReturnType a(e_);
        u_.SetIdentity();
        double ratio = 1.0;
        ElementType norm = a.Norm();
        for (int iter = 0; iter < max_iter_num && ratio > tol ; ++iter) {
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    ElementType cost;
                    ElementType sint;
                    if (Helper::CalcJacobRotCoeff<ElementType>(a(i, i), a(i, j), a(j, j), cost, sint)) {
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
            ElementType off_norm = 0;
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    off_norm += std::pow(a(i, j), 2);
                }
            }
            off_norm = std::sqrt(off_norm);
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
        std::sort(order.begin(), order.end(), [&](int i, int j) { return v_[i] < v_[j]; });
        VectorType vtmp = v_;
        for (int i = 0; i < n; ++i) {
            int ni = order[i];
            v_[i] = vtmp[ni];
            for (int j = 0; j < n; ++j) {
                u_(j, i) = a(j, ni);
            }
        }
    }
    HITNLLS_INLINE const EvalReturnType &U() const { return u_; }
    HITNLLS_INLINE const VectorType &V() const { return v_; }

    HITNLLS_INLINE int Rows() const { return shape_.Rows(); }
    HITNLLS_INLINE int Cols() const { return shape_.Cols(); }
    HITNLLS_INLINE int Size() const { return shape_.Size(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return shape_; }

private:
    EvalReturnType u_;
    VectorType v_;

    const E &e_;
    ShapeType shape_;
};

} // namespace matrix
} // namespace hitnlls