#pragma once

namespace hitnlls {
namespace matrix {

template <typename E> 
class SVD {
public:
    using EvalReturnType = typename ExpressionTypeTraits<E>::EvalReturnType;
    using ElementType = typename ExpressionTypeTraits<E>::ElementType;
    using ShapeType = typename ExpressionTypeTraits<E>::ShapeType;
    using UType = Matrix<ElementType, ShapeType::NROWS, MinDimTraits<ShapeType::NROWS, ShapeType::NCOLS>::N>;
    using VType = Matrix<ElementType, ShapeType::NCOLS, ShapeType::NCOLS>;
    using SType = Matrix<ElementType, MinDimTraits<ShapeType::NROWS, ShapeType::NCOLS>::N, 1>;

    HITNLLS_INLINE explicit SVD(const Expression<E> &e) : e_(e.Cast()) { Resize(e.Rows(), e.Cols()); }

    HITNLLS_INLINE void Resize(int m, int n) { shape_.Resize(m, n); u_.Resize(m, std::min(m, n)); v_.Resize(n, n); s_.Resize(std::min(m, n), 1); }
    void Compute() {
        const int max_iter_num = 10;
        const ElementType tol = 1e-7;
        const int m = Rows();
        const int n = Cols();
        const int d = std::min(m, n);
        EvalReturnType a(e_);
        v_.SetIdentity();
        double ratio = 1.0;
        ElementType norm = a.Norm();
        for (int iter = 0; iter < max_iter_num && ratio > tol ; ++iter) {
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    ElementType cost;
                    ElementType sint;
                    ElementType aii = 0;
                    ElementType aij = 0;
                    ElementType ajj = 0;
                    for (int k = 0; k < m; ++k) {
                        aii += a(k, i) * a(k, i);
                        aij += a(k, i) * a(k, j);
                        ajj += a(k, j) * a(k, j);
                    }
                    if (Helper::CalcJacobRotCoeff<ElementType>(aii, aij, ajj, cost, sint)) {
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
            ElementType off_norm = 0;
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    ElementType tmp = 0;
                    for (int k = 0; k < m; ++k) {
                        tmp += a(k, i) * a(k, j);
                    }
                    off_norm += std::pow(tmp, 2);
                }
            }
            off_norm = std::sqrt(off_norm);
            ratio = off_norm / norm;
        }
        Matrix<ElementType, ShapeType::NCOLS, 1> os(n);
        for (int i = 0; i < n; ++i) {
            ElementType tmp = 0;
            for (int k = 0; k < m; ++k) {
                tmp += a(k, i) * a(k, i);
            }
            os[i] = std::sqrt(tmp);
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
    }
    HITNLLS_INLINE const UType &U() const { return u_; }
    HITNLLS_INLINE const SType &S() const { return s_; }
    HITNLLS_INLINE const VType &V() const { return v_; }

    HITNLLS_INLINE int Rows() const { return shape_.Rows(); }
    HITNLLS_INLINE int Cols() const { return shape_.Cols(); }
    HITNLLS_INLINE int Size() const { return shape_.Size(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return shape_; }

private:
    UType u_;
    SType s_;
    VType v_;

    const E &e_;
    ShapeType shape_;
};

} // namespace matrix
} // namespace hitnlls