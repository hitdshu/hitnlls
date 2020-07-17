#pragma once

namespace nlls {
namespace internal {
#define MULT_DEFAULT_COST (20)
#define MULT_THRE_COST (3)
#define NLLS_MATRIX_EXTRACT_TYPES(E) \
    using ConstElementType = typename internal::ExpressionTraits<E>::ConstElementType; \
    using ElementType = typename internal::ExpressionTraits<E>::ElementType; \
    static constexpr int NR = internal::ExpressionTraits<E>::NR; \
    static constexpr int NC = internal::ExpressionTraits<E>::NC; \
    static constexpr bool A = internal::ExpressionTraits<E>::A; \
    static constexpr int C = internal::ExpressionTraits<E>::C;
#define NLLS_CHECK_DIM_COMPATIBLE(N1, N2) \
    static_assert(N1 == N2 || N1 * N2 == 0, "DIMENSION NOT COMPATIBLE");
#define NLLS_CHECK_DIM_ISVECTOR(N1, N2) \
    static_assert(N1 == 1 || N2 == 1, "NOT A VECTOR");
#define NLLS_CHECK_DIM_ISFIXVEC(N1, N2) \
    static_assert((N1 == 1 || N2 == 1) && N1 * N2 != 0, "NOT A FIXED SIZE VECTOR");
#define NLLS_CHECK_DIM_ISVEC1D(N1, N2) \
    static_assert(N1 * N2 == 1, "NOT A VECTOR1D");
#define NLLS_CHECK_DIM_ISVEC2D(N1, N2) \
    static_assert(N1 * N2 == 2, "NOT A VECTOR2D");
#define NLLS_CHECK_DIM_ISVEC3D(N1, N2) \
    static_assert(N1 * N2 == 3, "NOT A VECTOR3D");
#define NLLS_CHECK_DIM_ISDYNVEC(N1, N2) \
    static_assert((N1 == 0 && N2 == 1) || (N1 == 1 && N2 == 0), "NOT A DYNAMIC SIZE VECTOR");
#define NLLS_CHECK_DIM_ISFIXMAT(N1, N2) \
    static_assert(N1 != 0 && N2 != 0, "NOT A FIXED SIZE MATRIX");
#define NLLS_CHECK_DIM_ISDYNMAT(N1, N2) \
    static_assert(N1 == 0 && N2 == 0, "NOT A DYNAMIC SIZE MATRIX");
#define NLLS_DEFINE_SINGLE_MATRIX(name, T, M, N) \
    using Matrix##name = Matrix<T, M, N>;
#define NLLS_DEFINE_SINGLE_VECTOR(name, T, N) \
    using Vector##name = Matrix<T, N, 1>; \
    using RowVector##name = Matrix<T, 1, N>;  
#define NLLS_DEFINE_MATRIX(name, M, N) \
    NLLS_DEFINE_SINGLE_MATRIX(name##f, float, M, N) \
    NLLS_DEFINE_SINGLE_MATRIX(name##d, double, M, N) \
    NLLS_DEFINE_SINGLE_MATRIX(name##j, Jetf, M, N) \
    NLLS_DEFINE_SINGLE_MATRIX(name##i, int, M, N)
#define NLLS_DEFINE_VECTOR(name, N) \
    NLLS_DEFINE_SINGLE_VECTOR(name##f, float, N) \
    NLLS_DEFINE_SINGLE_VECTOR(name##d, double, N) \
    NLLS_DEFINE_SINGLE_VECTOR(name##j, Jetf, N) \
    NLLS_DEFINE_SINGLE_VECTOR(name##i, int, N)
} // namespace internal
} // namespace nlls