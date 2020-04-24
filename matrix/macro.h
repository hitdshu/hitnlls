#pragma once

namespace hitnlls {
namespace matrix {

#define HITNLLS_INLINE inline

#define HITNLLS_CRTP \
HITNLLS_INLINE const Derived &Cast() const { return static_cast<const Derived &>(*this); } \
HITNLLS_INLINE Derived &Cast() { return static_cast<Derived &>(*this); }

#define HITNLLS_CRTP_ENSURE(T) \
private: \
    HITNLLS_INLINE T() {} \
    friend Derived; \
public:

#define HITNLLS_EXPRESSION_DEFINE_ALIASING \
protected: \
    bool no_aliasing_ = false; \
public:

#define HITNLLS_MATRIX_INPUT_OVERLOAD \
private: \
    int input_idx_; \
public: \
    ThisType &operator<<(const ElementType &v) { \
        input_idx_ = 0; \
        blob_[input_idx_++] = v; \
        return *this; \
    } \
    ThisType &operator,(const ElementType &v) { \
        blob_[input_idx_++] = v; \
        return *this; \
    }

#define HITNLLS_EXPRESSION_DEFINE_NORM \
    HITNLLS_INLINE ElementType Norm() const { \
        ElementType norm = 0; \
        for (int i = 0; i < Rows(); ++i) { \
            for (int j = 0; j < Cols(); ++j) { \
                norm += std::pow(At(i, j), 2); \
            } \
        } \
        norm = std::sqrt(norm); \
        return norm; \
    }

#define HITNLLS_EXPRESSION_DEFINE_TRACE \
    HITNLLS_INLINE ElementType Trace() const { \
        ElementType sum = 0; \
        for (int i = 0; i < Rows() && i < Cols(); ++i) { \
            sum += At(i, i); \
        } \
        return sum; \
    }

#define HITNLLS_CHECK_FIXSIZE_MATRIX(T) \
static_assert(ExpressionDimTraits<T>::NROWS != 0 && ExpressionDimTraits<T>::NCOLS != 0, "MATRIX IS NOT FIXE SIZED");

#define HITNLLS_CHECK_VECTOR(T) \
static_assert(ExpressionDimTraits<T>::NROWS == 1 || ExpressionDimTraits<T>::NCOLS == 1, "MATRIX IS NOT A VECTOR");

#define HITNLLS_CHECK_VECTOR2(T) \
static_assert(ExpressionDimTraits<T>::NROWS * ExpressionDimTraits<T>::NCOLS == 2, "MATRIX IS NOT A VECTOR2");

#define HITNLLS_CHECK_NOT_VECTOR2(T) \
static_assert(ExpressionDimTraits<T>::NROWS * ExpressionDimTraits<T>::NCOLS != 2, "MATRIX IS A VECTOR2");

#define HITNLLS_CHECK_VECTOR3(T) \
static_assert(ExpressionDimTraits<T>::NROWS * ExpressionDimTraits<T>::NCOLS == 3, "MATRIX IS NOT A VECTOR3");

#define HITNLLS_CHECK_SINGLE_DIM_COMPATIBLE(N1, N2) \
static_assert(N1 == N2 || N1 == 0 || N2 == 0, "MATRIX DIMENSIONS ARE NOT COMPATIBLE");

#define HITNLLS_CHECK_DIM_COMPARATIBLE(T1, T2) \
HITNLLS_CHECK_SINGLE_DIM_COMPATIBLE(ExpressionDimTraits<T1>::NROWS, ExpressionDimTraits<T2>::NROWS) \
HITNLLS_CHECK_SINGLE_DIM_COMPATIBLE(ExpressionDimTraits<T1>::NCOLS, ExpressionDimTraits<T2>::NCOLS)

#define HITNLLS_CHECK_MULT_COMPATIBLE(T1, T2) \
HITNLLS_CHECK_SINGLE_DIM_COMPATIBLE(ExpressionDimTraits<T1>::NCOLS, ExpressionDimTraits<T2>::NROWS)

#define HITNLLS_CHECK_SQUARE_MATRIX(T) \
HITNLLS_CHECK_SINGLE_DIM_COMPATIBLE(T::NROWS, T::NCOLS)

#define HITNLLS_DEFINE_SINGLE_MATRIX(name, T, M, N) \
using Matrix##name = Matrix<T, M, N>;

#define HITNLLS_DEFINE_SINGLE_VECTOR(name, T, N) \
using Vector##name = Matrix<T, N, 1>;  \
using RowVector##name = Matrix<T, 1, N>;  \

#define HITNLLS_DEFINE_MATRIX(name, M, N) \
HITNLLS_DEFINE_SINGLE_MATRIX(name##f, float, M, N) \
HITNLLS_DEFINE_SINGLE_MATRIX(name##d, double, M, N) \
HITNLLS_DEFINE_SINGLE_MATRIX(name##i, int, M, N)

#define HITNLLS_DEFINE_VECTOR(name, N) \
HITNLLS_DEFINE_SINGLE_VECTOR(name##f, float, N) \
HITNLLS_DEFINE_SINGLE_VECTOR(name##d, double, N) \
HITNLLS_DEFINE_SINGLE_VECTOR(name##i, int, N)

enum { DYNAMIC = 0 };

} // namespace matrix
} // namespace hitnlls