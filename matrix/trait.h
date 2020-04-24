#pragma once

namespace hitnlls {
namespace matrix {

template <typename T, int N> class Blob;
template <typename E> class LUP;
template <typename E> class LLT;
template <typename E> class QR;
template <typename E> class EVD;
template <typename E> class SVD;
template <typename E> class ScalarExpression;
template <typename E> class TransposeExpression;
template <typename E> class InverseExpression;
template <typename E> class BlockExpression;
template <typename E> class MutableBlockExpression;
template <typename E> class SqrtExpression;
template <typename E> class NormalizedExpression;

template <class Layout> struct LayoutTraits;
template <int M, int N> class MatrixShape;
template <int M, int N> class LayoutCont;
template <int M, int N>
struct LayoutTraits<LayoutCont<M, N>> {
    static const int NROWS = M;
    static const int NCOLS = N;
    using ShapeType = MatrixShape<M, N>;
};
class LayoutView;
template <>
struct LayoutTraits<LayoutView> {
    static const int NROWS = DYNAMIC;
    static const int NCOLS = DYNAMIC;
    using ShapeType = MatrixShape<NROWS, NCOLS>;
};

template <class Derived> class Expression;
template <class T> using EnableIfExpression = typename std::enable_if<std::is_base_of<Expression<T>, T>::value>::type;
template <class T> using EnableIfNotExpression = typename std::enable_if<!std::is_base_of<Expression<T>, T>::value>::type;

template <class E> struct StorageTypeTraits { using Type = const E &; };
template <class T> struct StorageTypeTraits<ScalarExpression<T>> { using Type = ScalarExpression<T>; };

template <class E>
struct ExpressionTypeTraits {
    using EvalReturnType = typename E::EvalReturnType;
    using ElementType = typename E::ElementType;
    using ShapeType = typename E::ShapeType;
};
template <typename T, int M, int N> class Matrix;
template <typename T, int M, int N>
struct ExpressionTypeTraits<Matrix<T, M, N>> {
    using BlobType = Blob<T, M * N>;
    using LayoutType = LayoutCont<M, N>;
    using EvalReturnType = Matrix<T, M, N>;
    using ElementType = T;
    using ShapeType = typename LayoutTraits<LayoutType>::ShapeType;
};

template <typename E>
struct ExpressionDimTraits {
    static const int NROWS = ExpressionTypeTraits<E>::ShapeType::NROWS;
    static const int NCOLS = ExpressionTypeTraits<E>::ShapeType::NCOLS;
};

template <typename E1, typename E2, bool IS_MULT = false>
struct BinaryOpTypeTraits;
template <int N1, int N2>
struct CommonDimTraits {
    static const int N = (N1 > N2) ? N1 : N2;
};
template <typename E1, typename E2>
struct BinaryOpTypeTraits<E1, E2, false> {
private:
    using ShapeType1 = typename E1::ShapeType;
    using ShapeType2 = typename E2::ShapeType;
    static const int NROWS = CommonDimTraits<ShapeType1::NROWS, ShapeType2::NROWS>::N;
    static const int NCOLS = CommonDimTraits<ShapeType1::NCOLS, ShapeType2::NCOLS>::N;
public:
    using ElementType = typename std::common_type<typename E1::ElementType, typename E2::ElementType>::type;
    using ShapeType = MatrixShape<NROWS, NCOLS>;
    using EvalReturnType = Matrix<ElementType, NROWS, NCOLS>;
};
template <typename E1, typename E2>
struct BinaryOpTypeTraits<E1, E2, true> {
private:
    using ShapeType1 = typename E1::ShapeType;
    using ShapeType2 = typename E2::ShapeType;
    static const int NROWS = ShapeType1::NROWS;
    static const int NCOLS = ShapeType2::NCOLS;
public:
    using ElementType = typename std::common_type<typename E1::ElementType, typename E2::ElementType>::type;
    using ShapeType = MatrixShape<NROWS, NCOLS>;
    using EvalReturnType = Matrix<ElementType, NROWS, NCOLS>;
};

template <typename E>
struct TransposeTypeTraits {
private:
    using ShapeTypeE = typename E::ShapeType;
    static const int NROWS = ShapeTypeE::NCOLS;
    static const int NCOLS = ShapeTypeE::NROWS;
public:
    using ElementType = typename E::ElementType;
    using ShapeType = MatrixShape<NROWS, NCOLS>;
    using EvalReturnType = Matrix<ElementType, NROWS, NCOLS>;
};

template <int N1, int N2>
struct MinDimTraits {
    static const int N = (N1 > N2) ? N2 : N1;
};

template <typename E, bool = std::is_const<E>::value>
struct MapTypeTraits;
template <typename E>
struct MapTypeTraits<E, false> {
    using EType = E;
    using ElementType = typename ExpressionTypeTraits<EType>::ElementType;
    using ShapeType = typename ExpressionTypeTraits<EType>::ShapeType;
    using EvalReturnType = typename ExpressionTypeTraits<EType>::EvalReturnType;
    using LayoutType = typename ExpressionTypeTraits<EType>::LayoutType;
    using ElementPtrType = ElementType *;
    static const bool CONST_QUATIFIED = false;
};
template <typename E>
struct MapTypeTraits<E, true> {
    using EType = typename std::remove_cv<E>::type;
    using ElementType = typename ExpressionTypeTraits<EType>::ElementType;
    using ShapeType = typename ExpressionTypeTraits<EType>::ShapeType;
    using EvalReturnType = typename ExpressionTypeTraits<EType>::EvalReturnType;
    using LayoutType = typename ExpressionTypeTraits<EType>::LayoutType;
    using ElementPtrType = const ElementType *;
    static const bool CONST_QUATIFIED = true;
};
template <bool V> using EnableIfTrue = typename std::enable_if<V>::type;
template <bool V> using EnableIfFalse = typename std::enable_if<!V>::type;

} // namespace matrix
} // namespace hitnlls