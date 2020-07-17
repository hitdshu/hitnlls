#pragma once

namespace nlls {
template <class M> class MatrixBase;
template <typename T, int M, int N> class Matrix;
template <class E, bool C = std::is_const<E>::value> class Map;
template <typename E> class LUP;
template <typename E> class LLT;
template <typename E> class QR;
template <typename E> class EVD;
template <typename E> class SVD;
namespace internal {
enum { DYNAMIC = 0 };
template <class E>
struct ExpressionTraits {
    using ConstElementType = typename E::ConstElementType;
    using ElementType = typename E::ElementType;
    static constexpr int NR = E::NR;
    static constexpr int NC = E::NC;
    static constexpr bool A = E::A;
    static constexpr int C = E::C;
};
template <class E> struct ExpressionTraits<const E> : public ExpressionTraits<E> {};
template <class E> struct ExpressionTraits<E &> : public ExpressionTraits<E> {};
template <class E> class TransposeExp;
template <class E> class NormalizedExp;
template <class E> class InverseExp;
template <class E> class ElementwiseExp;
template <class E> class BlockExp;
template <class E> class MutableBlockExp;
template <class E> class NoAliasExp;
template <class E> class ScalarExp;
template <class E, class Op> class CoeffwiseUnaryOp;
template <class Op, class E1, class E2> class CoeffwiseBinaryOp;
template <class E1, class E2> class MultBinaryOp;
template <bool B> struct Bool2Type {};
template <class E> struct StorageTypeTraits { using Type = const E &; };
template <class T> struct StorageTypeTraits<ScalarExp<T>> { using Type = ScalarExp<T>; };
} // namespace internal
} // namespace nlls