#pragma once

namespace nlls {
namespace internal {
template <class E, class Op>
struct ExpressionTraits<CoeffwiseUnaryOp<E, Op>> {
    using ConstElementType = const typename E::ElementType;
    using ElementType = typename E::ElementType;
    static constexpr int NR = E::NR;
    static constexpr int NC = E::NC;
    static constexpr bool A = E::A;
    static constexpr int C = E::C;
};
template <class E, class Op>
class CoeffwiseUnaryOp : public Expression<CoeffwiseUnaryOp<E, Op>> {
public:
    NLLS_MATRIX_EXTRACT_TYPES(CoeffwiseUnaryOp)
    NLLS_INLINE explicit CoeffwiseUnaryOp(const Expression<E> &e) : e_(e.Cast()) {}
    NLLS_INLINE int Rows() const { return e_.Rows(); }
    NLLS_INLINE int Cols() const { return e_.Cols(); }
    NLLS_INLINE ElementType operator()(int r, int c) const { return Op::Apply(e_(r, c)); }
private:
    const E &e_;
};
template <class E, template <class E> class Op>
CoeffwiseUnaryOp<E, Op<E>> MakeCoeffwiseUnaryOp(const Expression<E> &e) {
    return CoeffwiseUnaryOp<E, Op<E>>(e);
}
template <class E>
struct CoeffwiseUnaryPlus {
    NLLS_INLINE static typename E::ElementType Apply(typename E::ConstElementType v) { return v; }
};
template <class E>
struct CoeffwiseUnaryMinus {
    NLLS_INLINE static typename E::ElementType Apply(typename E::ConstElementType v) { return -v; }
};
template <typename E>
CoeffwiseUnaryOp<E, CoeffwiseUnaryPlus<E>> operator+(const Expression<E> &e) {
    return MakeCoeffwiseUnaryOp<E, CoeffwiseUnaryPlus>(e);
}
template <typename E>
CoeffwiseUnaryOp<E, CoeffwiseUnaryMinus<E>> operator-(const Expression<E> &e) {
    return MakeCoeffwiseUnaryOp<E, CoeffwiseUnaryMinus>(e);
}
template <class Op, class E1, class E2>
struct ExpressionTraits<CoeffwiseBinaryOp<Op, E1, E2>> {
    using ElementType = typename std::common_type<typename E1::ElementType, typename E2::ElementType>::type;
    using ConstElementType = const ElementType;
    static constexpr int NR = (E1::NR > E2::NR) ? E1::NR : E2::NR;
    static constexpr int NC = (E1::NC > E2::NC) ? E1::NC : E2::NC;
    static constexpr bool A = E1::A || E2::A;
    static constexpr int C = E1::C + E2::C + 1;
};
template <class Op, class E1, class E2>
class CoeffwiseBinaryOp : public Expression<CoeffwiseBinaryOp<Op, E1, E2>> {
public:
    NLLS_MATRIX_EXTRACT_TYPES(CoeffwiseBinaryOp)
    NLLS_INLINE explicit CoeffwiseBinaryOp(const Expression<E1> &e1, const Expression<E2> &e2) : e1_(e1.Cast()), e2_(e2.Cast()) {}
    NLLS_INLINE int Rows() const { return e1_.Rows(); }
    NLLS_INLINE int Cols() const { return e1_.Cols(); }
    NLLS_INLINE ElementType operator()(int r, int c) const { return Op::Apply(e1_(r, c), e2_(r, c)); }
private:
    const E1 &e1_;
    typename StorageTypeTraits<E2>::Type e2_;
};
template <template <class E1, class E2> class Op, class E1, class E2>
CoeffwiseBinaryOp<Op<E1, E2>, E1, E2> MakeCoeffwiseBinaryOp(const Expression<E1> &e1, const Expression<E2> &e2) {
    return CoeffwiseBinaryOp<Op<E1, E2>, E1, E2>(e1, e2);
}
template <class E1, class E2>
struct CoeffwiseBinaryPlus {
    NLLS_INLINE static typename std::common_type<typename E1::ElementType, typename E2::ElementType>::type Apply(typename E1::ConstElementType &v1, typename E2::ConstElementType &v2) { return v1 + v2; }
};
template <class E1, class E2>
struct CoeffwiseBinaryMinus {
    NLLS_INLINE static typename std::common_type<typename E1::ElementType, typename E2::ElementType>::type Apply(typename E1::ConstElementType &v1, typename E2::ConstElementType &v2) { return v1 - v2; }
};
template <class E1, class E2>
struct CoeffwiseBinaryMult {
    NLLS_INLINE static typename std::common_type<typename E1::ElementType, typename E2::ElementType>::type Apply(typename E1::ConstElementType &v1, typename E2::ConstElementType &v2) { return v1 * v2; }
};
template <class E1, class E2>
struct CoeffwiseBinaryInvDivde {
    NLLS_INLINE static typename std::common_type<typename E1::ElementType, typename E2::ElementType>::type Apply(typename E1::ConstElementType &v1, typename E2::ConstElementType &v2) { return v2 / v1; }
};
template <typename E1, typename E2>
CoeffwiseBinaryOp<CoeffwiseBinaryPlus<E1, E2>, E1, E2> operator+(const Expression<E1> &e1, const Expression<E2> &e2) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryPlus, E1, E2>(e1, e2);
}
template <typename E>
CoeffwiseBinaryOp<CoeffwiseBinaryPlus<E, ScalarExp<typename E::ElementType>>, E, ScalarExp<typename E::ElementType>> operator+(const Expression<E> &e1, const typename E::ElementType &e2) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryPlus, E, ScalarExp<typename E::ElementType>>(e1, ScalarExp<typename E::ElementType>(e2));
}
template <typename E>
CoeffwiseBinaryOp<CoeffwiseBinaryPlus<E, ScalarExp<typename E::ElementType>>, E, ScalarExp<typename E::ElementType>> operator+(const typename E::ElementType &e2, const Expression<E> &e1) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryPlus, E, ScalarExp<typename E::ElementType>>(e1, ScalarExp<typename E::ElementType>(e2));
}
template <typename E1, typename E2>
CoeffwiseBinaryOp<CoeffwiseBinaryMinus<E1, E2>, E1, E2> operator-(const Expression<E1> &e1, const Expression<E2> &e2) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryMinus, E1, E2>(e1, e2);
}
template <typename E>
CoeffwiseBinaryOp<CoeffwiseBinaryMinus<E, ScalarExp<typename E::ElementType>>, E, ScalarExp<typename E::ElementType>> operator-(const Expression<E> &e1, const typename E::ElementType &e2) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryMinus, E, ScalarExp<typename E::ElementType>>(e1, ScalarExp<typename E::ElementType>(e2));
}
template <typename E>
CoeffwiseBinaryOp<CoeffwiseBinaryPlus<CoeffwiseUnaryOp<CoeffwiseUnaryMinus<E>, E>, ScalarExp<typename E::ElementType>>, CoeffwiseUnaryOp<CoeffwiseUnaryMinus<E>, E>, ScalarExp<typename E::ElementType>> operator-(const typename E::ElementType &e2, const Expression<E> &e1) { 
    return -e1 + e2;
}
template <typename E>
CoeffwiseBinaryOp<CoeffwiseBinaryMult<E, ScalarExp<typename E::ElementType>>, E, ScalarExp<typename E::ElementType>> operator*(const Expression<E> &e1, const typename E::ElementType &e2) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryMult, E, ScalarExp<typename E::ElementType>>(e1, ScalarExp<typename E::ElementType>(e2));
}
template <typename E>
CoeffwiseBinaryOp<CoeffwiseBinaryMult<E, ScalarExp<typename E::ElementType>>, E, ScalarExp<typename E::ElementType>> operator*(const typename E::ElementType &e2, const Expression<E> &e1) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryMult, E, ScalarExp<typename E::ElementType>>(e1, ScalarExp<typename E::ElementType>(e2));
}
template <typename E>
CoeffwiseBinaryOp<CoeffwiseBinaryMult<E, ScalarExp<typename E::ElementType>>, E, ScalarExp<typename E::ElementType>> operator/(const Expression<E> &e1, const typename E::ElementType &e2) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryMult, E, ScalarExp<typename E::ElementType>>(e1, ScalarExp<typename E::ElementType>(typename E::ElementType(1.0) / e2));
}
template <typename E>
CoeffwiseBinaryOp<CoeffwiseBinaryInvDivde<E, ScalarExp<typename E::ElementType>>, E, ScalarExp<typename E::ElementType>> operator/(const typename E::ElementType &e2, const Expression<E> &e1) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryInvDivde, E, ScalarExp<typename E::ElementType>>(e1, ScalarExp<typename E::ElementType>(e2));
}
template <class E1, class E2>
struct ExpressionTraits<MultBinaryOp<E1, E2>> {
    using ElementType = typename std::common_type<typename E1::ElementType, typename E2::ElementType>::type;
    using ConstElementType = const ElementType &;
    static constexpr int NR = E1::NR;
    static constexpr int NC = E2::NC;
    static constexpr bool A = false;
    static constexpr int C = 1;
};
template <class E1, class E2>
class MultBinaryOp : public Expression<MultBinaryOp<E1, E2>> {
public:
    NLLS_MATRIX_EXTRACT_TYPES(MultBinaryOp)
    using ReturnType = Matrix<ElementType, NR, NC>;
    NLLS_INLINE explicit MultBinaryOp(const Expression<E1> &e1, const Expression<E2> &e2) : e1_(e1.Cast()), e2_(e2.Cast()), result_(e1.Rows(), e2.Cols(), ElementType(0)) { 
        NLLS_CHECK_DIM_COMPATIBLE(E1::NC, E2::NR)
        Evaluate(Bool2Type<(E1::C <= MULT_THRE_COST)>(), Bool2Type<(E2::C <= MULT_THRE_COST)>());
    }
    NLLS_INLINE int Rows() const { return result_.Rows(); }
    NLLS_INLINE int Cols() const { return result_.Cols(); }
    NLLS_INLINE ConstElementType operator()(int r, int c) const { return result_(r, c); }
private:
    void Evaluate(Bool2Type<true>, Bool2Type<true>) {
        internal::Evaluator::Mult(result_, e1_, e2_);
    }
    void Evaluate(Bool2Type<true>, Bool2Type<false>) {
        Matrix<ElementType, E2::NR, E2::NC> e2_eval(e2_);
        internal::Evaluator::Mult(result_, e1_, e2_eval);
    }
    void Evaluate(Bool2Type<false>, Bool2Type<true>) {
        Matrix<ElementType, E1::NR, E1::NC> e1_eval(e1_);
        internal::Evaluator::Mult(result_, e1_eval, e2_);
    }
    void Evaluate(Bool2Type<false>, Bool2Type<false>) {
        Matrix<ElementType, E1::NR, E1::NC> e1_eval(e1_);
        Matrix<ElementType, E2::NR, E2::NC> e2_eval(e2_);
        internal::Evaluator::Mult(result_, e1_eval, e2_eval);
    }
    const E1 &e1_;
    const E2 &e2_;
    ReturnType result_;
};
template <typename E1, typename E2>
MultBinaryOp<E1, E2> operator*(const Expression<E1> &e1, const Expression<E2> &e2) {
    return MultBinaryOp<E1, E2>(e1, e2);
}
} // namespace internal
} // namespace nlls