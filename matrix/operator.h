#pragma once

namespace hitnlls {
namespace matrix {

template <class Op, class E>
class CoeffwiseUnaryOp : public Expression<CoeffwiseUnaryOp<Op, E>> {
public:
    using EvalReturnType = typename E::EvalReturnType;
    using ElementType = typename E::ElementType;
    using ShapeType = typename E::ShapeType;
    static const bool ALIASING = E::ALIASING;

    HITNLLS_INLINE explicit CoeffwiseUnaryOp(const Expression<E> &e) : e_(e.Cast()) { Resize(e.Rows(), e.Cols()); }

    HITNLLS_INLINE void Resize(int m, int n) { shape_.Resize(m, n); }

    HITNLLS_INLINE ElementType At(int i, int j) const { return Op::Apply(e_, i, j); }
    HITNLLS_INLINE int Rows() const { return shape_.Rows(); }
    HITNLLS_INLINE int Cols() const { return shape_.Cols(); }
    HITNLLS_INLINE int Size() const { return shape_.Size(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return shape_; }
    HITNLLS_EXPRESSION_DEFINE_NORM
    HITNLLS_EXPRESSION_DEFINE_TRACE
    HITNLLS_INLINE ElementType operator()(int i, int j) const { return At(i, j); }

private:
    const E &e_;
    ShapeType shape_;
};

template <template <class E> class Op, class E>
HITNLLS_INLINE CoeffwiseUnaryOp<Op<E>, E> MakeCoeffwiseUnaryOp(const Expression<E> &e) {
    return CoeffwiseUnaryOp<Op<E>, E>(e);
}

template <class E>
struct CoeffwiseUnaryPlus {
    HITNLLS_INLINE static typename E::ElementType Apply(const Expression<E> &e, int i, int j) {
        return e.Cast()(i, j);
    }
};

template <class E>
struct CoeffwiseUnaryMinus {
    HITNLLS_INLINE static typename E::ElementType Apply(const Expression<E> &e, int i, int j) {
        return -e.Cast()(i, j);
    }
};

template <typename E>
HITNLLS_INLINE CoeffwiseUnaryOp<CoeffwiseUnaryPlus<E>, E> operator+(const Expression<E> &e) {
    return MakeCoeffwiseUnaryOp<CoeffwiseUnaryPlus, E>(e);
}

template <typename E>
HITNLLS_INLINE CoeffwiseUnaryOp<CoeffwiseUnaryMinus<E>, E> operator-(const Expression<E> &e) {
    return MakeCoeffwiseUnaryOp<CoeffwiseUnaryMinus, E>(e);
}

template <class Op, class E1, class E2>
class CoeffwiseBinaryOp : public Expression<CoeffwiseBinaryOp<Op, E1, E2>> {
public:
    using EvalReturnType = typename BinaryOpTypeTraits<E1, E2>::EvalReturnType;
    using ElementType = typename BinaryOpTypeTraits<E1, E2>::ElementType;
    using ShapeType = typename BinaryOpTypeTraits<E1, E2>::ShapeType;
    static const bool ALIASING = E1::ALIASING || E2::ALIASING;

    HITNLLS_INLINE explicit CoeffwiseBinaryOp(const Expression<E1> &e1, const Expression<E2> &e2) : e1_(e1.Cast()), e2_(e2.Cast()) { Resize(e1.Rows(), e1.Cols()); }

    HITNLLS_INLINE void Resize(int m, int n) { shape_.Resize(m, n); }

    HITNLLS_INLINE ElementType At(int i, int j) const { return Op::Apply(e1_, e2_, i, j); }
    HITNLLS_INLINE int Rows() const { return shape_.Rows(); }
    HITNLLS_INLINE int Cols() const { return shape_.Cols(); }
    HITNLLS_INLINE int Size() const { return shape_.Size(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return shape_; }
    HITNLLS_EXPRESSION_DEFINE_NORM
    HITNLLS_EXPRESSION_DEFINE_TRACE
    HITNLLS_INLINE ElementType operator()(int i, int j) const { return At(i, j); }

private:
    const E1 &e1_;
    typename StorageTypeTraits<E2>::Type e2_;
    ShapeType shape_;
};

template <template <class E1, class E2> class Op, class E1, class E2>
HITNLLS_INLINE CoeffwiseBinaryOp<Op<E1, E2>, E1, E2> MakeCoeffwiseBinaryOp(const Expression<E1> &e1, const Expression<E2> &e2) {
    return CoeffwiseBinaryOp<Op<E1, E2>, E1, E2>(e1, e2);
}

template <class E1, class E2>
struct CoeffwiseBinaryPlus {
    HITNLLS_INLINE static typename BinaryOpTypeTraits<E1, E2>::ElementType Apply(const Expression<E1> &e1, const Expression<E2> &e2, int i, int j) {
        return e1.Cast()(i, j) + e2.Cast()(i, j);
    }
};

template <class E1, class E2>
struct CoeffwiseBinaryMinus {
    HITNLLS_INLINE static typename BinaryOpTypeTraits<E1, E2>::ElementType Apply(const Expression<E1> &e1, const Expression<E2> &e2, int i, int j) {
        return e1.Cast()(i, j) - e2.Cast()(i, j);
    }
};

template <class E1, class E2>
struct CoeffwiseBinaryMult {
    HITNLLS_INLINE static typename BinaryOpTypeTraits<E1, E2>::ElementType Apply(const Expression<E1> &e1, const Expression<E2> &e2, int i, int j) {
        return e1.Cast()(i, j) * e2.Cast()(i, j);
    }
};

template <typename E1, typename E2>
HITNLLS_INLINE CoeffwiseBinaryOp<CoeffwiseBinaryPlus<E1, E2>, E1, E2> operator+(const Expression<E1> &e1, const Expression<E2> &e2) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryPlus, E1, E2>(e1, e2);
}

template <typename E, typename T, EnableIfNotExpression<T> * = nullptr>
HITNLLS_INLINE CoeffwiseBinaryOp<CoeffwiseBinaryPlus<E, ScalarExpression<T>>, E, ScalarExpression<T>> operator+(const Expression<E> &e1, const T &e2) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryPlus, E, ScalarExpression<T>>(e1, ScalarExpression<T>(e2));
}

template <typename E, typename T, EnableIfNotExpression<T> * = nullptr>
HITNLLS_INLINE CoeffwiseBinaryOp<CoeffwiseBinaryPlus<E, ScalarExpression<T>>, E, ScalarExpression<T>> operator+(const T &e2, const Expression<E> &e1) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryPlus, E, ScalarExpression<T>>(e1, ScalarExpression<T>(e2));
}

template <typename E1, typename E2>
HITNLLS_INLINE CoeffwiseBinaryOp<CoeffwiseBinaryMinus<E1, E2>, E1, E2> operator-(const Expression<E1> &e1, const Expression<E2> &e2) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryMinus, E1, E2>(e1, e2);
}

template <typename E, typename T, EnableIfNotExpression<T> * = nullptr>
HITNLLS_INLINE CoeffwiseBinaryOp<CoeffwiseBinaryMinus<E, ScalarExpression<T>>, E, ScalarExpression<T>> operator-(const Expression<E> &e1, const T &e2) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryMinus, E, ScalarExpression<T>>(e1, ScalarExpression<T>(e2));
}

template <typename E, typename T, EnableIfNotExpression<T> * = nullptr>
HITNLLS_INLINE CoeffwiseBinaryOp<CoeffwiseBinaryPlus<CoeffwiseUnaryOp<CoeffwiseUnaryMinus<E>, E>, ScalarExpression<T>>, CoeffwiseUnaryOp<CoeffwiseUnaryMinus<E>, E>, ScalarExpression<T>> operator-(const T &e2, const Expression<E> &e1) {
    return -e1 + e2;
}

template <typename E, typename T, EnableIfNotExpression<T> * = nullptr>
HITNLLS_INLINE CoeffwiseBinaryOp<CoeffwiseBinaryMult<E, ScalarExpression<T>>, E, ScalarExpression<T>> operator*(const Expression<E> &e1, const T &e2) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryMult, E, ScalarExpression<T>>(e1, ScalarExpression<T>(e2));
}

template <typename E, typename T, EnableIfNotExpression<T> * = nullptr>
HITNLLS_INLINE CoeffwiseBinaryOp<CoeffwiseBinaryMult<E, ScalarExpression<T>>, E, ScalarExpression<T>> operator*(const T &e2, const Expression<E> &e1) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryMult, E, ScalarExpression<T>>(e1, ScalarExpression<T>(e2));
}

template <typename E, typename T, EnableIfNotExpression<T> * = nullptr>
HITNLLS_INLINE CoeffwiseBinaryOp<CoeffwiseBinaryMult<E, ScalarExpression<T>>, E, ScalarExpression<T>> operator/(const Expression<E> &e1, const T &e2) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryMult, E, ScalarExpression<T>>(e1, ScalarExpression<T>(1.0 / e2));
}

template <typename E, typename T, EnableIfNotExpression<T> * = nullptr>
HITNLLS_INLINE CoeffwiseBinaryOp<CoeffwiseBinaryMult<E, ScalarExpression<T>>, E, ScalarExpression<T>> operator/(const T &e2, const Expression<E> &e1) {
    return MakeCoeffwiseBinaryOp<CoeffwiseBinaryMult, E, ScalarExpression<T>>(e1, ScalarExpression<T>(1.0 / e2));
}

template <class E1, class E2>
class MultBinaryOp : public Expression<MultBinaryOp<E1, E2>> {
public:
    using EvalReturnType = typename BinaryOpTypeTraits<E1, E2, true>::EvalReturnType;
    using ElementType = typename BinaryOpTypeTraits<E1, E2, true>::ElementType;
    using ShapeType = typename BinaryOpTypeTraits<E1, E2, true>::ShapeType;
    static const bool ALIASING = false;

    HITNLLS_INLINE explicit MultBinaryOp(const Expression<E1> &e1, const Expression<E2> &e2) : e1_(e1.Cast()), e2_(e2.Cast()), result_(e1.Rows(), e2.Cols(), 0) { 
        HITNLLS_CHECK_MULT_COMPATIBLE(E1, E2)
        Resize(e1.Rows(), e2.Cols());
        Eval();
    }

    HITNLLS_INLINE void Resize(int m, int n) { shape_.Resize(m, n); }
    HITNLLS_INLINE const EvalReturnType &Eval() {
        Evaluator::Mult(result_, e1_, e2_);
    }

    HITNLLS_INLINE const ElementType &At(int i, int j) const { return result_(i, j); }
    HITNLLS_INLINE int Rows() const { return shape_.Rows(); }
    HITNLLS_INLINE int Cols() const { return shape_.Cols(); }
    HITNLLS_INLINE int Size() const { return shape_.Size(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return shape_; }
    HITNLLS_EXPRESSION_DEFINE_NORM
    HITNLLS_EXPRESSION_DEFINE_TRACE
    HITNLLS_INLINE const ElementType &operator()(int i, int j) const { return At(i, j); }

private:
    const E1 &e1_;
    const E2 &e2_;
    ShapeType shape_;
    EvalReturnType result_;
};

template <typename E1, typename E2>
HITNLLS_INLINE MultBinaryOp<E1, E2> operator*(const Expression<E1> &e1, const Expression<E2> &e2) {
    return MultBinaryOp<E1, E2>(e1, e2);
}

} // namespace matrix
} // namespace hitnlls