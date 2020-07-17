#pragma once

namespace nlls {
namespace internal {
template <class Derived>
class Expression {
public:
    NLLS_CRTP_REF
    NLLS_MATRIX_EXTRACT_TYPES(Derived)
    NLLS_INLINE int Rows() const { return Cast().Rows(); }
    NLLS_INLINE int Cols() const { return Cast().Cols(); }
    NLLS_INLINE ConstElementType operator()(int r, int c) const { return Cast()(r, c); }
    NLLS_INLINE int Size() const { return Cast().Rows() * Cast().Cols(); }
    ElementType Norm() const {
        ElementType result(0.0);
        int nrows = Rows();
        int ncols = Cols();
        for (int r = 0; r < nrows; ++r) {
            for (int c = 0; c < ncols; ++c) {
                result += nlls::pow(this->operator()(r, c), 2);
            }
        }
        return nlls::sqrt(result);
    }
    ElementType Trace() const {
        ElementType result(0.0);
        int nrows = Rows();
        int ncols = Cols();
        int n = (nrows > ncols) ? ncols : nrows;
        for (int i = 0; i < n; ++i) {
            result += this->operator()(i, i);
        }
        return result;
    }
    template <class EO>
    ElementType Dot(const Expression<EO> &eo) const {
        NLLS_CHECK_DIM_COMPATIBLE(NR, EO::NR)
        NLLS_CHECK_DIM_COMPATIBLE(NC, EO::NC)
        ElementType result(0.0);
        for (int i = 0; i < Rows(); ++i) {
            for (int j = 0; j < Cols(); ++j) {
                result += this->operator()(i, j) * eo(i, j);
            }
        }
        return result;
    }
    NLLS_INLINE NormalizedExp<Derived> Normalized() const { return NormalizedExp<Derived>(Cast()); }
    NLLS_INLINE TransposeExp<Derived> Transpose() const { return TransposeExp<Derived>(Cast()); }
    NLLS_INLINE InverseExp<Derived> Inverse() const { return InverseExp<Derived>(Cast()); }
    NLLS_INLINE ElementwiseExp<Derived> Transform(const std::function<ElementType(ConstElementType)> &op) const { return ElementwiseExp<Derived>(Cast(), op); }
    NLLS_INLINE ElementwiseExp<Derived> Sqrt() const { return ElementwiseExp<Derived>(Cast(), [](ConstElementType v) { return ElementType(nlls::sqrt(v)); }); }
    NLLS_INLINE ElementwiseExp<Derived> Power(int n) const { return ElementwiseExp<Derived>(Cast(), [n](ConstElementType v) { return ElementType(nlls::pow(v, n)); }); }
    NLLS_INLINE BlockExp<Derived> Block(int i, int j, int p, int q, int rs = 1, int cs = 1) const { return BlockExp<Derived>(Cast(), i, j, p, q, rs, cs); }
    NLLS_INLINE BlockExp<Derived> Row(int i) const { return BlockExp<Derived>(Cast(), i, 0, 1, Cols()); }
    NLLS_INLINE BlockExp<Derived> Col(int j) const { return BlockExp<Derived>(Cast(), 0, j, Rows(), 1); }
    NLLS_INLINE NoAliasExp<Derived> NoAlias() const { return NoAliasExp<Derived>(Cast()); }
    LUP<Derived> Lup() const { return LUP<Derived>(Cast()); }
    LLT<Derived> Llt() const { return LLT<Derived>(Cast()); }
    QR<Derived> Qr() const { return QR<Derived>(Cast()); }
    EVD<Derived> Evd() const { return EVD<Derived>(Cast()); }
    SVD<Derived> Svd() const { return SVD<Derived>(Cast()); }
    NLLS_CRTP_DEC(Expression)
};
template <class E>
struct ExpressionTraits<TransposeExp<E>> {
    using ConstElementType = typename E::ConstElementType;
    using ElementType = typename E::ElementType;
    static constexpr int NR = E::NC;
    static constexpr int NC = E::NR;
    static constexpr bool A = true;
    static constexpr int C = E::C + 1;
};
template <typename E>
class TransposeExp : public Expression<TransposeExp<E>> {
public:
    NLLS_MATRIX_EXTRACT_TYPES(TransposeExp)
    NLLS_INLINE explicit TransposeExp(const Expression<E> &e) : e_(e.Cast()) {}
    NLLS_INLINE int Rows() const { return e_.Cols(); }
    NLLS_INLINE int Cols() const { return e_.Rows(); }
    NLLS_INLINE ConstElementType operator()(int r, int c) const { return e_(c, r); }
private:
    const E &e_;
};
template <class E>
struct ExpressionTraits<NormalizedExp<E>> {
    using ConstElementType = const typename E::ElementType;
    using ElementType = typename E::ElementType;
    static constexpr int NR = E::NR;
    static constexpr int NC = E::NC;
    static constexpr bool A = E::A;
    static constexpr int C = E::C + 1;
};
template <typename E>
class NormalizedExp : public Expression<NormalizedExp<E>> {
public:
    NLLS_MATRIX_EXTRACT_TYPES(NormalizedExp)
    NLLS_INLINE explicit NormalizedExp(const Expression<E> &e) : e_(e.Cast()) { norm_ = e_.Norm(); }
    NLLS_INLINE int Rows() const { return e_.Rows(); }
    NLLS_INLINE int Cols() const { return e_.Cols(); }
    NLLS_INLINE ElementType operator()(int r, int c) const { return e_(r, c) / norm_; }
private:
    const E &e_;
    ElementType norm_;
};
template <class E>
struct ExpressionTraits<InverseExp<E>> {
    using ConstElementType = const typename E::ElementType &;
    using ElementType = typename E::ElementType;
    static constexpr int NR = E::NR;
    static constexpr int NC = E::NC;
    static constexpr bool A = false;
    static constexpr int C = 1;
};
template <typename E>
class InverseExp : public Expression<InverseExp<E>> {
public:
    NLLS_MATRIX_EXTRACT_TYPES(InverseExp)
    NLLS_INLINE explicit InverseExp(const Expression<E> &e) : lup_(e.Cast()) { NLLS_CHECK_DIM_COMPATIBLE(NR, NC) }
    NLLS_INLINE int Rows() const { return lup_.Rows(); }
    NLLS_INLINE int Cols() const { return lup_.Cols(); }
    NLLS_INLINE ConstElementType operator()(int r, int c) const { return lup_.Inverse()(r, c); }
private:
    LUP<E> lup_;
};
template <class E>
struct ExpressionTraits<ElementwiseExp<E>> {
    using ConstElementType = const typename E::ElementType;
    using ElementType = typename E::ElementType;
    static constexpr int NR = E::NR;
    static constexpr int NC = E::NC;
    static constexpr bool A = E::A;
    static constexpr int C = E::C + 1;
};
template <typename E>
class ElementwiseExp : public Expression<ElementwiseExp<E>> {
public:
    NLLS_MATRIX_EXTRACT_TYPES(ElementwiseExp)
    using OpType = std::function<ElementType(ConstElementType)>;
    NLLS_INLINE explicit ElementwiseExp(const Expression<E> &e, const OpType &op) : e_(e.Cast()), op_(op) {}
    NLLS_INLINE int Rows() const { return e_.Rows(); }
    NLLS_INLINE int Cols() const { return e_.Cols(); }
    NLLS_INLINE ElementType operator()(int r, int c) const { return op_(e_(r, c)); }
private:
    const E &e_;
    OpType op_;
};
template <class E>
struct ExpressionTraits<BlockExp<E>> {
    using ElementType = typename E::ElementType;
    using ConstElementType = const ElementType;
    static constexpr int NR = 0;
    static constexpr int NC = 0;
    static constexpr bool A = true;
    static constexpr int C = E::C + 1;
};
template <typename E>
class BlockExp : public Expression<BlockExp<E>> {
public:
    NLLS_MATRIX_EXTRACT_TYPES(BlockExp)
    NLLS_INLINE explicit BlockExp(const Expression<E> &e, int i, int j, int p, int q, int rs = 1, int cs = 1) : e_(e.Cast()), i_(i), j_(j), p_(p), q_(q), rs_(rs), cs_(cs) {}
    NLLS_INLINE int Rows() const { return p_; }
    NLLS_INLINE int Cols() const { return q_; }
    NLLS_INLINE ConstElementType operator()(int r, int c) const { return e_(r * rs_ + i_, c * cs_ + j_); }
private:
    const E &e_;
    int i_;
    int j_;
    int p_;
    int q_;
    int rs_;
    int cs_;
};
template <class T>
struct ExpressionTraits<ScalarExp<T>> {
    using ElementType = RemoveConstRefType<T>;
    using ConstElementType = const ElementType &;
    static constexpr int NR = 0;
    static constexpr int NC = 0;
    static constexpr bool A = false;
    static constexpr int C = 1;
};
template <typename T>
class ScalarExp : public Expression<ScalarExp<T>> {
public:
    NLLS_MATRIX_EXTRACT_TYPES(ScalarExp)
    NLLS_INLINE explicit ScalarExp(const T &v) : v_(v) {}
    NLLS_INLINE int Rows() const { return 0; }
    NLLS_INLINE int Cols() const { return 0; }
    NLLS_INLINE ConstElementType operator()(int r, int c) const { return v_; }
private:
    ElementType v_;
};
template <class E>
struct ExpressionTraits<NoAliasExp<E>> {
    using ElementType = typename E::ElementType;
    using ConstElementType = typename E::ConstElementType;
    static constexpr int NR = E::NR;
    static constexpr int NC = E::NC;
    static constexpr bool A = false;
    static constexpr int C = E::C;
};
template <typename E>
class NoAliasExp : public Expression<NoAliasExp<E>> {
public:
    NLLS_MATRIX_EXTRACT_TYPES(NoAliasExp)
    NLLS_INLINE explicit NoAliasExp(const Expression<E> &e) : e_(e.Cast()) {}
    NLLS_INLINE int Rows() const { return e_.Rows(); }
    NLLS_INLINE int Cols() const { return e_.Cols(); }
    NLLS_INLINE ConstElementType operator()(int r, int c) const { return e_(r, c); }
private:
    const E e_;
};
template <typename Derived> 
std::ostream &operator<<(std::ostream &os, const Expression<Derived> &m) {
    int nrows = m.Rows();
    int ncols = m.Cols();
    for (int r = 0; r < nrows; ++r) {
        for (int c = 0; c < ncols; ++c) {
            os << std::right << std::setw(12) << m(r, c);
            if (c != ncols - 1) {
                std::cout << " ";
            }
        }
        if (r != nrows - 1) {
            os << std::endl;
        }
    }
    return os;
}
} // namespace internal
} // namespace nlls