#pragma once

namespace nlls {
namespace internal {
template <class E, bool B>
struct ExpressionTraits<Map<E, B>> {
    using ElementType = typename E::ElementType;
    using ConstElementType = typename E::ConstElementType;
    static constexpr int NR = E::NR;
    static constexpr int NC = E::NC;
    static constexpr bool A = false;
    static constexpr int C = E::C;
};
} // namespace internal
template <typename E>
class Map<E, false> : public MatrixBase<Map<E>> {
public:
    NLLS_MATRIX_EXTRACT_TYPES(Map)
    using BaseType = MatrixBase<Map<E>>;
    using LayoutType = typename internal::Layout<NR, NC>;
    NLLS_INLINE explicit Map(ElementType *data, int m, int n) : data_(data), layout_(m, n) {}
    NLLS_INLINE explicit Map(ElementType *data, int m) : data_(data), layout_(m, m) { NLLS_CHECK_DIM_ISVECTOR(NR, NC) }
    NLLS_INLINE explicit Map(ElementType *data) : data_(data) { NLLS_CHECK_DIM_ISFIXMAT(NR, NC) }
    NLLS_INLINE int Rows() const { return layout_.Rows(); }
    NLLS_INLINE int Cols() const { return layout_.Cols(); }
    NLLS_INLINE ConstElementType operator()(int r, int c) const { return data_[layout_.Offset(r, c)]; }
    NLLS_INLINE ElementType &operator()(int r, int c) { return data_[layout_.Offset(r, c)]; }
    NLLS_INLINE ConstElementType operator[](int i) const { NLLS_CHECK_DIM_ISVECTOR(NR, NC) return data_[i]; }
    NLLS_INLINE ElementType &operator[](int i) { NLLS_CHECK_DIM_ISVECTOR(NR, NC) return data_[i]; }
    template <class Derived>
    NLLS_INLINE Map &operator=(const internal::Expression<Derived> &e) { 
        NLLS_CHECK_DIM_COMPATIBLE(NR, Derived::NR)
        NLLS_CHECK_DIM_COMPATIBLE(NC, Derived::NC)
        BaseType::Assign(e);
        return *this;
    }
    NLLS_INLINE Map &operator=(const Map &m) { data_ = m.data_; layout_ = m.layout_; return *this; }
    NLLS_INLINE Map &operator=(Map &&m) { data_ = std::move(m.data_); layout_ = std::move(m.layout_); return *this; }
    NLLS_INLINE Map &operator=(const ElementType &v) { internal::ScalarExp<ElementType> e(v); BaseType::Assign(e); return *this; }
    template <class EO> Map &operator+=(const internal::Expression<EO> &e) { *this = *this + e; return *this; }
    template <class EO> Map &operator-=(const internal::Expression<EO> &e) { *this = *this - e; return *this; }
    template <class EO> Map &operator*=(const internal::Expression<EO> &e) { *this = *this * e; return *this; }
    NLLS_INLINE Map &operator+=(const ElementType &v) { *this = *this + v; return *this; }
    NLLS_INLINE Map &operator-=(const ElementType &v) { *this = *this - v; return *this; }
    NLLS_INLINE Map &operator*=(const ElementType &v) { *this = *this * v; return *this; }
    NLLS_INLINE Map &operator/=(const ElementType &v) { *this = *this / v; return *this; }
private:
    ElementType *data_;
    LayoutType layout_;
};
template <typename E>
class Map<E, true> : public internal::Expression<Map<E>> {
public:
    NLLS_MATRIX_EXTRACT_TYPES(Map)
    using BaseType = internal::Expression<Map<E>>;
    using LayoutType = typename internal::Layout<NR, NC>;
    NLLS_INLINE explicit Map(const ElementType *data, int m, int n) : data_(data), layout_(m, n) {}
    NLLS_INLINE explicit Map(const ElementType *data, int m) : data_(data), layout_(m, m) { NLLS_CHECK_DIM_ISVECTOR(NR, NC) }
    NLLS_INLINE explicit Map(const ElementType *data) : data_(data) { NLLS_CHECK_DIM_ISFIXMAT(NR, NC) }
    NLLS_INLINE int Rows() const { return layout_.Rows(); }
    NLLS_INLINE int Cols() const { return layout_.Cols(); }
    NLLS_INLINE ConstElementType operator()(int r, int c) const { return data_[layout_.Offset(r, c)]; }
    NLLS_INLINE ConstElementType operator[](int i) const { NLLS_CHECK_DIM_ISVECTOR(NR, NC) return data_[i]; }
    NLLS_INLINE const ElementType *Data() const { return data_; }
    NLLS_INLINE internal::BlockExp<Map> Block(int i, int j, int p, int q, int rs = 1, int cs = 1) const { return BaseType::Block(*this, i, j, p, q, rs, cs); }
    NLLS_INLINE internal::BlockExp<Map> Row(int i) const { return internal::BlockExp<Map>(*this, i, 0, 1, Cols()); }
    NLLS_INLINE internal::BlockExp<Map> Col(int j) const { return internal::BlockExp<Map>(*this, 0, j, Rows(), 1); }
private:
    const ElementType *data_;
    LayoutType layout_;
};
} // namespace nlls