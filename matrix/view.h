#pragma once

namespace hitnlls {
namespace matrix {

template <typename E>
class BlockExpression : public Expression<BlockExpression<E>> {
public:
    using ElementType = typename ExpressionTypeTraits<E>::ElementType;
    using ShapeType = MatrixShape<0, 0>;
    using EvalReturnType = Matrix<ElementType, 0, 0>;
    using LayoutType = LayoutView;
    static const bool ALIASING = true;

    HITNLLS_INLINE explicit BlockExpression(const Expression<E> &e, int i, int j, int p, int q, int r = 1, int c = 1) : e_(e.Cast()), layout_(i, j, p, q, r, c) {}

    HITNLLS_INLINE void Resize(int m, int n) { layout_.Resize(m, n); }

    HITNLLS_INLINE ElementType At(int i, int j) const { return e_.At(OffsetRow(i, j), OffsetCol(i, j)); }
    HITNLLS_INLINE int Rows() const { return layout_.Rows(); }
    HITNLLS_INLINE int Cols() const { return layout_.Cols(); }
    HITNLLS_INLINE int Size() const { return layout_.Size(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return layout_.Shape(); }
    HITNLLS_INLINE int OffsetRow(int i, int j) const { return layout_.OffsetRow(i, j); }
    HITNLLS_INLINE int OffsetCol(int i, int j) const { return layout_.OffsetCol(i, j); }

    HITNLLS_INLINE ElementType operator()(int i, int j) const { return At(i, j); }
    HITNLLS_EXPRESSION_DEFINE_NORM
    HITNLLS_EXPRESSION_DEFINE_TRACE
private:
    const E &e_;
    LayoutType layout_;
};

template <typename E>
class MutableBlockExpression : public Expression<MutableBlockExpression<E>> {
public:
    using ThisType = MutableBlockExpression<E>;
    using ElementType = typename ExpressionTypeTraits<E>::ElementType;
    using ShapeType = MatrixShape<0, 0>;
    using EvalReturnType = Matrix<ElementType, 0, 0>;
    using LayoutType = LayoutView;
    static const bool ALIASING = true;

    HITNLLS_INLINE explicit MutableBlockExpression(Expression<E> &e, int i, int j, int p, int q, int r = 1, int c = 1) : e_(e.Cast()), layout_(i, j, p, q, r, c) {}

    HITNLLS_INLINE void Resize(int m, int n) { layout_.Resize(m, n); }
    HITNLLS_INLINE BlockExpression<ThisType> Block(int i, int j, int p, int q, int r = 1, int c = 1) const { return BlockExpression<ThisType>(*this, i, j, p, q, r, c); }
    HITNLLS_INLINE MutableBlockExpression<ThisType> Block(int i, int j, int p, int q, int r = 1, int c = 1) { 
        return MutableBlockExpression<ThisType>(*this, i, j, p, q, r, c);
    }
    HITNLLS_INLINE BlockExpression<ThisType> Row(int i) const { return BlockExpression<ThisType>(*this, i, 0, 1, Cols()); }
    HITNLLS_INLINE BlockExpression<ThisType> Col(int j) const { return BlockExpression<ThisType>(*this, 0, j, Rows(), 1); }
    HITNLLS_INLINE MutableBlockExpression<ThisType> Row(int i) { return MutableBlockExpression<ThisType>(*this, i, 0, 1, Cols()); }
    HITNLLS_INLINE MutableBlockExpression<ThisType> Col(int j) { return MutableBlockExpression<ThisType>(*this, 0, j, Rows(), 1); }
    template <class EO>
    HITNLLS_INLINE void Assign(const Expression<EO> &e, bool aliasing) {
        Resize(e.Rows(), e.Cols());
        if (aliasing) {
            EvalReturnType tmp(e);
            Evaluator::Assign(*this, e.Cast());
        } else {
            Evaluator::Assign(*this, e.Cast());
        }
    }

    HITNLLS_INLINE const ElementType &At(int i, int j) const { return e_.At(OffsetRow(i, j), OffsetCol(i, j)); }
    HITNLLS_INLINE ElementType &At(int i, int j) { return e_.At(OffsetRow(i, j), OffsetCol(i, j)); }
    HITNLLS_INLINE int Rows() const { return layout_.Rows(); }
    HITNLLS_INLINE int Cols() const { return layout_.Cols(); }
    HITNLLS_INLINE int Size() const { return layout_.Size(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return layout_.Shape(); }
    HITNLLS_INLINE int OffsetRow(int i, int j) const { return layout_.OffsetRow(i, j); }
    HITNLLS_INLINE int OffsetCol(int i, int j) const { return layout_.OffsetCol(i, j); }

    HITNLLS_INLINE const ElementType &operator()(int i, int j) const { return At(i, j); }
    HITNLLS_INLINE ElementType &operator()(int i, int j) { return At(i, j); }

    template <class EO>
    HITNLLS_INLINE ThisType &operator=(const Expression<EO> &e) {
        Assign(e, e.Aliasing() && !Expression<ThisType>::no_aliasing_);
        return *this;
    }
    template <class T, EnableIfNotExpression<T> * = nullptr>
    HITNLLS_INLINE ThisType &operator=(const T &v) {
        ScalarExpression<T> e(v); 
        Assign(e, e.Aliasing() && !Expression<ThisType>::no_aliasing_);
        return *this;
    }
    template <class EO> ThisType &operator+=(const Expression<EO> &e) { return *this = *this + e; }
    template <class EO> ThisType &operator-=(const Expression<EO> &e) { return *this = *this - e; }
    template <class EO> ThisType &operator*=(const Expression<EO> &e) { return *this = *this * e; }
    template <class T, EnableIfNotExpression<T> * = nullptr> ThisType &operator+=(const T &v) { return *this = *this + v; }
    template <class T, EnableIfNotExpression<T> * = nullptr> ThisType &operator-=(const T &v) { return *this = *this - v; }
    template <class T, EnableIfNotExpression<T> * = nullptr> ThisType &operator*=(const T &v) { return *this = *this * v; }

    HITNLLS_EXPRESSION_DEFINE_NORM
    HITNLLS_EXPRESSION_DEFINE_TRACE
    HITNLLS_INLINE ThisType &Normalize() {
        ElementType norm = Norm();
        for (int i = 0; i < Rows(); ++i) {
            for (int j = 0; j < Cols(); ++j) {
                At(i, j) /= norm;
            }
        }
        return *this;
    }
private:
    E &e_;
    LayoutType layout_;
};

template <typename E, bool CONST_QUALIFIED = MapTypeTraits<E>::CONST_QUATIFIED>
class Map;
template <typename E>
class Map<E, false> : public Expression<Map<E>> {
public:
    using ThisType = Map<E>;
    using EType = typename MapTypeTraits<E>::EType;
    using ElementType = typename MapTypeTraits<E>::ElementType;
    using ShapeType = typename MapTypeTraits<E>::ShapeType;
    using EvalReturnType = typename MapTypeTraits<E>::EvalReturnType;
    using LayoutType = typename MapTypeTraits<E>::LayoutType;
    using ElementPtrType = typename MapTypeTraits<E>::ElementPtrType;
    static const bool ALIASING = false;

    HITNLLS_INLINE explicit Map(ElementPtrType data, int m, int n) : data_(data), layout_(m, n) {}
    HITNLLS_INLINE explicit Map(ElementPtrType data, int m) : data_(data), layout_(m, m) { HITNLLS_CHECK_VECTOR(ThisType) }
    HITNLLS_INLINE explicit Map(ElementPtrType data) : data_(data) { HITNLLS_CHECK_FIXSIZE_MATRIX(ThisType) }

    HITNLLS_INLINE void Resize(int m, int n) { layout_.Resize(m, n); }
    HITNLLS_INLINE BlockExpression<ThisType> Block(int i, int j, int p, int q, int r = 1, int c = 1) const { return BlockExpression<ThisType>(*this, i, j, p, q, r, c); }
    HITNLLS_INLINE MutableBlockExpression<ThisType> Block(int i, int j, int p, int q, int r = 1, int c = 1) { 
        return MutableBlockExpression<ThisType>(*this, i, j, p, q, r, c);
    }
    HITNLLS_INLINE BlockExpression<ThisType> Row(int i) const { return BlockExpression<ThisType>(*this, i, 0, 1, Cols()); }
    HITNLLS_INLINE BlockExpression<ThisType> Col(int j) const { return BlockExpression<ThisType>(*this, 0, j, Rows(), 1); }
    HITNLLS_INLINE MutableBlockExpression<ThisType> Row(int i) { return MutableBlockExpression<ThisType>(*this, i, 0, 1, Cols()); }
    HITNLLS_INLINE MutableBlockExpression<ThisType> Col(int j) { return MutableBlockExpression<ThisType>(*this, 0, j, Rows(), 1); }
    template <class EO>
    HITNLLS_INLINE void Assign(const Expression<EO> &e, bool aliasing) {
        Resize(e.Rows(), e.Cols());
        if (aliasing) {
            EvalReturnType tmp(e);
            Evaluator::Assign(*this, e.Cast());
        } else {
            Evaluator::Assign(*this, e.Cast());
        }
    }

    HITNLLS_INLINE const ElementType &At(int i, int j) const { return data_[Offset(i, j)];; }
    HITNLLS_INLINE ElementType &At(int i, int j) { return data_[Offset(i, j)]; }
    HITNLLS_INLINE int Rows() const { return layout_.Rows(); }
    HITNLLS_INLINE int Cols() const { return layout_.Cols(); }
    HITNLLS_INLINE int Size() const { return layout_.Size(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return layout_.Shape(); }
    HITNLLS_INLINE int Offset(int i, int j) const { return layout_.Offset(i, j); }

    HITNLLS_INLINE const ElementType &operator()(int i, int j) const { return At(i, j); }
    HITNLLS_INLINE ElementType &operator()(int i, int j) { return At(i, j); }

    template <class EO>
    HITNLLS_INLINE ThisType &operator=(const Expression<EO> &e) {
        Assign(e, e.Aliasing() && !Expression<ThisType>::no_aliasing_);
        return *this;
    }
    template <class T, EnableIfNotExpression<T> * = nullptr>
    HITNLLS_INLINE ThisType &operator=(const T &v) {
        ScalarExpression<T> e(v); 
        Assign(e, e.Aliasing() && !Expression<ThisType>::no_aliasing_);
        return *this;
    }
    template <class EO> ThisType &operator+=(const Expression<EO> &e) { return *this = *this + e; }
    template <class EO> ThisType &operator-=(const Expression<EO> &e) { return *this = *this - e; }
    template <class EO> ThisType &operator*=(const Expression<EO> &e) { return *this = *this * e; }
    template <class T, EnableIfNotExpression<T> * = nullptr> ThisType &operator+=(const T &v) { return *this = *this + v; }
    template <class T, EnableIfNotExpression<T> * = nullptr> ThisType &operator-=(const T &v) { return *this = *this - v; }
    template <class T, EnableIfNotExpression<T> * = nullptr> ThisType &operator*=(const T &v) { return *this = *this * v; }
    HITNLLS_INLINE ThisType &Normalize() {
        ElementType norm = Norm();
        for (int i = 0; i < Rows(); ++i) {
            for (int j = 0; j < Cols(); ++j) {
                At(i, j) /= norm;
            }
        }
        return *this;
    }
    HITNLLS_EXPRESSION_DEFINE_NORM
    HITNLLS_EXPRESSION_DEFINE_TRACE
private:
    ElementPtrType data_;
    LayoutType layout_;
};
template <typename E>
class Map<E, true> : public Expression<Map<E>> {
public:
    using ThisType = Map<E>;
    using EType = typename MapTypeTraits<E>::EType;
    using ElementType = typename MapTypeTraits<E>::ElementType;
    using ShapeType = typename MapTypeTraits<E>::ShapeType;
    using EvalReturnType = typename MapTypeTraits<E>::EvalReturnType;
    using LayoutType = typename MapTypeTraits<E>::LayoutType;
    using ElementPtrType = typename MapTypeTraits<E>::ElementPtrType;
    static const bool ALIASING = false;

    HITNLLS_INLINE explicit Map(ElementPtrType data, int m, int n) : data_(data), layout_(m, n) {}
    HITNLLS_INLINE explicit Map(ElementPtrType data, int m) : data_(data), layout_(m, m) { HITNLLS_CHECK_VECTOR(ThisType) }
    HITNLLS_INLINE explicit Map(ElementPtrType data) : data_(data) { HITNLLS_CHECK_FIXSIZE_MATRIX(ThisType) }

    HITNLLS_INLINE void Resize(int m, int n) { layout_.Resize(m, n); }
    HITNLLS_INLINE BlockExpression<ThisType> Block(int i, int j, int p, int q, int r = 1, int c = 1) const { return BlockExpression<ThisType>(*this, i, j, p, q, r, c); }
    HITNLLS_INLINE BlockExpression<ThisType> Row(int i) const { return BlockExpression<ThisType>(*this, i, 0, 1, Cols()); }
    HITNLLS_INLINE BlockExpression<ThisType> Col(int j) const { return BlockExpression<ThisType>(*this, 0, j, Rows(), 1); }

    HITNLLS_INLINE const ElementType &At(int i, int j) const { return data_[Offset(i, j)];; }
    HITNLLS_INLINE int Rows() const { return layout_.Rows(); }
    HITNLLS_INLINE int Cols() const { return layout_.Cols(); }
    HITNLLS_INLINE int Size() const { return layout_.Size(); }
    HITNLLS_INLINE const ShapeType &Shape() const { return layout_.Shape(); }
    HITNLLS_INLINE int Offset(int i, int j) const { return layout_.Offset(i, j); }

    HITNLLS_EXPRESSION_DEFINE_NORM
    HITNLLS_EXPRESSION_DEFINE_TRACE

    HITNLLS_INLINE const ElementType &operator()(int i, int j) const { return At(i, j); }

    HITNLLS_INLINE ThisType &Normalize() {
        ElementType norm = Norm();
        for (int i = 0; i < Rows(); ++i) {
            for (int j = 0; j < Cols(); ++j) {
                At(i, j) /= norm;
            }
        }
        return *this;
    }
private:
    ElementPtrType data_;
    LayoutType layout_;
};

} // namespace matrix
} // namespace hitnlls