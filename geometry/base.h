#pragma once

#include <iostream>

#include "../common/macros.h"
#include "../matrix/dense.h"

#include "traits.h"

namespace nlls {
namespace internal {
template <class Derived>
struct LieGroupBase {
    static constexpr int DIM = LieGroupTraits<Derived>::DIM;
    static constexpr int DOF = LieGroupTraits<Derived>::DOF;
    using Scalar = typename LieGroupTraits<Derived>::Scalar;
    using LieGroup = typename LieGroupTraits<Derived>::LieGroup;
    using Tangent = typename LieGroupTraits<Derived>::Tangent;
    using Adjoint = Matrix<Scalar, DOF, DOF>;
    using Vector = Matrix<Scalar, DIM, 1>;
    using Rotation = Matrix<Scalar, DIM, DIM>;
    using Transform = Matrix<Scalar, DIM+1, DIM+1>;
    using Dvdl = Matrix<Scalar, DIM, DOF>;
    using Dvdv = Matrix<Scalar, DIM, DIM>;
    using Dldl = Matrix<Scalar, DOF, DOF>;

    NLLS_CRTP_REF

    Tangent Log() const { return Cast().Log(); }
    Adjoint Adj() const { return Cast().Adj(); } 
    Rotation ToRotation() const { return Cast().ToRotation(); }
    Transform ToTransform() const { return Cast().ToTransform(); }

    LieGroup InverseWithDeriv(Dldl *d = nullptr) const {
        LieGroup result = Cast().Inverse();
        if (d) { *d = -result.Adj(); } 
        return result;
    }
    LieGroup ComposeWithDeriv(const LieGroup &g, Dldl *dl = nullptr, Dldl *dr = nullptr) const { 
        if (dl) { *dl = Dldl::Identity(); }
        if (dr) { *dr = Adj(); }
        return Cast().Compose(g);
    }
    Vector ActWithDeriv(const Vector &v, Dvdl *dl = nullptr, Dvdv *dv = nullptr) const {
        if (dv) { *dv = ToRotation(); }
        return Cast().Act(v, dl);
    }
    LieGroup RightPlus(const Tangent &t) const { return Compose(t.Exp()); }
    LieGroup LeftPlus(const Tangent &t) const { return t.Exp().Compose(Cast()); }
    Tangent RightMinus(const LieGroup &g) const { return g.Inverse().Compose(Cast()).Log(); }
    Tangent LeftMinus(const LieGroup &g) const { return Compose(g.Inverse()).Log(); }
    Derived &SetIdentity() { Cast() = Tangent::Zero().Exp(); return Cast(); }

    static const Derived &Identity() { const static Derived e(Derived().SetIdentity()); return e; }
    static Derived Random() { return Tangent::Random().Exp(); }

    NLLS_CRTP_DEC(LieGroupBase)
};

template <class Derived>
struct TangentBase {
    static constexpr int DIM = TangentTraits<Derived>::DIM;
    static constexpr int DOF = TangentTraits<Derived>::DOF;
    using Scalar = typename TangentTraits<Derived>::Scalar;
    using LieGroup = typename TangentTraits<Derived>::LieGroup;
    using Tangent = typename TangentTraits<Derived>::Tangent;
    using LieAlg = Matrix<Scalar, DOF, DOF>;
    using Vector = Matrix<Scalar, DIM, 1>;
    using Storage = Matrix<Scalar, DOF, 1>;

    NLLS_CRTP_REF

    LieGroup Exp() const { return Cast().Exp(); }
    LieAlg Generator(int i) const { return Cast().Generator(i); }
    LieAlg Hat() const { return Cast().Hat(); }

    const Storage &Data() const { return data_; }
    Storage &Data() { return data_; }
    Derived &SetZero() { data_ = Scalar(0.0); return Cast(); }
    Derived &SetRandom() { data_ = Storage::Random(); return Cast(); }

    static const Derived &Zero() { const static Derived e(Derived().SetZero()); return e; }
    static Derived Random() { Derived t; return t.SetRandom(); }

private:
    Storage data_;

    NLLS_CRTP_DEC(TangentBase)
};

template <class Derived>
Derived operator*(const LieGroupBase<Derived> &g1, const typename LieGroupBase<Derived>::LieGroup &g2) {
    return g1.Compose(g2);
}
template <class Derived>
typename LieGroupBase<Derived>::Vector operator*(const LieGroupBase<Derived> &g, const typename LieGroupBase<Derived>::Vector &v) {
    return g.Act(v);
}
template <class Derived>
Derived operator*(const LieGroupBase<Derived> &g, const typename LieGroupBase<Derived>::Tangent &t) {
    return g.RightPlus(t);
}
template <class Derived>
Derived operator*(const typename LieGroupBase<Derived>::Tangent &t, const LieGroupBase<Derived> &g) {
    return g.LeftPlus(t);
}

template <class Derived> std::ostream &operator<<(std::ostream &os, const LieGroupBase<Derived> &l) {
    os << l.ToTransform();
    return os;
}
template <class Derived> std::ostream &operator<<(std::ostream &os, const TangentBase<Derived> &t) {
    os << t.Data().Transpose();
    return os;
}

} // namespace internal
} // namespace nlls