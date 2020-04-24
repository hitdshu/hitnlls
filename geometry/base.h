#pragma once

#include "geometry/trait.h"

namespace hitnlls {
namespace geometry {

template <class Derived>
struct LieGroupBase {
    static const int DIM = LieGroupTraits<Derived>::DIM;
    static const int DOF = LieGroupTraits<Derived>::DOF;
    using Scalar = typename LieGroupTraits<Derived>::Scalar;
    using Vector = typename LieGroupTraits<Derived>::Vector;
    using Storage = typename LieGroupTraits<Derived>::Storage;
    using Adjoint = typename LieGroupTraits<Derived>::Adjoint;
    using LieGroup = typename LieGroupTraits<Derived>::LieGroup;
    using Tangent = typename LieGroupTraits<Derived>::Tangent;

    const Derived &Cast() const { return static_cast<const Derived &>(*this); }
    Derived &Cast() { return static_cast<Derived &>(*this); }
    LieGroup Inverse() const { return Cast().Inverse(); }
    Tangent Log() const { return Cast().Log(); }
    Adjoint Adj() const { return Cast().Adj(); } 
    LieGroup Compose(const LieGroup &g) const { return Cast().Compose(g); }
    Vector Act(const Vector &v) const { return Cast().Act(v); }

    LieGroup RightPlus(const Tangent &t) const { return Compose(t.Exp()); }
    LieGroup LeftPlus(const Tangent &t) const { return t.Exp().Compose(Cast()); }
    Tangent RightMinus(const LieGroup &g) const { return g.Inverse().Compose(Cast()).Log(); }
    Tangent LeftMinus(const LieGroup &g) const { return Compose(g.Inverse()).Log(); }
    Derived &SetIdentity() { Cast() = Tangent::Zero().Exp(); return Cast(); }

    static Derived Identity() { const static LieGroup e(LieGroup().SetIdentity()); return e; }

protected:
    Storage data_;
};

template <class Derived>
struct TangentBase {
    static const int DIM = TangentTraits<Derived>::DIM;
    static const int DOF = TangentTraits<Derived>::DOF;
    using Scalar = typename TangentTraits<Derived>::Scalar;
    using Storage = typename TangentTraits<Derived>::Storage;
    using LieGroup = typename TangentTraits<Derived>::LieGroup;
    using Tangent = typename TangentTraits<Derived>::Tangent;

    const Derived &Cast() const { return static_cast<const Derived &>(*this); }
    Derived &Cast() { return static_cast<Derived &>(*this); }
    LieGroup Exp() const { return Cast().Exp(); }

    Derived &SetZero() { data_ = 0; return Cast(); }
    const Storage &ToVector() const { return data_; }

    static Derived Zero() { const static Tangent t(Tangent().SetZero()); return t; }

protected:
    Storage data_;
};

template <class Derived>
Derived operator*(const LieGroupBase<Derived> &g1, const LieGroupBase<Derived> &g2) {
    return g1.RightPlus(g2);
}

} // namespace geometry
} // namespace hitnlls