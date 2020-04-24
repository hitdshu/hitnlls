#pragma once

#include "geometry/macro.h"
#include "geometry/base.h"
#include "geometry/jet.h"

namespace hitnlls {
namespace geometry {

template <typename T> class SO2;
template <typename T> class SO2Tangent;

template <typename T>
struct LieGroupTraits<SO2<T>> {
    static const int DIM = 2;
    static const int DOF = 1;
    using Scalar = T;
    using Vector = matrix::Matrix<Scalar, DIM, 1>;
    using Storage = matrix::Matrix<Scalar, DIM, 1>;
    using Adjoint = T;
    using LieGroup = SO2<T>;
    using Tangent = SO2Tangent<T>;

    using Rotation = matrix::Matrix<Scalar, DIM, DIM>;
    using Transformation = matrix::Matrix<Scalar, DIM + 1, DIM + 1>;
};

template <typename T>
struct TangentTraits<SO2Tangent<T>> {
    static const int DIM = 2;
    static const int DOF = 1;
    using Scalar = T;
    using Storage = T;
    using LieGroup = SO2<T>;
    using Tangent = SO2Tangent<T>;
};

template <typename T>
class SO2 : public LieGroupBase<SO2<T>> {
public:
    using ThisType = SO2<T>;
    using BaseType = LieGroupBase<SO2<T>>;
    HITNLLS_INHERIT_GROUP_PROPERTIES
    using Rotation = typename LieGroupTraits<ThisType>::Rotation;
    using Transformation = typename LieGroupTraits<ThisType>::Transformation;

    SO2(const Scalar &real, const Scalar &imag) { Real() = real; Imag() = imag; }
    SO2(const Scalar &theta) { Real() = cos(theta); Imag() = sin(theta); }
    SO2(const Rotation &rot) { Real() = rot(0, 0); Imag() = rot(1, 0); }
    SO2() { Real() = 1.0; Imag() = 0.0; }

    LieGroup Inverse() const { return LieGroup(Real(), -Imag()); }
    Tangent Log() const { return Tangent(Angle()); }
    Adjoint Adj() const { return Adjoint(1.0); } 
    LieGroup Compose(const LieGroup &g) const {
        Scalar cr = Real() * g.Real() - Imag() * g.Imag();
        Scalar ci = Real() * g.Imag() + Imag() * g.Real();
        LieGroup result(cr, ci);
        return result;
    }
    Vector Act(const Vector &v) const { return ToRotation() * v; }

    Rotation ToRotation() const {
        Rotation rot;
        rot << Real(), -Imag(), Imag(), Real();
        return rot;
    }
    Transformation ToTransform() const {
        Transformation tf = Transformation::Identity();
        tf.Block(0, 0, 2, 2) = ToRotation();
        return tf;
    }
    void Normalize() { BaseType::data_.Normalize(); }
    Scalar Angle() const { return atan2(Imag(), Real()); }
    Scalar &Real() { return BaseType::data_[0]; }
    const Scalar &Real() const { return BaseType::data_[0]; }
    Scalar &Imag() { return BaseType::data_[1]; }
    const Scalar &Imag() const { return BaseType::data_[1]; }
};

template <typename T>
class SO2Tangent : public TangentBase<SO2Tangent<T>> {
public:
    using ThisType = SO2Tangent<T>;
    using BaseType = TangentBase<SO2Tangent<T>>;
    HITNLLS_INHERIT_TANGENT_PROPERTIES

    SO2Tangent(const Scalar &theta) { Angle() = theta; }
    SO2Tangent() { Angle() = 0.0; }

    LieGroup Exp() const { return LieGroup(cos(Angle()), sin(Angle())); }

    const Scalar &Angle() const { return BaseType::data_; }
    Scalar &Angle() { return BaseType::data_; }
};

using SO2d = SO2<double>;
using SO2f = SO2<float>;
using SO2dTangent = SO2Tangent<double>;
using SO2fTangent = SO2Tangent<float>;

} // namespace geometry
} // namespace hitnlls