#pragma once

#include "../matrix/dense.h"

#include "macros.h"
#include "traits.h"
#include "utils.h"
#include "base.h"

namespace nlls {
namespace internal {
template <typename T>
struct LieGroupTraits<SO2<T>> {
    static constexpr int DIM = 2;
    static constexpr int DOF = 1;
    using Scalar = T;
    using LieGroup = SO2<T>;
    using Tangent = SO2Tangent<T>;
    using Storage = Matrix<Scalar, 2, 1>;
};
template <typename T>
struct TangentTraits<SO2Tangent<T>> {
    static const int DIM = 2;
    static const int DOF = 1;
    using Scalar = T;
    using LieGroup = SO2<T>;
    using Tangent = SO2Tangent<T>;
};
} // namespace internal

template <typename T>
class SO2 : public internal::LieGroupBase<SO2<T>> {
public:
    using ThisType = SO2<T>;
    using BaseType = internal::LieGroupBase<SO2<T>>;
    NLLS_INHERIT_GROUP_PROPERTIES
    NLLS_EXTRACT_GROUP_STORAGE(SO2)

    explicit SO2(const Scalar &real, const Scalar &imag) { Real() = real; Imag() = imag; }
    explicit SO2(const Scalar &theta) { Real() = nlls::cos(theta); Imag() = nlls::sin(theta); }
    explicit SO2(const Rotation &rot) { Real() = rot(0, 0); Imag() = rot(1, 0); }
    explicit SO2(const Scalar *data) { Real() = data[0]; Imag() = data[1]; }
    explicit SO2() { Real() = Scalar(1.0); Imag() = Scalar(0.0); }

    LieGroup Inverse() const { return LieGroup(Real(), -Imag()); }
    LieGroup Compose(const LieGroup &g) const {
        Scalar cr = Real() * g.Real() - Imag() * g.Imag();
        Scalar ci = Real() * g.Imag() + Imag() * g.Real();
        LieGroup result(cr, ci);
        return result;
    }
    Vector Act(const Vector &v, Dvdl *dl = nullptr) const {
        Rotation rot = ToRotation();
        Vector result = rot * v; 
        if (dl) { dl[0] = -result[1]; dl[1] = result[0]; }
        return result;
    }
    Tangent Log() const { Tangent result(Angle()); return result; }
    Adjoint Adj() const { return Adjoint(1.0); } 
    Rotation ToRotation() const { Rotation rot; rot << Real(), -Imag(), Imag(), Real(); return rot; }
    Transform ToTransform() const { Transform tf = Transform::Identity(); tf.Block(0, 0, 2, 2) = ToRotation(); return tf; }

    void Normalize() { data_.Normalize(); }
    Scalar Angle() const { return nlls::atan2(Imag(), Real()); }
    Scalar &Real() { return data_[0]; }
    const Scalar &Real() const { return data_[0]; }
    Scalar &Imag() { return data_[1]; }
    const Scalar &Imag() const { return data_[1]; }

private:
    Storage data_;
};

template <typename T>
class SO2Tangent : public internal::TangentBase<SO2Tangent<T>> {
public:
    using ThisType = SO2Tangent<T>;
    using BaseType = internal::TangentBase<SO2Tangent<T>>;
    NLLS_INHERIT_TANGENT_PROPERTIES

    SO2Tangent(const Scalar &theta) { Angle() = theta; }
    SO2Tangent(const Scalar *data) { Angle() = *data; }
    SO2Tangent() { Angle() = Scalar(0.0); }

    LieGroup Exp() const { return LieGroup(nlls::cos(Angle()), nlls::sin(Angle())); }
    LieAlg Generator(int i) const { return Skew2<Scalar>(Scalar(1.0)); }
    LieAlg Hat() const { return Skew2<Scalar>(Angle()); }

    const Scalar &Angle() const { return BaseType::Data()[0]; }
    Scalar &Angle() { return BaseType::Data()[0]; }
};

using SO2d = SO2<double>;
using SO2f = SO2<float>;
using SO2dTangent = SO2Tangent<double>;
using SO2fTangent = SO2Tangent<float>;

} // namespace nlls