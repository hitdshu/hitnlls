#pragma once

#include "geometry/macro.h"
#include "geometry/base.h"
#include "geometry/jet.h"

namespace hitnlls {
namespace geometry {

template <typename T> class SO3;
template <typename T> class SO3Tangent;

template <typename T>
struct LieGroupTraits<SO3<T>> {
    static const int DIM = 3;
    static const int DOF = 3;
    using Scalar = T;
    using Vector = matrix::Matrix<Scalar, DIM, 1>;
    using Storage = matrix::Matrix<Scalar, 4, 1>;
    using Adjoint = matrix::Matrix<Scalar, DOF, DOF>;
    using LieGroup = SO3<T>;
    using Tangent = SO3Tangent<T>;

    using Rotation = matrix::Matrix<Scalar, 3, 3>;
    using Transformation = matrix::Matrix<Scalar, 4, 4>;
};

template <typename T>
struct TangentTraits<SO3Tangent<T>> {
    static const int DIM = 3;
    static const int DOF = 3;
    using Scalar = T;
    using Storage = matrix::Matrix<Scalar, DOF, 1>;
    using LieGroup = SO3<T>;
    using Tangent = SO3Tangent<T>;
};

template <typename T>
class SO3 : public LieGroupBase<SO3<T>> {
public:
    using ThisType = SO3<T>;
    using BaseType = LieGroupBase<SO3<T>>;
    HITNLLS_INHERIT_GROUP_PROPERTIES
    using Rotation = typename LieGroupTraits<ThisType>::Rotation;
    using Transformation = typename LieGroupTraits<ThisType>::Transformation;

    SO3(const Scalar &w, const Scalar &x, const Scalar &y, const Scalar &z) { W() = w; X() = x; Y() = y; Z() = z; Normalize(); }
    SO3(const Vector &ax) {
        Scalar theta = ax.Norm();
        if (theta < Constants<Scalar>::EPS) {
            W() = 1; X() = ax[0] / 2; Y() = ax[1] / 2; Z() = ax[2] / 2; 
        } else {
            Vector axis = ax.Normalized();
            Scalar s = sin(theta / 2.0);
            W() = cos(theta / 2.0); X() = axis[0] * s; Y() = axis[1] * s; Z() = axis[2] * s;
        }
    }
    SO3(const Rotation &rot) {
        W() = sqrt(rot.Trace() + 1.0) / 2.0;
        X() = (rot(2, 1) - rot(1, 2)) / (4 * W());
        Y() = (rot(0, 2) - rot(2, 0)) / (4 * W());
        Z() = (rot(1, 0) - rot(0, 1)) / (4 * W());
    }
    SO3() { W() = 1.0; X() = 0.0; Y() = 0.0; Z() = 0.0; }

    LieGroup Inverse() const { return LieGroup(W(), -X(), -Y(), -Z()); }
    Tangent Log() const {
        Scalar s = sqrt(X()*X() + Y()*Y() + Z()*Z());
        if (s < Constants<Scalar>::EPS) {
            return Tangent(X()*2, Y()*2, Z()*2);
        } else {
            Scalar theta = 2.0 * atan2(s, W());
            return Tangent(X()/s*theta, Y()/s*theta, Z()/s*theta);
        }
    }
    Adjoint Adj() const { return ToRotation(); }
    LieGroup Compose(const LieGroup &g) const {
        Scalar w1 = W(); Scalar x1 = X(); Scalar y1 = Y(); Scalar z1 = Z();
        Scalar w2 = g.W(); Scalar x2 = g.X(); Scalar y2 = g.Y(); Scalar z2 = g.Z();
        Scalar w = w1*w2 - x1*x2 - y1*y2 - z1*z2;
        Scalar x = w1*x2 + x1*w2 + y1*z2 - z1*y2;
        Scalar y = w1*y2 - x1*z2 + y1*w2 + z1*w2;
        Scalar z = w1*z2 + x1*y2 - y1*x2 + z1*w2;
        return LieGroup(w, x, y, z);
    }
    Vector Act(const Vector &v) const { return ToRotation() * v; }

    Rotation ToRotation() const {
        Rotation rot;
        Scalar w = W(); Scalar x = X(); Scalar y = Y(); Scalar z = Z();
        rot(0, 0) = w * w + x * x - y * y - z * z;
        rot(0, 1) = 2 * (x * y - w * z);
        rot(0, 2) = 2 * (x * z + w * y);
        rot(1, 0) = 2 * (x * y + w * z);
        rot(1, 1) = w * w - x * x + y * y - z * z;
        rot(1, 2) = 2 * (y * z - w * x);
        rot(2, 0) = 2 * (x * z - w * y);
        rot(2, 1) = 2 * (y * z + w * x);
        rot(2, 2) = w * w - x * x - y * y + z * z;
        return rot;
    }
    Transformation ToTransform() const {
        Transformation tf = Transformation::Identity();
        tf.Block(0, 0, 3, 3) = ToRotation();
        return tf;
    }
    void Normalize() { BaseType::data_.Normalize(); }
    Scalar &W() { return BaseType::data_[0]; }
    const Scalar &W() const { return BaseType::data_[0]; }
    Scalar &X() { return BaseType::data_[1]; }
    const Scalar &X() const { return BaseType::data_[1]; }
    Scalar &Y() { return BaseType::data_[2]; }
    const Scalar &Y() const { return BaseType::data_[2]; }
    Scalar &Z() { return BaseType::data_[3]; }
    const Scalar &Z() const { return BaseType::data_[3]; }
};

template <typename T>
class SO3Tangent : public TangentBase<SO3Tangent<T>> {
public:
    using ThisType = SO3Tangent<T>;
    using BaseType = TangentBase<SO3Tangent<T>>;
    HITNLLS_INHERIT_TANGENT_PROPERTIES
    using Vector = Storage;

    SO3Tangent(const Scalar &x, const Scalar &y, const Scalar &z) { X() = x; Y() = y; Z() = z; }
    SO3Tangent(const Vector &ax) { X() = ax[0]; Y() = ax[1]; Z() = ax[2]; }
    SO3Tangent() { X() = 0.0; Y() = 0.0; Z() = 0.0; }

    LieGroup Exp() const { Vector ax; ax << X(), Y(), Z(); return LieGroup(ax); }

    Scalar &X() { return BaseType::data_[0]; }
    const Scalar &X() const { return BaseType::data_[0]; }
    Scalar &Y() { return BaseType::data_[1]; }
    const Scalar &Y() const { return BaseType::data_[1]; }
    Scalar &Z() { return BaseType::data_[2]; }
    const Scalar &Z() const { return BaseType::data_[2]; }
};

using SO3d = SO3<double>;
using SO3f = SO3<float>;
using SO3dTangent = SO3Tangent<double>;
using SO3fTangent = SO3Tangent<float>;

} // namespace geometry
} // namespace hitnlls