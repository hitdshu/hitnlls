#pragma once

#include "../matrix/dense.h"

#include "macros.h"
#include "traits.h"
#include "utils.h"
#include "base.h"

namespace nlls {
namespace internal {
template <typename T>
struct LieGroupTraits<SO3<T>> {
    static const int DIM = 3;
    static const int DOF = 3;
    using Scalar = T;
    using LieGroup = SO3<T>;
    using Tangent = SO3Tangent<T>;
    using Storage = Matrix<Scalar, 4, 1>;
};
template <typename T>
struct TangentTraits<SO3Tangent<T>> {
    static const int DIM = 3;
    static const int DOF = 3;
    using Scalar = T;
    using LieGroup = SO3<T>;
    using Tangent = SO3Tangent<T>;
};
} // namespace internal

template <typename T>
class SO3 : public internal::LieGroupBase<SO3<T>> {
public:
    using ThisType = SO3<T>;
    using BaseType = internal::LieGroupBase<SO3<T>>;
    NLLS_INHERIT_GROUP_PROPERTIES
    NLLS_EXTRACT_GROUP_STORAGE(SO3)

    SO3(const Scalar &w, const Scalar &x, const Scalar &y, const Scalar &z) { W() = w; X() = x; Y() = y; Z() = z; Normalize(); }
    SO3(const Vector &ax) {
        Scalar theta = ax.Norm();
        if (theta < internal::Constants<Scalar>::EPS) {
            W() = Scalar(1); X() = ax[0] / Scalar(2); Y() = ax[1] / Scalar(2); Z() = ax[2] / Scalar(2); 
        } else {
            Vector axis = ax.Normalized();
            Scalar s = nlls::sin(theta / Scalar(2.0));
            W() = nlls::cos(theta / Scalar(2.0)); X() = axis[0] * s; Y() = axis[1] * s; Z() = axis[2] * s;
        }
    }
    SO3(const Rotation &rot) {
        W() = nlls::sqrt(rot.Trace() + Scalar(1.0)) / Scalar(2.0);
        X() = (rot(2, 1) - rot(1, 2)) / (Scalar(4) * W());
        Y() = (rot(0, 2) - rot(2, 0)) / (Scalar(4) * W());
        Z() = (rot(1, 0) - rot(0, 1)) / (Scalar(4) * W());
    }
    SO3(const Scalar *data) { W() = data[0]; X() = data[1]; Y() = data[2]; Z() = data[3]; }
    SO3() { W() = Scalar(1.0); X() = Scalar(0.0); Y() = Scalar(0.0); Z() = Scalar(0.0); }

    LieGroup Inverse() const { return LieGroup(W(), -X(), -Y(), -Z()); }
    LieGroup Compose(const LieGroup &g) const {
        Scalar w1 = W(); Scalar x1 = X(); Scalar y1 = Y(); Scalar z1 = Z();
        Scalar w2 = g.W(); Scalar x2 = g.X(); Scalar y2 = g.Y(); Scalar z2 = g.Z();
        Scalar w = w1*w2 - x1*x2 - y1*y2 - z1*z2;
        Scalar x = w1*x2 + x1*w2 + y1*z2 - z1*y2;
        Scalar y = w1*y2 - x1*z2 + y1*w2 + z1*x2;
        Scalar z = w1*z2 + x1*y2 - y1*x2 + z1*w2;
        return LieGroup(w, x, y, z);
    }
    Vector Act(const Vector &v, Dvdl *dl = nullptr) const {
        Rotation rot = ToRotation();
        Vector result = rot * v;
        if (dl) { *dl = -Skew3(result); }
        return result;
    }
    Tangent Log() const {
        Scalar s = nlls::sqrt(X()*X() + Y()*Y() + Z()*Z());
        Vector ax;
        if (s < internal::Constants<Scalar>::EPS) {
            ax << X()*Scalar(2), Y()*Scalar(2), Z()*Scalar(2);
        } else {
            Scalar theta = Scalar(2.0) * nlls::atan2(s, W());
            ax << X()/s*theta, Y()/s*theta, Z()/s*theta;
        }
        Tangent result(ax);
        return result;
    }
    Adjoint Adj() const { return ToRotation(); }
    Rotation ToRotation() const {
        Rotation rot;
        Scalar w = W(); Scalar x = X(); Scalar y = Y(); Scalar z = Z();
        rot(0, 0) = w * w + x * x - y * y - z * z;
        rot(0, 1) = Scalar(2) * (x * y - w * z);
        rot(0, 2) = Scalar(2) * (x * z + w * y);
        rot(1, 0) = Scalar(2) * (x * y + w * z);
        rot(1, 1) = w * w - x * x + y * y - z * z;
        rot(1, 2) = Scalar(2) * (y * z - w * x);
        rot(2, 0) = Scalar(2) * (x * z - w * y);
        rot(2, 1) = Scalar(2) * (y * z + w * x);
        rot(2, 2) = w * w - x * x - y * y + z * z;
        return rot;
    }
    Transform ToTransform() const {
        Transform tf = Transform::Identity();
        tf.Block(0, 0, 3, 3) = ToRotation();
        return tf;
    }
    void Normalize() { data_.Normalize(); }
    Scalar &W() { return data_[0]; }
    const Scalar &W() const { return data_[0]; }
    Scalar &X() { return data_[1]; }
    const Scalar &X() const { return data_[1]; }
    Scalar &Y() { return data_[2]; }
    const Scalar &Y() const { return data_[2]; }
    Scalar &Z() { return data_[3]; }
    const Scalar &Z() const { return data_[3]; }

private:
    Storage data_;
};

template <typename T>
class SO3Tangent : public internal::TangentBase<SO3Tangent<T>> {
public:
    using ThisType = SO3Tangent<T>;
    using BaseType = internal::TangentBase<SO3Tangent<T>>;
    NLLS_INHERIT_TANGENT_PROPERTIES

    SO3Tangent(const Scalar &x, const Scalar &y, const Scalar &z) { X() = x; Y() = y; Z() = z; }
    SO3Tangent(const Vector &ax) { X() = ax[0]; Y() = ax[1]; Z() = ax[2]; }
    SO3Tangent(const Scalar *d) { X() = d[0]; Y() = d[1]; Z() = d[2]; }
    SO3Tangent() { X() = Scalar(0.0); Y() = Scalar(0.0); Z() = Scalar(0.0); }

    LieGroup Exp() const {
        Vector ax = AngleAxis();
        return LieGroup(ax);
    }
    LieAlg Generator(int i) const {
        LieAlg result;
        result.SetZero();
        if (i == 0) {
            result(1, 2) = Scalar(-1);
            result(2, 1) = Scalar(1);
        } else if (i == 1) {
            result(0, 2) = Scalar(1);
            result(2, 0) = Scalar(-1);
        } else {
            result(0, 1) = Scalar(-1);
            result(1, 0) = Scalar(1);
        }
        return result;
    }
    LieAlg Hat() const { return Skew3(BaseType::Data()); }

    Vector AngleAxis() const { return Vector(X(), Y(), Z()); }
    Scalar Angle() const { return nlls::sqrt(X()*X() + Y()*Y() + Z()*Z()); }
    Scalar &X() { return BaseType::Data()[0]; }
    const Scalar &X() const { return BaseType::Data()[0]; }
    Scalar &Y() { return BaseType::Data()[1]; }
    const Scalar &Y() const { return BaseType::Data()[1]; }
    Scalar &Z() { return BaseType::Data()[2]; }
    const Scalar &Z() const { return BaseType::Data()[2]; }
};

using SO3d = SO3<double>;
using SO3f = SO3<float>;
using SO3dTangent = SO3Tangent<double>;
using SO3fTangent = SO3Tangent<float>;

} // namespace nlls