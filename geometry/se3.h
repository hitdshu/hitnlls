#pragma once

#include "../matrix/dense.h"

#include "macros.h"
#include "traits.h"
#include "utils.h"
#include "base.h"

namespace nlls {
namespace internal {
template <typename T>
struct LieGroupTraits<SE3<T>> {
    static const int DIM = 3;
    static const int DOF = 6;
    using Scalar = T;
    using LieGroup = SE3<T>;
    using Tangent = SE3Tangent<T>;
    using Storage = Matrix<Scalar, 7, 1>;
};
template <typename T>
struct TangentTraits<SE3Tangent<T>> {
    static const int DIM = 3;
    static const int DOF = 6;
    using Scalar = T;
    using LieGroup = SE3<T>;
    using Tangent = SE3Tangent<T>;
};
} // namespace internal

template <typename T>
class SE3 : public internal::LieGroupBase<SE3<T>> {
public:
    using ThisType = SE3<T>;
    using BaseType = internal::LieGroupBase<SE3<T>>;
    NLLS_INHERIT_GROUP_PROPERTIES
    NLLS_EXTRACT_GROUP_STORAGE(SE3)

    SE3(const Vector &t, const Scalar &w, const Scalar &x, const Scalar &y, const Scalar &z) {
        W() = w; X() = x; Y() = y; Z() = z; Normalize(); 
        TX() = t[0]; TY() = t[1]; TZ() = t[2];
    }
    SE3(const Vector &t, const Vector &ax) {
        Scalar theta = ax.Norm();
        if (theta < internal::Constants<Scalar>::EPS) {
            W() = Scalar(1); X() = ax[0] / Scalar(2); Y() = ax[1] / Scalar(2); Z() = ax[2] / Scalar(2); 
        } else {
            Vector axis = ax.Normalized();
            Scalar s = nlls::sin(theta / Scalar(2.0));
            W() = nlls::cos(theta / Scalar(2.0)); X() = axis[0] * s; Y() = axis[1] * s; Z() = axis[2] * s;
        }
        TX() = t[0]; TY() = t[1]; TZ() = t[2];
    }
    SE3(const Vector &t, const Rotation &rot) {
        W() = sqrt(rot.Trace() + Scalar(1.0)) / Scalar(2.0);
        X() = (rot(2, 1) - rot(1, 2)) / (Scalar(4) * W());
        Y() = (rot(0, 2) - rot(2, 0)) / (Scalar(4) * W());
        Z() = (rot(1, 0) - rot(0, 1)) / (Scalar(4) * W());
        TX() = t[0]; TY() = t[1]; TZ() = t[2];
    }
    SE3(const Transform &tf) {
        Rotation rot = tf.Block(0, 0, 3, 3);
        Vector t = tf.Block(0, 3, 3, 1);
        W() = sqrt(rot.Trace() + Scalar(1.0)) / Scalar(2.0);
        X() = (rot(2, 1) - rot(1, 2)) / (Scalar(4) * W());
        Y() = (rot(0, 2) - rot(2, 0)) / (Scalar(4) * W());
        Z() = (rot(1, 0) - rot(0, 1)) / (Scalar(4) * W());
        TX() = t[0]; TY() = t[1]; TZ() = t[2];
    }
    SE3(const Scalar *d) { for (int i = 0; i < data_.Size(); ++i) { data_[i] = d[i]; } }
    SE3() { W() = Scalar(1.0); X() = Scalar(0.0); Y() = Scalar(0.0); Z() = Scalar(0.0); TX() = Scalar(0); TY() = Scalar(0); TZ() = Scalar(0); }

    LieGroup Inverse() const { return LieGroup(-ToRotation().Transpose()*TB(), W(), -X(), -Y(), -Z()); }
    LieGroup Compose(const LieGroup &g) const {
        Rotation r1 = ToRotation();
        Vector t1 = ToTranslation();
        Rotation r2 = g.ToRotation();
        Vector t2 = g.ToTranslation();
        return LieGroup(r1*t2+t1, Rotation(r1*r2));
    }
    Vector Act(const Vector &v, Dvdl *dl = nullptr) const {
        Rotation rot = ToRotation();
        Vector trans = ToTranslation();
        Vector result = rot * v + trans;
        if (dl) { *dl.Block(0, 0, 3, 3).SetIdentity(); *dl.Block(0, 3, 3, 3) = -Skew3(result); }
        return result;
    }
    Tangent Log() const {
        Scalar s = nlls::sqrt(X()*X() + Y()*Y() + Z()*Z());
        Vector ax;
        Scalar theta = Scalar(2.0) * nlls::atan2(s, W());
        if (nlls::abs(theta) < internal::Constants<Scalar>::EPS) {
            ax << X()*Scalar(2), Y()*Scalar(2), Z()*Scalar(2);
        } else {
            ax << X()/s*theta, Y()/s*theta, Z()/s*theta;
        }
        Rotation V;
        Rotation hat = Skew3(ax);
        if (nlls::abs(theta) < internal::Constants<Scalar>::EPS) {
            V = Rotation::Identity() - Scalar(0.5) * hat;
        } else {
            V = Rotation::Identity() + (Scalar(1) - nlls::cos(theta)) / (theta * theta) * hat + (theta - nlls::sin(theta)) / (theta * theta * theta) * hat * hat;
        }
        Tangent result(V.Inverse() * TB(), ax);
        return result;
    }
    Adjoint Adj() const {
        Adjoint adj;
        adj.SetZero();
        Rotation rot = ToRotation();
        adj.Block(0, 0, 3, 3) = rot;
        adj.Block(3, 3, 3, 3) = rot;
        Vector t = ToTranslation();
        Rotation that = Skew3(t);
        adj.Block(0, 3, 3, 3) = that * rot;
        return adj;
    }
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
        tf.Block(0, 3, 3, 1) = ToTranslation();
        return tf;
    }
    Vector ToTranslation() const { return TB(); }
    void Normalize() { data_.Block(0, 0, 4, 1).Normalize(); }
    Scalar &W() { return data_[0]; }
    const Scalar &W() const { return data_[0]; }
    Scalar &X() { return data_[1]; }
    const Scalar &X() const { return data_[1]; }
    Scalar &Y() { return data_[2]; }
    const Scalar &Y() const { return data_[2]; }
    Scalar &Z() { return data_[3]; }
    const Scalar &Z() const { return data_[3]; }
    Scalar &TX() { return data_[4]; }
    const Scalar &TX() const { return data_[4]; }
    Scalar &TY() { return data_[5]; }
    const Scalar &TY() const { return data_[5]; }
    Scalar &TZ() { return data_[6]; }
    const Scalar &TZ() const { return data_[6]; }
    Vector TB() const { return data_.Block(4, 0, 3, 1); }

private:
    Storage data_;
};

template <typename T>
class SE3Tangent : public internal::TangentBase<SE3Tangent<T>> {
public:
    using ThisType = SE3Tangent<T>;
    using BaseType = internal::TangentBase<SE3Tangent<T>>;
    NLLS_INHERIT_TANGENT_PROPERTIES
    using Rotation = typename LieGroup::Rotation;

    SE3Tangent(const Vector &t, const Scalar &x, const Scalar &y, const Scalar &z) {
        X() = x; Y() = y; Z() = z;
        TX() = t[0]; TY() = t[1]; TZ() = t[2];
    }
    SE3Tangent(const Vector &t, const Vector &ax) {
        X() = ax[0]; Y() = ax[1]; Z() = ax[2];
        TX() = t[0]; TY() = t[1]; TZ() = t[2];
    }
    SE3Tangent(const Scalar *d) {
        TX() = d[0]; TY() = d[1]; TZ() = d[2];
        X() = d[3]; Y() = d[4]; Z() = d[5];
    }
    SE3Tangent() { X() = Scalar(0.0); Y() = Scalar(0.0); Z() = Scalar(0.0); TX() = Scalar(0); TY() = Scalar(0); TZ() = Scalar(0); }
    LieGroup Exp() const {
        Vector ax; 
        ax << X(), Y(), Z(); 
        Vector t;
        t << TX(), TY(), TZ();
        Scalar theta_sq = X()*X() + Y()*Y() + Z()*Z();
        Scalar theta = nlls::sqrt(theta_sq);
        Rotation V;
        Rotation hat = Skew3(ax);
        if (nlls::abs(theta) < internal::Constants<Scalar>::EPS) {
            V = Rotation::Identity() - Scalar(0.5) * hat;
        } else {
            V = Rotation::Identity() + (Scalar(1) - nlls::cos(theta)) / nlls::pow(theta, 2) * hat + (theta - nlls::sin(theta)) / nlls::pow(theta, 3) * hat * hat;
        }
        return LieGroup(V * t, ax);
    }
    LieAlg Generator(int i) const {
        LieAlg result;
        result.SetZero();
        if (i == 0) {
            result(0, 3) = Scalar(1);
        } else if (i == 1) {
            result(1, 3) = Scalar(1);
        } else if (i == 2) {
            result(2, 3) = Scalar(1);
        } else if (i == 3) {
            result(1, 2) = Scalar(-1);
            result(2, 1) = Scalar(1);
        } else if (i == 4) {
            result(0, 2) = Scalar(1);
            result(2, 0) = Scalar(-1);
        } else if (i == 5) {
            result(0, 1) = Scalar(-1);
            result(1, 0) = Scalar(1);
        }
        return result;
    }
    LieAlg Hat() const {
        LieAlg result;
        result.SetZero();
        Vector ax; 
        ax << X(), Y(), Z(); 
        Vector t;
        t << TX(), TY(), TZ();
        result.Block(0, 0, 3, 3) = Skew3(ax);
        result.Block(0, 3, 3, 1) = t;
        return result;
    }

    Vector AngleAxis() const { return Vector(X(), Y(), Z()); }
    Vector TB() const { return Vector(TX(), TY(), TZ()); }
    Scalar Angle() const { return nlls::sqrt(X()*X() + Y()*Y() + Z()*Z()); }
    Scalar &TX() { return BaseType::Data()[0]; }
    const Scalar &TX() const { return BaseType::Data()[0]; }
    Scalar &TY() { return BaseType::Data()[1]; }
    const Scalar &TY() const { return BaseType::Data()[1]; }
    Scalar &TZ() { return BaseType::Data()[2]; }
    const Scalar &TZ() const { return BaseType::Data()[2]; }
    Scalar &X() { return BaseType::Data()[3]; }
    const Scalar &X() const { return BaseType::Data()[3]; }
    Scalar &Y() { return BaseType::Data()[4]; }
    const Scalar &Y() const { return BaseType::Data()[4]; }
    Scalar &Z() { return BaseType::Data()[5]; }
    const Scalar &Z() const { return BaseType::Data()[5]; }
};

using SE3d = SE3<double>;
using SE3f = SE3<float>;
using SE3dTangent = SE3Tangent<double>;
using SE3fTangent = SE3Tangent<float>;

} // namespace nlls