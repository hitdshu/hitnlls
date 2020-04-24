#pragma once

#include "geometry/macro.h"
#include "geometry/base.h"
#include "geometry/helper.h"
#include "geometry/jet.h"

namespace hitnlls {
namespace geometry {

template <typename T> class SE3;
template <typename T> class SE3Tangent;

template <typename T>
struct LieGroupTraits<SE3<T>> {
    static const int DIM = 3;
    static const int DOF = 6;
    using Scalar = T;
    using Vector = matrix::Matrix<Scalar, DIM, 1>;
    using Storage = matrix::Matrix<Scalar, 7, 1>;
    using Adjoint = matrix::Matrix<Scalar, DOF, DOF>;
    using LieGroup = SE3<T>;
    using Tangent = SE3Tangent<T>;

    using Rotation = matrix::Matrix<Scalar, 3, 3>;
    using Transformation = matrix::Matrix<Scalar, 4, 4>;
};

template <typename T>
struct TangentTraits<SE3Tangent<T>> {
    static const int DIM = 3;
    static const int DOF = 6;
    using Scalar = T;
    using Storage = matrix::Matrix<Scalar, DOF, 1>;
    using Vector = matrix::Matrix<Scalar, 3, 1>;
    using LieGroup = SE3<T>;
    using Tangent = SE3Tangent<T>;
    using Rotation = matrix::Matrix<Scalar, 3, 3>;
};

template <typename T>
class SE3 : public LieGroupBase<SE3<T>> {
public:
    using ThisType = SE3<T>;
    using BaseType = LieGroupBase<SE3<T>>;
    HITNLLS_INHERIT_GROUP_PROPERTIES
    using Rotation = typename LieGroupTraits<ThisType>::Rotation;
    using Transformation = typename LieGroupTraits<ThisType>::Transformation;

    SE3(const Vector &t, const Scalar &w, const Scalar &x, const Scalar &y, const Scalar &z) {
        W() = w; X() = x; Y() = y; Z() = z; Normalize(); 
        TX() = t[0]; TY() = t[1]; TZ() = t[2];
    }
    SE3(const Vector &t, const Vector &ax) {
        Scalar theta = ax.Norm();
        if (theta < Scalar(Constants<Scalar>::EPS)) {
            W() = Scalar(1); X() = ax[0] / Scalar(2); Y() = ax[1] / Scalar(2); Z() = ax[2] / Scalar(2); 
        } else {
            Vector axis = ax.Normalized();
            Scalar s = sin(theta / Scalar(2.0));
            W() = cos(theta / Scalar(2.0)); X() = axis[0] * s; Y() = axis[1] * s; Z() = axis[2] * s;
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
    SE3() { W() = 1.0; X() = 0.0; Y() = 0.0; Z() = 0.0; TX() = 0; TY() = 0; TZ() = 0; }

    LieGroup Inverse() const { return LieGroup(-ToRotation().Transpose()*TB(), W(), -X(), -Y(), -Z()); }
    Tangent Log() const {
        Scalar s = sqrt(X()*X() + Y()*Y() + Z()*Z());
        Vector ax;
        Scalar theta = Scalar(2.0) * atan2(s, W());
        if (abs(theta) < Constants<Scalar>::EPS) {
            ax << X()*Scalar(2), Y()*Scalar(2), Z()*Scalar(2);
        } else {
            ax << X()/s*theta, Y()/s*theta, Z()/s*theta;
        }
        Rotation V;
        Rotation hat = Skew(ax);
        if (abs(theta) < Constants<Scalar>::EPS) {
            V = Rotation::Identity() - Scalar(0.5) * hat;
        } else {
            V = Rotation::Identity() + (Scalar(1) - cos(theta)) / (theta * theta) * hat + (theta - sin(theta)) / (theta * theta * theta) * hat * hat;
        }
        return Tangent(V.Inverse() * TB(), ax);
    }
    Adjoint Adj() const {
        Adjoint adj;
        adj.SetZero();
        Rotation rot = ToRotation();
        adj.Block(0, 0, 3, 3) = rot;
        adj.Block(3, 3, 3, 3) = rot;
        Vector t = ToTranslation();
        Rotation that = Skew(t);
        adj.Block(3, 0, 3, 3) = that * rot;
        return adj;
    }
    LieGroup Compose(const LieGroup &g) const {
        Rotation r1 = ToRotation();
        Vector t1 = ToTranslation();
        Rotation r2 = g.ToRotation();
        Vector t2 = g.ToTranslation();
        return LieGroup(r1*t2+t1, Rotation(r1*r2));
    }
    Vector Act(const Vector &v) const { return ToRotation() * v + ToTranslation(); }

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
    Transformation ToTransform() const {
        Transformation tf = Transformation::Identity();
        tf.Block(0, 0, 3, 3) = ToRotation();
        tf.Block(0, 3, 3, 1) = ToTranslation();
        return tf;
    }
    Vector ToTranslation() const { return TB(); }
    void Normalize() { BaseType::data_.Block(0, 0, 4, 1).Normalize(); }
    Scalar &W() { return BaseType::data_[0]; }
    const Scalar &W() const { return BaseType::data_[0]; }
    Scalar &X() { return BaseType::data_[1]; }
    const Scalar &X() const { return BaseType::data_[1]; }
    Scalar &Y() { return BaseType::data_[2]; }
    const Scalar &Y() const { return BaseType::data_[2]; }
    Scalar &Z() { return BaseType::data_[3]; }
    const Scalar &Z() const { return BaseType::data_[3]; }
    Scalar &TX() { return BaseType::data_[4]; }
    const Scalar &TX() const { return BaseType::data_[4]; }
    Scalar &TY() { return BaseType::data_[5]; }
    const Scalar &TY() const { return BaseType::data_[5]; }
    Scalar &TZ() { return BaseType::data_[6]; }
    const Scalar &TZ() const { return BaseType::data_[6]; }
    Vector TB() const { return BaseType::data_.Block(4, 0, 3, 1); }
};

template <typename T>
class SE3Tangent : public TangentBase<SE3Tangent<T>> {
public:
    using ThisType = SE3Tangent<T>;
    using BaseType = TangentBase<SE3Tangent<T>>;
    HITNLLS_INHERIT_TANGENT_PROPERTIES
    using Vector = typename TangentTraits<ThisType>::Vector;
    using Rotation = typename TangentTraits<ThisType>::Rotation;

    SE3Tangent(const Vector &t, const Scalar &x, const Scalar &y, const Scalar &z) {
        X() = x; Y() = y; Z() = z;
        TX() = t[0]; TY() = t[1]; TZ() = t[2];
    }
    SE3Tangent(const Vector &t, const Vector &ax) {
        X() = ax[0]; Y() = ax[1]; Z() = ax[2];
        TX() = t[0]; TY() = t[1]; TZ() = t[2];
    }
    SE3Tangent() { X() = 0.0; Y() = 0.0; Z() = 0.0; TX() = 0; TY() = 0; TZ() = 0; }

    LieGroup Exp() const {
        Vector ax; 
        ax << X(), Y(), Z(); 
        Vector t;
        t << TX(), TY(), TZ();
        Scalar theta_sq = X()*X() + Y()*Y() + Z()*Z();
        Scalar theta = sqrt(theta_sq);
        Rotation V;
        Rotation hat = Skew(ax);
        if (abs(theta) < Scalar(Constants<Scalar>::EPS)) {
            V = Rotation::Identity() - Scalar(0.5) * hat;
        } else {
            V = Rotation::Identity() + (Scalar(1) - cos(theta)) / (theta * theta) * hat + (theta - sin(theta)) / (theta * theta * theta) * hat * hat;
        }
        return LieGroup(V * t, ax);
    }

    Scalar &X() { return BaseType::data_[0]; }
    const Scalar &X() const { return BaseType::data_[0]; }
    Scalar &Y() { return BaseType::data_[1]; }
    const Scalar &Y() const { return BaseType::data_[1]; }
    Scalar &Z() { return BaseType::data_[2]; }
    const Scalar &Z() const { return BaseType::data_[2]; }
    Scalar &TX() { return BaseType::data_[3]; }
    const Scalar &TX() const { return BaseType::data_[3]; }
    Scalar &TY() { return BaseType::data_[4]; }
    const Scalar &TY() const { return BaseType::data_[4]; }
    Scalar &TZ() { return BaseType::data_[5]; }
    const Scalar &TZ() const { return BaseType::data_[5]; }
};

using SE3d = SE3<double>;
using SE3f = SE3<float>;
using SE3dTangent = SE3Tangent<double>;
using SE3fTangent = SE3Tangent<float>;

} // namespace geometry
} // namespace hitnlls