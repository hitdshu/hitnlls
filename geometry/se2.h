#pragma once

#include "geometry/macro.h"
#include "geometry/base.h"
#include "geometry/jet.h"

namespace hitnlls {
namespace geometry {

template <typename T> class SE2;
template <typename T> class SE2Tangent;

template <typename T>
struct LieGroupTraits<SE2<T>> {
    static const int DIM = 2;
    static const int DOF = 3;
    using Scalar = T;
    using Vector = matrix::Matrix<Scalar, DIM, 1>;
    using Storage = matrix::Matrix<Scalar, 4, 1>;
    using Adjoint = matrix::Matrix<Scalar, DOF, DOF>;
    using LieGroup = SE2<T>;
    using Tangent = SE2Tangent<T>;

    using Rotation = matrix::Matrix<Scalar, 2, 2>;
    using Transformation = matrix::Matrix<Scalar, 3, 3>;
    using Translation = matrix::Matrix<Scalar, 2, 1>;
};

template <typename T>
struct TangentTraits<SE2Tangent<T>> {
    static const int DIM = 2;
    static const int DOF = 3;
    using Scalar = T;
    using Storage = matrix::Matrix<Scalar, DOF, 1>;
    using LieGroup = SE2<T>;
    using Tangent = SE2Tangent<T>;
};

template <typename T>
class SE2 : public LieGroupBase<SE2<T>> {
public:
    using ThisType = SE2<T>;
    using BaseType = LieGroupBase<SE2<T>>;
    HITNLLS_INHERIT_GROUP_PROPERTIES
    using Rotation = typename LieGroupTraits<ThisType>::Rotation;
    using Transformation = typename LieGroupTraits<ThisType>::Transformation;
    using Translation = typename LieGroupTraits<ThisType>::Translation;

    SE2(const Scalar &x, const Scalar &y, const Scalar &real, const Scalar &imag) { Real() = real; Imag() = imag; X() = x; Y() = y; }
    SE2(const Scalar &x, const Scalar &y, const Scalar &theta) { Real() = cos(theta); Imag() = sin(theta); X() = x; Y() = y; }
    SE2() { Real() = 1.0; Imag() = 0.0; X() = 0; Y() = 0; }

    LieGroup Inverse() const { return LieGroup(-Real()*X()-Imag()*Y(), Imag()*X()-Real()*Y(), Real(), -Imag()); }
    Tangent Log() const {
        Scalar theta = Angle();
        Scalar A;
        Scalar B;
        if (abs(theta) < Constants<Scalar>::EPS) {
            A = 1.0 - theta * theta / 6.0;
            B = theta / 2.0 - theta * theta * theta / 24.0;
        } else {
            A = Imag() / theta;
            B = (1.0 - Real()) / theta;
        }
        Scalar nz = 1.0 / (A * A + B * B);
        A *= nz;
        B *= nz;
        return Tangent(A*X()+B*Y(), -B*X()+A*Y(), theta);
    }
    Adjoint Adj() const {
        Adjoint adj;
        adj.SetIdentity();
        adj.block(0, 0, 2, 2) = ToRotation();
        adj(0, 2) = Y();
        adj(1, 2) = -X();
        return adj;
    }
    LieGroup Compose(const LieGroup &g) const {
        Scalar cr = Real() * g.Real() - Imag() * g.Imag();
        Scalar ci = Real() * g.Imag() + Imag() * g.Real();
        Vector t = ToRotation() * g.Translation() + Translation();
        LieGroup result(t[0], t[1], cr, ci);
        return result;
    }
    Vector Act(const Vector &v) const { return ToRotation() * v + ToTranslation(); }

    Rotation ToRotation() const {
        Rotation rot;
        rot << Real(), -Imag(), Imag(), Real();
        return rot;
    }
    Transformation ToTransform() const {
        Transformation tf = Transformation::Identity();
        tf.Block(0, 0, 2, 2) = ToRotation();
        tf(0, 2) = X();
        tf(1, 2) = Y();
        return tf;
    }
    Translation ToTranslation() const {
        Translation tl;
        tl << X(), Y();
        return tl;
    }
    void Normalize() { BaseType::data_.Block(0 ,0, 2, 1).Normalize(); }
    Scalar Angle() const { return atan2(Imag(), Real()); }
    Scalar &X() { return BaseType::data_[2]; }
    const Scalar &X() const { return BaseType::data_[2]; }
    Scalar &Y() { return BaseType::data_[3]; }
    const Scalar &Y() const { return BaseType::data_[3]; }
    Scalar &Real() { return BaseType::data_[0]; }
    const Scalar &Real() const { return BaseType::data_[0]; }
    Scalar &Imag() { return BaseType::data_[1]; }
    const Scalar &Imag() const { return BaseType::data_[1]; }
};

template <typename T>
class SE2Tangent : public TangentBase<SE2Tangent<T>> {
public:
    using ThisType = SE2Tangent<T>;
    using BaseType = TangentBase<SE2Tangent<T>>;
    HITNLLS_INHERIT_TANGENT_PROPERTIES

    SE2Tangent(const Scalar &x, const Scalar &y, const Scalar &theta) { X() = x; Y() = y; Angle() = theta; }
    SE2Tangent() { Angle() = 0.0; X() = 0; Y() = 0; }

    LieGroup Exp() const {
        Scalar theta = Angle();
        Scalar real = cos(theta);
        Scalar imag = sin(theta);
        Scalar A;
        Scalar B;
        if (abs(theta) < Scalar(Constants<Scalar>::EPS)) {
            A = Scalar(1.0) - theta * theta / Scalar(6.0);
            B = theta / Scalar(2.0) - theta * theta * theta / Scalar(24.0);
        } else {
            A = imag / theta;
            B = (Scalar(1.0) - real) / theta;
        }
        return LieGroup(A*X()-B*Y(), B*X()+A*Y(), theta);
    }

    const Scalar &Angle() const { return BaseType::data_[0]; }
    Scalar &Angle() { return BaseType::data_[0]; }
    const Scalar &X() const { return BaseType::data_[1]; }
    Scalar &X() { return BaseType::data_[1]; }
    const Scalar &Y() const { return BaseType::data_[2]; }
    Scalar &Y() { return BaseType::data_[2]; }
};

using SE2d = SE2<double>;
using SE2f = SE2<float>;
using SE2dTangent = SE2Tangent<double>;
using SE2fTangent = SE2Tangent<float>;

} // namespace geometry
} // namespace hitnlls