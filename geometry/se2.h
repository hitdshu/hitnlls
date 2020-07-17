#pragma once

#include "../matrix/dense.h"

#include "macros.h"
#include "traits.h"
#include "utils.h"
#include "base.h"

namespace nlls {
namespace internal {
template <typename T>
struct LieGroupTraits<SE2<T>> {
    static const int DIM = 2;
    static const int DOF = 3;
    using Scalar = T;
    using LieGroup = SE2<T>;
    using Tangent = SE2Tangent<T>;
    using Storage = Matrix<Scalar, 4, 1>;
};
template <typename T>
struct TangentTraits<SE2Tangent<T>> {
    static const int DIM = 2;
    static const int DOF = 3;
    using Scalar = T;
    using LieGroup = SE2<T>;
    using Tangent = SE2Tangent<T>;
};
} // namespace internal

template <typename T>
class SE2 : public internal::LieGroupBase<SE2<T>> {
public:
    using ThisType = SE2<T>;
    using BaseType = internal::LieGroupBase<SE2<T>>;
    NLLS_INHERIT_GROUP_PROPERTIES
    NLLS_EXTRACT_GROUP_STORAGE(SE2)

    SE2(const Scalar &x, const Scalar &y, const Scalar &r, const Scalar &i) { Real() = r; Imag() = i; X() = x; Y() = y; }
    SE2(const Scalar &x, const Scalar &y, const Scalar &t) { Real() = nlls::cos(t); Imag() = nlls::sin(t); X() = x; Y() = y; }
    SE2(const Matrix<Scalar, 3, 1> &xyt) { Real() = nlls::cos(xyt[2]); Imag() = nlls::sin(xyt[2]); X() = xyt[0]; Y() = xyt[1]; }
    SE2(const Transform &tf) { Real() = tf(0, 0); Imag() = tf(1, 0); X() = tf(0, 2); Y() = tf(1, 2); }
    SE2(const Scalar *data) { for (int i = 0; i < data_.Size(); ++i) { data_[i] = data[i]; } }
    SE2() { Real() = Scalar(1.0); Imag() = Scalar(0.0); X() = Scalar(0.0); Y() = Scalar(0.0); }

    LieGroup Inverse() const { return LieGroup(-Real()*X()-Imag()*Y(), Imag()*X()-Real()*Y(), Real(), -Imag()); }
    LieGroup Compose(const LieGroup &g) const {
        Scalar cr = Real() * g.Real() - Imag() * g.Imag();
        Scalar ci = Imag() * g.Real() + Real() * g.Imag();
        Vector t = ToRotation() * g.ToTranslation() + ToTranslation();
        return LieGroup(t[0], t[1], cr, ci);
    }
    Vector Act(const Vector &v, Dvdl *dl = nullptr) const {
        Rotation rot = ToRotation();
        Vector result = rot * v + ToTranslation();
        if (dl) {
            *dl << Scalar(1), Scalar(0), -result[1], 
                Scalar(0), Scalar(1), result[0];
        }
        return result;
    }
    Tangent Log() const {
        Scalar theta = Angle(); 
        Scalar A; 
        Scalar B;
        if (nlls::abs(theta) < internal::Constants<Scalar>::EPS) {
            A = Scalar(1.0) - theta * theta / Scalar(6.0);
            B = theta / Scalar(2.0) - theta * theta * theta / Scalar(24.0);
        } else {
            A = Imag() / theta;
            B = (Scalar(1.0) - Real()) / theta;
        }
        Scalar nz = Scalar(1.0) / (A * A + B * B);
        A *= nz;
        B *= nz;
        return Tangent(A*X()+B*Y(), -B*X()+A*Y(), theta);
    }
    Adjoint Adj() const {
        Adjoint adj;
        adj.SetIdentity();
        adj.Block(0, 0, 2, 2) = ToRotation();
        adj(0, 2) = Y();
        adj(1, 2) = -X();
        return adj;
    }
    Rotation ToRotation() const { Rotation rot; rot << Real(), -Imag(), Imag(), Real(); return rot; }
    Transform ToTransform() const { Transform tf = Transform::Identity(); tf.Block(0, 0, 2, 2) = ToRotation(); tf(0, 2) = X(); tf(1, 2) = Y(); return tf; }
    Vector ToTranslation() const { Vector tl; tl << X(), Y(); return tl; }

    void Normalize() { data_.Block(0 ,0, 2, 1).Normalize(); }
    Scalar Angle() const { return atan2(Imag(), Real()); }
    Scalar &X() { return data_[0]; }
    const Scalar &X() const { return data_[0]; }
    Scalar &Y() { return data_[1]; }
    const Scalar &Y() const { return data_[1]; }
    Scalar &Real() { return data_[2]; }
    const Scalar &Real() const { return data_[2]; }
    Scalar &Imag() { return data_[3]; }
    const Scalar &Imag() const { return data_[3]; }

private:
    Storage data_;
};

template <typename T>
class SE2Tangent : public internal::TangentBase<SE2Tangent<T>> {
public:
    using ThisType = SE2Tangent<T>;
    using BaseType = internal::TangentBase<SE2Tangent<T>>;
    NLLS_INHERIT_TANGENT_PROPERTIES
    SE2Tangent(const Scalar &x, const Scalar &y, const Scalar &theta) { X() = x; Y() = y; Angle() = theta; }
    SE2Tangent(const Matrix<Scalar, 3, 1> &xyt) { X() = xyt[0]; Y() = xyt[1]; Angle() = xyt[2]; }
    SE2Tangent(const Scalar *d) { X() = d[0]; Y() = d[1]; Angle() = d[2]; }
    SE2Tangent() { Angle() = Scalar(0.0); X() = Scalar(0.0); Y() = Scalar(0.0); }

    LieGroup Exp() const {
        Scalar theta = Angle();
        Scalar real = nlls::cos(theta);
        Scalar imag = nlls::sin(theta);
        Scalar A;
        Scalar B;
        if (nlls::abs(theta) < internal::Constants<Scalar>::EPS) {
            A = Scalar(1.0) - theta * theta / Scalar(6.0);
            B = theta / Scalar(2.0) - theta * theta * theta / Scalar(24.0);
        } else {
            A = imag / theta;
            B = (Scalar(1.0) - real) / theta;
        }
        return LieGroup(A*X()-B*Y(), B*X()+A*Y(), theta);
    }
    LieAlg Generator(int i) const {
        LieAlg result;
        result.SetZero();
        if (i == 0) {
            result(0, 2) = Scalar(1);
        } else if (i == 1) {
            result(1, 2) = Scalar(1);
        } else {
            result(0, 1) = Scalar(-1);
            result(1, 0) = Scalar(1);
        }
        return result;
    }
    LieAlg Hat() const {
        LieAlg result;
        result << Scalar(0), -Angle(), X(),
            Angle(), Scalar(0), Y(), 
            Scalar(0), Scalar(0), Scalar(0);
        return result;
    }

    const Scalar &X() const { return BaseType::Data()[0]; }
    Scalar &X() { return BaseType::Data()[0]; }
    const Scalar &Y() const { return BaseType::Data()[1]; }
    Scalar &Y() { return BaseType::Data()[1]; }
    const Scalar &Angle() const { return BaseType::Data()[2]; }
    Scalar &Angle() { return BaseType::Data()[2]; }
};

using SE2d = SE2<double>;
using SE2f = SE2<float>;
using SE2dTangent = SE2Tangent<double>;
using SE2fTangent = SE2Tangent<float>;

} // namespace nlls