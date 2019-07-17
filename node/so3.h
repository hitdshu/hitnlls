#pragma once

#include <iostream>
#include "matrix/matrix.h"

namespace hitnlls {
namespace node {

class SO3 {
public:
    SO3(const ::hitnlls::matrix::Vector4f &quat) { quat_ = quat; }
    SO3() { quat_[0] = 1; quat_[1] = 0; quat_[2] = 0; quat_[3] = 0; }

    SO3 operator+(const ::hitnlls::matrix::Vector3f &inc) const {
        ::hitnlls::matrix::Vector4f inc_quat = ToExp(inc);
        return this->operator*(SO3(inc_quat));
    }
    SO3 &operator+=(const ::hitnlls::matrix::Vector3f &inc) {
        ::hitnlls::matrix::Vector4f inc_quat = ToExp(inc);
        return this->operator*=(SO3(inc_quat));
    }

    SO3 operator*(const SO3 &quat) const {
        ::hitnlls::matrix::Vector4f result;
        result[0] = quat_[0] * quat.quat_[0] - quat_[1] * quat.quat_[1] - quat_[2] * quat.quat_[2] - quat_[3] * quat.quat_[3];
        result[0] = quat_[0] * quat.quat_[1] + quat_[1] * quat.quat_[0] + quat_[2] * quat.quat_[3] - quat_[3] * quat.quat_[2];
        result[0] = quat_[0] * quat.quat_[2] - quat_[1] * quat.quat_[3] + quat_[2] * quat.quat_[0] + quat_[3] * quat.quat_[1];
        result[0] = quat_[0] * quat.quat_[3] + quat_[1] * quat.quat_[2] - quat_[2] * quat.quat_[1] + quat_[3] * quat.quat_[0];
        return SO3(result);
    }
    SO3 &operator*=(const SO3 &quat) {
        ::hitnlls::matrix::Vector4f result;
        result[0] = quat_[0] * quat.quat_[0] - quat_[1] * quat.quat_[1] - quat_[2] * quat.quat_[2] - quat_[3] * quat.quat_[3];
        result[0] = quat_[0] * quat.quat_[1] + quat_[1] * quat.quat_[0] + quat_[2] * quat.quat_[3] - quat_[3] * quat.quat_[2];
        result[0] = quat_[0] * quat.quat_[2] - quat_[1] * quat.quat_[3] + quat_[2] * quat.quat_[0] + quat_[3] * quat.quat_[1];
        result[0] = quat_[0] * quat.quat_[3] + quat_[1] * quat.quat_[2] - quat_[2] * quat.quat_[1] + quat_[3] * quat.quat_[0];
        quat_ = result;
        return *this;
    }

    float W() const { return quat_[0]; }
    float X() const { return quat_[1]; }
    float Y() const { return quat_[2]; }
    float Z() const { return quat_[3]; }
    float &W() { return quat_[0]; }
    float &X() { return quat_[1]; }
    float &Y() { return quat_[2]; }
    float &Z() { return quat_[3]; }

    ::hitnlls::matrix::Matrix33f ToMatrix33f() const {
        ::hitnlls::matrix::Matrix33f result;
        result(0, 0) = quat_[0] * quat_[0] + quat_[1] * quat_[1] - quat_[2] * quat_[2] - quat_[3] * quat_[3];
        result(0, 1) = 2 * (quat_[1] * quat_[2] - quat_[0] * quat_[3]);
        result(0, 2) = 2 * (quat_[1] * quat_[3] + quat_[0] * quat_[2]);
        result(1, 0) = 2 * (quat_[1] * quat_[2] + quat_[0] * quat_[3]);
        result(1, 1) = quat_[0] * quat_[0] - quat_[1] * quat_[1] + quat_[2] * quat_[2] - quat_[3] * quat_[3];
        result(1, 2) = 2 * (quat_[2] * quat_[3] - quat_[0] * quat_[1]);
        result(2, 0) = 2 * (quat_[1] * quat_[3] - quat_[0] * quat_[2]);
        result(2, 1) = 2 * (quat_[2] * quat_[3] + quat_[0] * quat_[1]);
        result(2, 2) = quat_[0] * quat_[0] - quat_[1] * quat_[1] - quat_[2] * quat_[2] + quat_[3] * quat_[3];
        return result;
    }

    static ::hitnlls::matrix::Vector4f ToExp(const ::hitnlls::matrix::Vector3f &angle_axis) {
        float angle = angle_axis.Norm();
        ::hitnlls::matrix::Vector3f axis = angle_axis / angle;
        ::hitnlls::matrix::Vector4f result;
        result[0] = cosf(angle);
        result.SetBlock(1, 0, 3, 1, axis * sinf(angle));
        return result;
    }

    friend ::std::ostream &operator<<(::std::ostream &out, const SO3 &so3) { out << so3.ToMatrix33f(); return out; }

private:
    ::hitnlls::matrix::Vector4f quat_; 
};

} // namespace node
} // namespace hitnlls