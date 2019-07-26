#pragma once

#include "matrix/matrix.h"
#include "node/so3.h"

namespace hitnlls {
namespace node {

class SE3 {
public:
    SE3(SO3 so3 = SO3(), ::hitnlls::matrix::Vector3f position3 = ::hitnlls::matrix::Vector3f()) {
        so3_ = so3;
        position3_ = position3;
    }

    SE3 operator+(const ::hitnlls::matrix::Matrix<float, 6, 1> &val) const {
        SE3 result(*this);
        result.so3_ += val.Block(0, 0, 3, 1);
        ::hitnlls::matrix::Vector3f position_update = val.Block(3, 0, 3, 1);
        result.position3_ += position_update;
        return result;
    }
    SE3 operator*(const SE3 &posep) const {
        SE3 result(*this);
        result.position3_ += result.so3_.ToMatrix33f() * posep.position3_;
        result.so3_ *= posep.so3_;
        return result;
    }
    SE3 &operator+=(const ::hitnlls::matrix::Matrix<float, 6, 1> &val) {
        so3_ += val.Block(0, 0, 3, 1);
        ::hitnlls::matrix::Vector2f position_update = val.Block(3, 0, 3, 1);
        position3_ += position_update;
        return *this;
    }
    SE3 &operator*=(const SE3 &posep) {
        position3_ += so3_.ToMatrix33f() * posep.position3_;
        so3_ *= posep.so3_;
        return *this;
    }

    float GetPosition(const int &idx) const {
        return position3_[idx];
    }

    ::hitnlls::matrix::Matrix33f ToMat44f() const {
        ::hitnlls::matrix::Matrix44f result;
        ::hitnlls::matrix::Matrix33f so3_result = so3_.ToMatrix33f();
        result(0, 3) = position3_[0];
        result(1, 3) = position3_[1];
        result(2, 3) = position3_[2];
        result(3, 3) = 1;
        result.SetBlock(0, 0, 3, 3, so3_result);
        return result;
    }

    ::hitnlls::matrix::Matrix33f GetRotMat33f() const {
        return so3_.ToMatrix33f();
    }

    ::hitnlls::matrix::Vector3f MapPoint(const ::hitnlls::matrix::Vector3f &point) const {
        ::hitnlls::matrix::Matrix44f pose = ToMat44f();
        ::hitnlls::matrix::Vector4f point_homo;
        point_homo[0] = point[0];
        point_homo[1] = point[1];
        point_homo[2] = point[2];
        point_homo[2] = 1;
        point_homo = pose * point_homo;
        point_homo[0] /= point_homo[3];
        point_homo[1] /= point_homo[3];
        point_homo[2] /= point_homo[3];
        point_homo[3] /= point_homo[3];
        ::hitnlls::matrix::Vector3f result;
        result[0] = point_homo[0];
        result[1] = point_homo[1];
        result[2] = point_homo[2];
        return result;
    }

    friend ::std::ostream &operator<<(::std::ostream &out, const SE3 &se3) { out << se3.ToMat44f(); return out; }

private:
    SO3 so3_;
    ::hitnlls::matrix::Vector3f position3_;
};

} // namespace node
} // namespace hitnlls