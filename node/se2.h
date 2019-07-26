#pragma once

#include "matrix/matrix.h"
#include "node/so2.h"

namespace hitnlls {
namespace node {

class SE2 {
public:
    SE2(SO2 so2 = SO2(), ::hitnlls::matrix::Vector2f position2 = ::hitnlls::matrix::Vector2f()) {
        so2_ = so2;
        position2_ = position2;
    }

    SE2 operator+(const ::hitnlls::matrix::Vector3f &val) const {
        SE2 result(*this);
        result.so2_ += val(0);
        ::hitnlls::matrix::Vector2f position_update;
        position_update[0] = val[1];
        position_update[1] = val[2];
        result.position2_ += position_update;
        return result;
    }
    SE2 operator-(const SE2 &posep) const {
        SE2 result(*this);
        result.so2_ -= posep.so2_;
        result.position2_ -= posep.position2_;
        return result;
    }
    SE2 operator*(const SE2 &posep) const {
        SE2 result(*this);
        result.position2_ += result.so2_.ToMat22f() * posep.position2_;
        result.so2_ *= posep.so2_;
        return result;
    }
    SE2 &operator+=(const ::hitnlls::matrix::Vector3f &val) {
        so2_ += val(0);
        ::hitnlls::matrix::Vector2f position_update;
        position_update[0] = val[1];
        position_update[1] = val[2];
        position2_ += position_update;
        return *this;
    }
    SE2 &operator-=(const SE2 &posep) {
        so2_ -= posep.so2_;
        position2_ -= posep.position2_;
        return *this;
    }
    SE2 &operator*=(const SE2 &posep) {
        position2_ += so2_.ToMat22f() * posep.position2_;
        so2_ *= posep.so2_;
        return *this;
    }

    SE2 Inverse() const {
        SE2 result;
        result.so2_ = SO2(- so2_.GetTheta());
        result.position2_ = - result.so2_.ToMat22f() * position2_;
        return result;
    }

    float GetPosition(const int &idx) const {
        return position2_[idx];
    }

    ::hitnlls::matrix::Matrix33f ToMat33f() const {
        ::hitnlls::matrix::Matrix33f result;
        ::hitnlls::matrix::Matrix22f so2_result = so2_.ToMat22f();
        result(0, 0) = so2_result(0, 0);
        result(0, 1) = so2_result(0, 1);
        result(0, 2) = position2_[0];
        result(1, 0) = so2_result(1, 0);
        result(1, 1) = so2_result(1, 1);
        result(1, 2) = position2_[1];
        result(2, 2) = 1;
        return result;
    }

    ::hitnlls::matrix::Matrix<float, 3, 1> ToVector3f() const {
        ::hitnlls::matrix::Matrix<float, 3, 1> result;
        result[0] = so2_.GetTheta();
        result[1] = position2_[0];
        result[2] = position2_[1];
        return result;
    }

    ::hitnlls::matrix::Matrix22f GetRotMat22f() const {
        return so2_.ToMat22f();
    }

    ::hitnlls::matrix::Vector2f MapPoint(const ::hitnlls::matrix::Vector2f &point) const {
        ::hitnlls::matrix::Matrix33f pose = ToMat33f();
        ::hitnlls::matrix::Vector3f point_homo;
        point_homo[0] = point[0];
        point_homo[1] = point[1];
        point_homo[2] = 1;
        point_homo = pose * point_homo;
        point_homo[0] /= point_homo[2];
        point_homo[1] /= point_homo[2];
        point_homo[2] /= point_homo[2];
        ::hitnlls::matrix::Vector2f result;
        result[0] = point_homo[0];
        result[1] = point_homo[1];
        return result;
    }

    friend ::std::ostream &operator<<(::std::ostream &out, const SE2 &se2) { out << se2.ToMat33f(); return out; }

private:
    SO2 so2_;
    ::hitnlls::matrix::Vector2f position2_;
};

} // namespace node
} // namespace hitnlls