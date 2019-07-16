#pragma once

#include "matrix/matrix.h"

namespace hitnlls {
namespace node {

namespace {
float Mod2PI(const float &theta) {
    const float pi = 3.1415926;
    float result = theta;
    while (result < - pi) {
        result += pi;
    }
    while (result > pi) {
        result -= pi;
    }
    return result;
}
}

class SO2 {
public:
    SO2(float theta = 0) { theta_ = theta; theta_ = Mod2PI(theta_); }

    SO2 operator+(const float &val) const {
        SO2 result(*this);
        result.theta_ += val;
        result.theta_ = Mod2PI(result.theta_);
        return result;
    }
    SO2 operator+(const ::hitnlls::matrix::Vector1f &val) const {
        SO2 result(*this);
        result.theta_ += val(0, 0);
        result.theta_ = Mod2PI(result.theta_);
        return result;
    }
    SO2 operator-(const SO2 &rotp) const {
        SO2 result(*this);
        result.theta_ -= rotp.theta_;
        result.theta_ = Mod2PI(result.theta_);
        return result;
    }
    SO2 operator*(const SO2 &roti) const {
        SO2 result(*this);
        result.theta_ += roti.theta_;
        result.theta_ = Mod2PI(result.theta_);
        return result;
    }
    SO2 &operator+=(const float &val) {
        theta_ += val;
        theta_ = Mod2PI(theta_);
        return *this;
    }
    SO2 &operator+=(const ::hitnlls::matrix::Vector1f &val) {
        theta_ += val(0, 0);
        theta_ = Mod2PI(theta_);
        return *this;
    }
    SO2 &operator-=(const SO2 &rotp) {
        theta_ -= rotp.theta_;
        theta_ = Mod2PI(theta_);
        return *this;
    }
    SO2 &operator*=(const SO2 &roti) {
        theta_ += roti.theta_;
        theta_ = Mod2PI(theta_);
        return *this;
    }

    ::hitnlls::matrix::Matrix22f ToMat22f() const {
        ::hitnlls::matrix::Matrix22f result;
        result(0, 0) = cosf(theta_);
        result(0, 1) = - sinf(theta_);
        result(1, 0) = sinf(theta_);
        result(1, 1) = cosf(theta_);
        return result;
    }

    float GetTheta() const {
        return theta_;
    }

    friend ::std::ostream &operator<<(::std::ostream &out, const SO2 &so2) { out << so2.ToMat22f(); return out; }

private:
    float theta_;
};

} // namespace node
} // namespace hitnlls