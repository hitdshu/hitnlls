#pragma once

#include "../matrix/dense.h"

namespace nlls {

struct Equidist {
    static constexpr int DOF = 4;
    static Vector3f Undistort(const Vector3f &pt, const float *coeff = nullptr);
    static Vector3f Distort(const Vector3f &pt, const float *coeff, Matrix33f *dddp = nullptr, Matrix34f *dddc = nullptr);
    Equidist() = delete;
};

} // namespace nlls