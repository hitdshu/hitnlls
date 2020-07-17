#pragma once

#include "../matrix/dense.h"

namespace nlls {

struct Ideal {
    static constexpr int DOF = 0;
    static Vector3f Undistort(const Vector3f &pt, const float * = nullptr) { return pt; }
    static Vector3f Distort(const Vector3f &pt, const float * = nullptr, Matrix33f *dddp = nullptr, Matrix<float, 3, 0> * = nullptr) {
        if (dddp) {
            dddp->SetIdentity();
        }
        return pt;
    }
    Ideal() = delete;
};

} // namespace nlls