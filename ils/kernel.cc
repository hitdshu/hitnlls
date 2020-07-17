#include "kernel.h"

namespace nlls {

Vector2f HuberKernel::Robustify(float chi) const {
    float err = sqrt(chi);
    Vector2f result;
    if (err < delta_) {
        result[0] = chi / 2.0;
        result[1] = 1.0;
    } else {
        result[0] = err * delta_ - delta_ * delta_ / 2;
        result[1] = delta_ / err;
    }
    return result;
}

Vector2f CauchyKernel::Robustify(float chi) const {
    Vector2f result;
    result[0] = c_ * c_ / 2 * log(1 + chi / c_ / c_);
    result[1] = 1 / (1 + chi / c_ / c_);
    return result;
}

Vector2f GemanMcClureKernel::Robustify(float chi) const {
    Vector2f result;
    result[0] = chi / 2.0 / (sigma_ + chi);
    result[1] = sigma_ / pow(sigma_ + chi, 2);
    return result;
}

Vector2f ScDcsKernel::Robustify(float chi) const {
    Vector2f result;
    if (chi < phi_) {
        result[0] = chi / 2.0;
        result[1] = 1;
    } else {
        result[0] = 2 * phi_ * chi / (phi_ + chi) - phi_ / 2;
        result[1] = 4 * pow(phi_, 2) / pow(phi_ + chi, 2);
    }
    return result;
}

} // namespace nlls