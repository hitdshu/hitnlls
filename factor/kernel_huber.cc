#include "factor/kernel_huber.h"
#include <cmath>

namespace hitnlls {
namespace factor {

::hitnlls::matrix::Vector3f HuberKernel::Robustify(float error) {
    ::hitnlls::matrix::Vector3f result;
    float delta_sq = delta_ * delta_;
    if (error < delta_sq) {
        result[0] = error;
        result[1] = 1;
        result[2] = 0;
    } else {
        float error_sqrt = sqrt(error);
        result[0] = 2 * error_sqrt * delta_ - delta_sq;
        result[1] = delta_ / error_sqrt;
        result[2] = - 0.5 * result[1] / error;
    }
    return result;
}

HITNLLS_REGISTER_KERNEL(HuberKernel);

} // namespace factor
} // namespace hitnlls