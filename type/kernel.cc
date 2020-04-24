#include "type/kernel.h"

namespace hitnlls {

matrix::Vector3d HuberKernel::Robustify(double chi) {
    matrix::Vector3d result;
    double delta_sq = delta_ * delta_;
    if (chi < delta_sq) {
        result[0] = chi;
        result[1] = 1;
        result[2] = 0;
    } else {
        double error_sqrt = std::sqrt(chi);
        result[0] = 2 * error_sqrt * delta_ - delta_sq;
        result[1] = delta_ / error_sqrt;
        result[2] = - 0.5 * result[1] / chi;
    }
    return result;
}

HITNLLS_REGISTER_KERNEL(HuberKernel)

} // namespace hitnlls