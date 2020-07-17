#include "factor.h"

namespace nlls {

void FactorBase::Robustify() {
    if (kernel_ == nullptr) {
        chi_ /= 2.0;
        scale_ = 1.0;
    } else {
        Vector2f result = kernel_->Robustify(chi_);
        chi_ = result[0];
        scale_ = sqrt(result[1]);
    }
}

} // namespace nlls