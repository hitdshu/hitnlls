#include "type/factor.h"

namespace hitnlls {

void FactorBase::Robustify() {
    if (kernel_ == nullptr) {
        robust_scales_ << this->GetChi(), 1.0, 0.0;
        return;
    } else {
        robust_scales_ = kernel_->Robustify(this->GetChi());
    }
}

} // namespace hitnlls