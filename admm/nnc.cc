#include "nnc.h"

namespace nlls {

void NncX::Project() {
    using std::max;
    VectorXf val = BaseType::GetValue();
    for (int idx = 0; idx < val.Size(); ++idx) {
        val[idx] = max<float>(0, val[idx]);
    }
    BaseType::SetValue(val);
}

} // namespace nlls