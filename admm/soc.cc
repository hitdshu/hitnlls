#include "soc.h"

namespace nlls {

void SocX::Project() {
    int n = BaseType::GetDim();
    VectorXf val = BaseType::GetValue();
    VectorXf v = val.Block(0, 0, n - 1, 1);
    float t = val[n - 1];
    float vn = v.Norm();
    if (t <= -vn) {
        val.SetZero();
    } else if (t >= vn) {
    } else {
        float a = (vn + t) / 2;
        val.Block(0, 0, n - 1, 1) *= a / vn;
        val[n - 1] = a;
    }
    BaseType::SetValue(val);
}

} // namespace nlls