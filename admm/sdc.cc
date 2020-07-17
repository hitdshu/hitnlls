#include "sdc.h"

namespace nlls {

void SdcX::Project() {
    using std::max;
    MatrixXf val = BaseType::GetValue();
    int n = val.Rows();
    EVD<MatrixXf> solver(val);
    VectorXf eigen_vals = solver.V();
    MatrixXf eigen_vecs = solver.U();
    for (int idx = 0; idx < n; ++idx) {
        eigen_vals[idx] = max<float>(0.0, eigen_vals[idx]);
    }
    MatrixXf eigen_vals_mat(n, n);
    eigen_vals_mat.SetIdentity();
    for (int idx = 0; idx < n; ++idx) {
        eigen_vals_mat(idx, idx) = eigen_vals[idx];
    }
    val = eigen_vecs * eigen_vals_mat * eigen_vecs.Transpose();
    BaseType::SetValue(val);
}

} // namespace nlls