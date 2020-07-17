#include "radtan.h"

namespace nlls {

Vector3f RadTan::Undistort(const Vector3f &pt, const float *coeff) {
    Vector3f target = pt / pt[2];
    Vector3f result = target;
    const int n = 10;
    const double eps = 1e-10;
    for (int i = 0; i < n; ++i) {
        Matrix33f dddp;
        Vector3f tmp = Distort(result, coeff, &dddp);
        Vector2f err((target - tmp).Block(0, 0, 2, 1));
        if (err.Norm() < eps) {
            break;
        }
        Matrix22f dddp_tmp = dddp.Block(0, 0, 2, 2);
        QR<Matrix22f> qr(dddp_tmp);
        Vector2f inc = qr.Solve(err);
        result.Block(0, 0, 2, 2) += inc;
    }
    return result;
}

Vector3f RadTan::Distort(const Vector3f &pt, const float *coeff, Matrix33f *dddp, Matrix34f *dddc) {
    const float k1 = coeff[0];
    const float k2 = coeff[1];
    const float p1 = coeff[2];
    const float p2 = coeff[3];
    float m = pt[0] / pt[2];
    float n = pt[1] / pt[2];
    float r = sqrt(m * m + n * n);
    float r2 = r * r;
    float r4 = r2 * r2;
    float rm = 1 + k1 * r2 + k2 * r4;
    float nm = m * rm + 2 * p1 * m * n + p2 * (r2 + 2 * m * m);
    float nn = n * rm + p1 * (r2 + 2 * n * n) + 2 * p2 * m * n;
    Vector3f result;
    result << nm, nn, 1;
    if (dddp) {
        Matrix23f dmndp;
        dmndp << 1 / pt[2], 0, -pt[0] / pow(pt[2], 2), 
            0, 1 / pt[2], -pt[1] / pow(pt[2], 2);
        Vector2f dr2dmn(2 * m, 2 * n);
        Vector2f dr4dmn = r * r * dr2dmn;
        Vector2f drmdmn = k1 * dr2dmn + k2 * dr4dmn;
        Matrix22f dnmndmn;
        dnmndmn << rm + m * drmdmn[0] + 2 * p1 * n + p2 * (dr2dmn[0] + 4 * m), m * drmdmn[1] + 2 * p1 * m + p2 * dr2dmn[1], 
            n * drmdmn[1] + p1 * dr2dmn[0] + 2 * p2 * n, rm + n * drmdmn[1] + p1 * (dr2dmn[1] + 4 * n) + 2 * p2 * m;
        dddp->Block(0, 0, 2, 3) = dnmndmn * dmndp;
        dddp->Row(2) = 0;
    }
    if (dddc) {
        (*dddc) << m * r2, m * r4, 2 * m * n, r2 + 2 * m * m, 
            n * r2, n * r4, r2 + 2 * n * n, 2 * m * n, 
            0, 0, 0, 0;
    }
    return result;
}

} // namespace nlls