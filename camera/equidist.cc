#include "equidist.h"

namespace nlls {

Vector3f Equidist::Undistort(const Vector3f &pt, const float *coeff) {
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

Vector3f Equidist::Distort(const Vector3f &pt, const float *coeff, Matrix33f *dddp, Matrix34f *dddc) {
    const float k1 = coeff[0];
    const float k2 = coeff[1];
    const float k3 = coeff[2];
    const float k4 = coeff[3];
    float m = pt[0] / pt[2];
    float n = pt[1] / pt[2];
    float r = sqrt(m * m + n * n);
    float r2 = r * r;
    float theta = atan2(r2, 1.0);
    float theta2 = theta * theta;
    float theta3 = theta2 * theta;
    float theta4 = theta2 * theta2;
    float theta5 = theta4 * theta;
    float theta6 = theta2 * theta4;
    float theta7 = theta6 * theta;
    float theta8 = theta4 * theta4;
    float theta9 = theta8 * theta;
    float theta_d = theta * (1 + k1 * theta2 + k2 * theta4 + k3 * theta6 + k4 * theta8);
    Vector3f result;
    float cos_phi = m / r;
    float sin_phi = n / r;
    result << theta_d * cos_phi, theta_d * sin_phi, 1;
    if (dddp) {
        Matrix23f dmndp;
        dmndp << 1 / pt[2], 0, - pt[0] / pow(pt[2], 2), 
            0, 1 / pt[2], - pt[1] / pow(pt[2], 2);
        Vector2f dr2dmn(2 * m, 2 * n);
        Vector2f drdmn(m / r, n / r);
        float dtdr2 = 1.0 / (1.0 + r2 * r2);
        float dtddt = 1 + k1 * 3 * theta2 + k2 * 5 * theta4 + k3 * 7 * theta6 + k4 * 9 * theta8;
        Vector2f dtddmn = dtddt * dtdr2 * dr2dmn;
        Vector2f dcdmn;
        dcdmn << 1 / r - m / r2 * drdmn[0], -m / r2 * drdmn[1];
        Vector2f dsdmn;
        dsdmn << -n / r2 * drdmn[0], 1 / r - n / r2 * drdmn[1];
        Matrix22f dnmndmn;
        dnmndmn << dtddmn[0] * cos_phi + theta_d * dcdmn[0], dtddmn[1] * cos_phi + theta_d * dcdmn[1], 
            dtddmn[0] * sin_phi + theta_d * dsdmn[0], dtddmn[1] * sin_phi + theta_d * dsdmn[1];
        dddp->Block(0, 0, 2, 3) = dnmndmn * dmndp;
        dddp->Row(2) = 0;
    }
    if (dddc) {
        (*dddc) << theta3 * cos_phi,  theta5 * cos_phi, theta7 * cos_phi, theta9 * cos_phi, 
            theta3 * sin_phi,  theta5 * sin_phi, theta7 * sin_phi, theta9 * sin_phi, 
            0, 0, 0, 0;
    }
    return result;
}

} // namespace nlls