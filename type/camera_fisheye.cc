#include "type/camera_fisheye.h"

namespace hitnlls {
namespace {
template<typename ValueType>
double EvalPoly(ValueType *poly, int order, double value) {
    double result;
    result = 0;
    for (int idx = 0; idx < order; ++idx) {
        result *= value;
        result += poly[order - 1 - idx];
    }
    return result;
}
} // namespace 

matrix::Vector3d CameraFisheye::Pixel2Ray(const matrix::Vector2d &pixel) const {
    matrix::Vector2d pixel_cc;
    pixel_cc[0] = pixel[0] - param_.cx;
    pixel_cc[1] = pixel[1] - param_.cy;
    matrix::Matrix22d mat_k;
    mat_k(0, 0) = param_.c;
    mat_k(0, 1) = param_.d;
    mat_k(1, 0) = param_.e;
    mat_k(1, 1) = 1;
    matrix::Vector2d sc = mat_k.Inverse() * pixel_cc;
    double xs = sc[0];
    double ys = sc[1];
    double pho = std::sqrt(xs * xs + ys * ys);
    matrix::Vector3d ray;
    ray[0] = xs;
    ray[1] = ys;
    ray[2] = EvalPoly(param_.poly.data(), param_.poly.size(), pho);
    ray[2] = std::abs(ray[2]);
    ray[0] = ray[0] / ray[2];
    ray[1] = ray[1] / ray[2];
    ray[2] = ray[2] / ray[2];
    return ray;
}

std::vector<matrix::Vector3d> CameraFisheye::Pixels2Rays(const std::vector<matrix::Vector2d> &pixels) const {
    std::vector<matrix::Vector3d> rays;
    for (size_t idx = 0; idx < pixels.size(); ++idx) {
        rays.push_back(Pixel2Ray(pixels[idx]));
    }
    return rays;
}

matrix::Vector2d CameraFisheye::Project(const matrix::Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double theta = std::atan(z / std::sqrt(x * x + y * y));
    double cosphi = x / std::sqrt(x * x + y * y);
    double sinphi = y / std::sqrt(x * x + y * y);
    double pho = EvalPoly(param_.invp.data(), param_.invp.size(), -theta);
    matrix::Vector2d pixel;
    pixel[0] = param_.c * pho * cosphi + param_.d * pho * sinphi + param_.cx;
    pixel[1] = param_.e * pho * cosphi + pho * sinphi + param_.cy;
    return pixel;
}

std::vector<matrix::Vector2d> CameraFisheye::Project(const std::vector<matrix::Vector3d> &points) const {
    std::vector<matrix::Vector2d> pixels;
    for (int idx = 0; idx < points.size(); ++idx) {
        pixels.push_back(Project(points[idx]));
    }
    return pixels;
}

matrix::Matrix23d CameraFisheye::Derivative(const matrix::Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double theta = std::atan(z / std::sqrt(x * x + y * y));
    double cosphi = x / std::sqrt(x * x + y * y);
    double sinphi = y / std::sqrt(x * x + y * y);
    double pho = EvalPoly(param_.invp.data(), param_.invp.size(), -theta);
    double fpho = EvalPoly(param_.poly.data(), param_.poly.size(), pho);
    double dpho = (pho * pho + fpho * fpho) / 
        (- param_.poly[0] + param_.poly[1] * pho * pho + 2 * param_.poly[3] * pho * pho * pho
        + 3 * param_.poly[4] * pho * pho * pho * pho);
    double dadp11 = (x * z) / ((sqrt(x * x + y * y)) * (x * x + y * y) + z * z * sqrt(x * x + y * y));
    double dadp12 = (y * z) / ((sqrt(x * x + y * y)) * (x * x + y * y) + z * z * sqrt(x * x + y * y));
    double dadp13 = - (sqrt(x * x + y * y)) / (x * x + y * y + z * z);
    double dadp21 = - y / ( x * x + y * y);
    double dadp22 = x / (x * x + y * y);
    double dadp23 = 0;
    matrix::Matrix23d deriv;
    deriv(0, 0) = dpho * cosphi * dadp11 - pho * sinphi * dadp21;
    deriv(0, 1) = dpho * cosphi * dadp12 - pho * sinphi * dadp22;
    deriv(0, 2) = dpho * cosphi * dadp13 - pho * sinphi * dadp23;
    deriv(1, 0) = dpho * sinphi * dadp11 + pho * cosphi * dadp21;
    deriv(1, 1) = dpho * sinphi * dadp12 + pho * cosphi * dadp22;
    deriv(1, 2) = dpho * sinphi * dadp13 + pho * cosphi * dadp23;
    matrix::Matrix22d aff;
    aff << param_.c, param_.d, param_.e, 1;
    return aff * deriv;
}

HITNLLS_REGISTER_CAMERA(CameraFisheye)

} // namespace hitnlls