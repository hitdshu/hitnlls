#include "factor/camera_fisheye.h"
#include <cmath>

namespace hitnlls {
namespace factor {
namespace {
template<typename ValueType>
float EvalPoly(ValueType *poly, int order, float value) {
    float result;
    result = 0;
    for (int idx = 0; idx < order; ++idx) {
        result *= value;
        result += poly[order - 1 - idx];
    }
    return result;
}
}

::hitnlls::matrix::Vector3f CameraFisheye::Pixel2Ray(const ::hitnlls::matrix::Vector2f &pixel) const {
    ::hitnlls::matrix::Vector2f pixel_cc;
    pixel_cc[0] = pixel[0] - param_.cx;
    pixel_cc[1] = pixel[1] - param_.cy;
    ::hitnlls::matrix::Matrix22f mat_k;
    mat_k(0, 0) = param_.c;
    mat_k(0, 1) = param_.d;
    mat_k(1, 0) = param_.e;
    mat_k(1, 1) = 1;
    ::hitnlls::matrix::Vector2f sc = mat_k.Inverse() * pixel_cc;
    float xs = sc[0];
    float ys = sc[1];
    float pho = sqrt(xs * xs + ys * ys);
    ::hitnlls::matrix::Vector3f ray;
    ray[0] = xs;
    ray[1] = ys;
    ray[2] = EvalPoly(param_.poly.get(), param_.poly_deg, pho);
    ray[2] = fabs(ray[2]);
    ray[0] = ray[0] / ray[2];
    ray[1] = ray[1] / ray[2];
    ray[2] = ray[2] / ray[2];
    return ray;
}

::std::vector<::hitnlls::matrix::Vector3f> CameraFisheye::Pixels2Rays(const std::vector<::hitnlls::matrix::Vector2f> &pixels) const {
    ::std::vector<::hitnlls::matrix::Vector3f> rays;
    for (size_t idx = 0; idx < pixels.size(); ++idx) {
        rays.push_back(Pixel2Ray(pixels[idx]));
    }
    return rays;
}

::hitnlls::matrix::Vector2f CameraFisheye::Project(const ::hitnlls::matrix::Vector3f &point) const {
    float x = point[0];
    float y = point[1];
    float z = point[2];
    float theta = atan(z / sqrt(x * x + y * y));
    float cosphi = x / sqrt(x * x + y * y);
    float sinphi = y / sqrt(x * x + y * y);
    float pho = EvalPoly(param_.invp.get(), param_.invp_deg, -theta);
    ::hitnlls::matrix::Vector2f pixel;
    pixel[0] = param_.c * pho * cosphi + param_.d * pho * sinphi + param_.cx;
    pixel[1] = param_.e * pho * cosphi + pho * sinphi + param_.cy;
    return pixel;
}

::std::vector<::hitnlls::matrix::Vector2f> CameraFisheye::Project(const std::vector<::hitnlls::matrix::Vector3f> &points) const {
    ::std::vector<::hitnlls::matrix::Vector2f> pixels;
    for (int idx = 0; idx < points.size(); ++idx) {
        pixels.push_back(Project(points[idx]));
    }
    return pixels;
}

::hitnlls::matrix::Matrix<float, 2, 3> CameraFisheye::Derivative(const ::hitnlls::matrix::Vector3f &point) const {
    float x = point[0];
    float y = point[1];
    float z = point[2];
    float theta = atan(z / sqrt(x * x + y * y));
    float cosphi = x / sqrt(x * x + y * y);
    float sinphi = y / sqrt(x * x + y * y);
    float pho = EvalPoly(param_.invp.get(), param_.invp_deg, -theta);
    float fpho = EvalPoly(param_.poly.get(), param_.poly_deg, pho);
    float dpho = (pho * pho + fpho * fpho) / 
        (- param_.poly.get()[0] + param_.poly.get()[1] * pho * pho + 2 * param_.poly.get()[3] * pho * pho * pho
        + 3 * param_.poly.get()[4] * pho * pho * pho * pho);
    float dadp11 = (x * z) / ((sqrt(x * x + y * y)) * (x * x + y * y) + z * z * sqrt(x * x + y * y));
    float dadp12 = (y * z) / ((sqrt(x * x + y * y)) * (x * x + y * y) + z * z * sqrt(x * x + y * y));
    float dadp13 = - (sqrt(x * x + y * y)) / (x * x + y * y + z * z);
    float dadp21 = - y / ( x * x + y * y);
    float dadp22 = x / (x * x + y * y);
    float dadp23 = 0;
    ::hitnlls::matrix::Matrix<float, 2, 3> deriv;
    deriv(0, 0) = dpho * cosphi * dadp11 - pho * sinphi * dadp21;
    deriv(0, 1) = dpho * cosphi * dadp12 - pho * sinphi * dadp22;
    deriv(0, 2) = dpho * cosphi * dadp13 - pho * sinphi * dadp23;
    deriv(1, 0) = dpho * sinphi * dadp11 + pho * cosphi * dadp21;
    deriv(1, 1) = dpho * sinphi * dadp12 + pho * cosphi * dadp22;
    deriv(1, 2) = dpho * sinphi * dadp13 + pho * cosphi * dadp23;
    return deriv;
}

HITNLLS_REGISTER_CAMERA(CameraFisheye);

} // namespace factor
} // namespace hitnlls