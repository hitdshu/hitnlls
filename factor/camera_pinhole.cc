#include "factor/camera_pinhole.h"

namespace hitnlls {
namespace factor {

::hitnlls::matrix::Vector3f CameraPinhole::Pixel2Ray(const ::hitnlls::matrix::Vector2f &pixel) const {
    ::hitnlls::matrix::Vector3f ray;
    ray[0] = (pixel[0] - param_.cx) / param_.fx;
    ray[1] = (pixel[1] - param_.cy) / param_.fy;
    ray[2] = 1;
    return ray;
}

::std::vector<::hitnlls::matrix::Vector3f> CameraPinhole::Pixels2Rays(const std::vector<::hitnlls::matrix::Vector2f> &pixels) const {
    ::std::vector<::hitnlls::matrix::Vector3f> rays;
    for (size_t idx = 0; idx < pixels.size(); ++idx) {
        rays.push_back(Pixel2Ray(pixels[idx]));
    }
    return rays;
}

::hitnlls::matrix::Vector2f CameraPinhole::Project(const ::hitnlls::matrix::Vector3f &point) const {
    ::hitnlls::matrix::Vector2f pixel;
    pixel[0] = param_.fx * point[0] / point[2] + param_.cx;
    pixel[1] = param_.fy * point[1] / point[2] + param_.cy;
    return pixel;
}

::std::vector<::hitnlls::matrix::Vector2f> CameraPinhole::Project(const std::vector<::hitnlls::matrix::Vector3f> &points) const {
    ::std::vector<::hitnlls::matrix::Vector2f> pixels;
    for (int idx = 0; idx < points.size(); ++idx) {
        pixels.push_back(Project(points[idx]));
    }
    return pixels;
}

::hitnlls::matrix::Matrix<float, 2, 3> CameraPinhole::Derivative(const ::hitnlls::matrix::Vector3f &point) const {
    ::hitnlls::matrix::Matrix<float, 2, 3> deriv;
    deriv(0, 0) = param_.fx / point[2];
    deriv(0, 1) = 0;
    deriv(0, 2) = - param_.fx * point[0] / (point[2] * point[2]);
    deriv(1, 0) = 0;
    deriv(1, 1) = param_.fy / point[2];
    deriv(1, 2) = - param_.fy * point[1] / (point[2] * point[2]);
    return deriv;
}

HITNLLS_REGISTER_CAMERA(CameraPinhole);

} // namespace factor
} // namespace hitnlls