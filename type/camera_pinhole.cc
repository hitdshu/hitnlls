#include "type/camera_pinhole.h"

namespace hitnlls {

matrix::Vector3d CameraPinhole::Pixel2Ray(const matrix::Vector2d &pixel) const {
    matrix::Vector3d ray;
    ray[0] = (pixel[0] - param_.cx) / param_.fx;
    ray[1] = (pixel[1] - param_.cy) / param_.fy;
    ray[2] = 1;
    return ray;
}

std::vector<matrix::Vector3d> CameraPinhole::Pixels2Rays(const std::vector<matrix::Vector2d> &pixels) const {
    std::vector<matrix::Vector3d> rays;
    for (size_t idx = 0; idx < pixels.size(); ++idx) {
        rays.push_back(Pixel2Ray(pixels[idx]));
    }
    return rays;
}

matrix::Vector2d CameraPinhole::Project(const matrix::Vector3d &point) const {
    matrix::Vector2d pixel;
    pixel[0] = param_.fx * point[0] / point[2] + param_.cx;
    pixel[1] = param_.fy * point[1] / point[2] + param_.cy;
    return pixel;
}

std::vector<matrix::Vector2d> CameraPinhole::Project(const std::vector<matrix::Vector3d> &points) const {
    std::vector<matrix::Vector2d> pixels;
    for (int idx = 0; idx < points.size(); ++idx) {
        pixels.push_back(Project(points[idx]));
    }
    return pixels;
}

JetVector2d CameraPinhole::Project(const JetVector3d &point) const {
    using geometry::Jetd;
    JetVector2d pixel;
    pixel[0] = Jetd(param_.fx) * point[0] / point[2] + Jetd(param_.cx);
    pixel[1] = Jetd(param_.fy) * point[1] / point[2] + Jetd(param_.cy);
    return pixel;
}

matrix::Matrix23d CameraPinhole::Derivative(const matrix::Vector3d &point) const {
    matrix::Matrix23d deriv;
    deriv(0, 0) = param_.fx / point[2];
    deriv(0, 1) = 0;
    deriv(0, 2) = - param_.fx * point[0] / (point[2] * point[2]);
    deriv(1, 0) = 0;
    deriv(1, 1) = param_.fy / point[2];
    deriv(1, 2) = - param_.fy * point[1] / (point[2] * point[2]);
    return deriv;
}

HITNLLS_REGISTER_CAMERA(CameraPinhole)

} // namespace hitnlls