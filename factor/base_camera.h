#pragma once

#include <vector>
#include <string>
#include <memory>
#include "common/register.h"
#include "matrix/matrix.h"

namespace hitnlls {
namespace factor {

struct CameraParam {
    int type;
    float fx;
    float fy;
    float cx;
    float cy;
    ::std::unique_ptr<float []> poly;
    ::std::unique_ptr<float []> invp;
    int poly_deg;
    int invp_deg;
    float c;
    float d;
    float e;
    CameraParam() {
        type = -1;
        poly_deg = -1;
        invp_deg = -1;
    }
};

class BaseCamera {
public:
    BaseCamera() = default;
    virtual ~BaseCamera() = default;

    virtual bool Init(const CameraParam &param) = 0;

    virtual ::hitnlls::matrix::Vector3f Pixel2Ray(const ::hitnlls::matrix::Vector2f &pixel) const = 0;
    virtual ::std::vector<::hitnlls::matrix::Vector3f> Pixels2Rays(const std::vector<::hitnlls::matrix::Vector2f> &pixels) const = 0;

    virtual ::hitnlls::matrix::Vector2f Project(const ::hitnlls::matrix::Vector3f &point) const = 0;
    virtual ::std::vector<::hitnlls::matrix::Vector2f> Project(const std::vector<::hitnlls::matrix::Vector3f> &points) const = 0;

    virtual ::hitnlls::matrix::Matrix<float, 2, 3> Derivative(const ::hitnlls::matrix::Vector3f &point) const = 0;

    virtual void Print() const = 0;
    virtual ::std::string Name() const = 0;

    BaseCamera(const BaseCamera &) = delete;
    BaseCamera &operator=(const BaseCamera &) = delete;
};

HITNLLS_REGISTER_REGISTERER(BaseCamera);
#define HITNLLS_REGISTER_CAMERA(name) \
    HITNLLS_REGISTER_CLASS(BaseCamera, name)

} // namespace factor
} // namespace hitnlls