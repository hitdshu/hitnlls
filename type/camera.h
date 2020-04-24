#pragma once

#include <vector>
#include <string>
#include <memory>

#include "common/register.h"
#include "matrix/dense.h"

namespace hitnlls {
    
struct CameraParam {
    int type = -1;
    float fx;
    float fy;
    float cx;
    float cy;
    std::vector<double> poly;
    std::vector<double> invp;
    float c;
    float d;
    float e;
};

class CameraBase {
public:
    CameraBase() = default;
    virtual ~CameraBase() = default;

    virtual bool Init(const CameraParam &param) = 0;

    virtual matrix::Vector3d Pixel2Ray(const matrix::Vector2d &pixel) const = 0;
    virtual std::vector<matrix::Vector3d> Pixels2Rays(const std::vector<matrix::Vector2d> &pixels) const = 0;

    virtual matrix::Vector2d Project(const matrix::Vector3d &point) const = 0;
    virtual std::vector<matrix::Vector2d> Project(const std::vector<matrix::Vector3d> &points) const = 0;

    virtual matrix::Matrix23d Derivative(const matrix::Vector3d &point) const = 0;

    virtual void Print() const = 0;
    virtual std::string Name() const = 0;

    CameraBase(const CameraBase &) = delete;
    CameraBase &operator=(const CameraBase &) = delete;
    CameraBase &operator=(CameraBase &&) = delete;
};

HITNLLS_REGISTER_REGISTER(CameraBase)
#define HITNLLS_REGISTER_CAMERA(name) \
    HITNLLS_REGISTER_CLASS(CameraBase, name)

} // namespace hitnlls