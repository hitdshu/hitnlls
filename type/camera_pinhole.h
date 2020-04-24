#pragma once

#include "type/camera.h"
#include "geometry/jet.h"

namespace hitnlls {

using JetVector2d = matrix::Matrix<geometry::Jetd, 2, 1>;
using JetVector3d = matrix::Matrix<geometry::Jetd, 3, 1>;

class CameraPinhole : public CameraBase {
public:
    virtual bool Init(const CameraParam &param) override final {
        if (0 != param.type) {
            return false;
        } else {
            param_.type = param.type;
            param_.fx = param.fx;
            param_.fy = param.fy;
            param_.cx = param.cx;
            param_.cy = param.cy;
            return true;
        }
    }

    virtual matrix::Vector3d Pixel2Ray(const matrix::Vector2d &pixel) const override final;
    virtual std::vector<matrix::Vector3d> Pixels2Rays(const std::vector<matrix::Vector2d> &pixels) const override final;

    virtual matrix::Vector2d Project(const matrix::Vector3d &point) const override final;
    virtual std::vector<matrix::Vector2d> Project(const std::vector<matrix::Vector3d> &points) const override final;
    JetVector2d Project(const JetVector3d &point) const;

    virtual matrix::Matrix23d Derivative(const matrix::Vector3d &point) const override final;

    virtual void Print() const override final {
        std::cout << "***********************************************" << std::endl;
        std::cout << "Pinhole camera with parameters: " << std::endl;
        std::cout << "    fx " << param_.fx << ", fy " << param_.fy << ", cx " << param_.cx << ", cy " << param_.cy << std::endl;
        std::cout << "***********************************************" << std::endl;
    }
    virtual std::string Name() const override final {
        return std::string("CameraPinhole");
    }

protected:
    CameraParam param_;
};

} // namespace hitnlls