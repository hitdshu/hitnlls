#pragma once

#include "type/camera.h"

namespace hitnlls {

class CameraFisheye : public CameraBase {
public:
    virtual bool Init(const CameraParam &param) override final {
        if (1 != param.type) {
            return false;
        } else {
            param_.type = param.type;
            param_.cx = param.cx;
            param_.cy = param.cy;
            param_.c = param.c;
            param_.d = param.d;
            param_.e = param.e;
            param_.poly = param.poly;
            param_.invp = param.invp;
            return true;
        }
    }

    virtual matrix::Vector3d Pixel2Ray(const matrix::Vector2d &pixel) const override final;
    virtual std::vector<matrix::Vector3d> Pixels2Rays(const std::vector<matrix::Vector2d> &pixels) const override final;

    virtual matrix::Vector2d Project(const matrix::Vector3d &point) const override final;
    virtual std::vector<matrix::Vector2d> Project(const std::vector<matrix::Vector3d> &points) const override final;

    virtual matrix::Matrix23d Derivative(const matrix::Vector3d &point) const override final;

    virtual void Print() const override final {
        std::cout << "***********************************************" << std::endl;
        std::cout << "Fisheye camera with parameters: " << std::endl;
        std::cout << "    cx " << param_.cx << ", cy " << param_.cy << ", c " << param_.c << ", d " << param_.d << ", e " << param_.e << std::endl;
        std::cout << "    poly";
        for (int pidx = 0; pidx < param_.poly.size(); ++pidx) {
            std::cout << " " << param_.poly[pidx];
        }
        std::cout << std::endl;
        std::cout << "    inv poly";
        for (int pidx = 0; pidx < param_.invp.size(); ++pidx) {
            std::cout << " " << param_.invp[pidx];
        }
        std::cout << std::endl;
        std::cout << "***********************************************" << ::std::endl;
    }
    virtual std::string Name() const override final {
        return std::string("CameraFisheye");
    }

protected:
    CameraParam param_;
};

} // namespace hitnlls