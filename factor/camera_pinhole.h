#pragma one

#include "factor/base_camera.h"
#include <iostream>

namespace hitnlls {
namespace factor {

class CameraPinhole : public BaseCamera {
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

    virtual ::hitnlls::matrix::Vector3f Pixel2Ray(const ::hitnlls::matrix::Vector2f &pixel) const override final;
    virtual ::std::vector<::hitnlls::matrix::Vector3f> Pixels2Rays(const std::vector<::hitnlls::matrix::Vector2f> &pixels) const override final;

    virtual ::hitnlls::matrix::Vector2f Project(const ::hitnlls::matrix::Vector3f &point) const override final;
    virtual ::std::vector<::hitnlls::matrix::Vector2f> Project(const std::vector<::hitnlls::matrix::Vector3f> &points) const override final;

    virtual ::hitnlls::matrix::Matrix<float, 2, 3> Derivative(const ::hitnlls::matrix::Vector3f &point) const override final;

    virtual void Print() const override final {
        ::std::cout << "***********************************************" << ::std::endl;
        ::std::cout << "Pinhole camera with parameters: " << ::std::endl;
        ::std::cout << "    fx " << param_.fx << ", fy " << param_.fy << ", cx " << param_.cx << ", cy " << param_.cy << ::std::endl;
        ::std::cout << "***********************************************" << ::std::endl;
    }
    virtual ::std::string Name() const final {
        return ::std::string("CameraPinhole");
    }

protected:
    CameraParam param_;
};

} // namesapce factor
} // namespace hitnlls