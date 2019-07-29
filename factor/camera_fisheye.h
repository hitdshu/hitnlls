#pragma one

#include "factor/base_camera.h"
#include <iostream>

namespace hitnlls {
namespace factor {

class CameraFisheye : public BaseCamera {
    virtual bool Init(CameraParam &param) override final {
        if (1 != param.type) {
            return false;
        } else {
            param_.type = param.type;
            param_.cx = param.cx;
            param_.cy = param.cy;
            param_.c = param.c;
            param_.d = param.d;
            param_.e = param.e;
            param_.poly_deg = param.poly_deg;
            param_.invp_deg = param.invp_deg;
            param_.poly.reset(param.poly.release());
            param_.invp.reset(param.invp.release());
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
        ::std::cout << "Fisheye camera with parameters: " << ::std::endl;
        ::std::cout << "    cx " << param_.cx << ", cy " << param_.cy << ", c " << param_.c << ", d " << param_.d << ", e " << param_.e << ::std::endl;
        ::std::cout << "    poly";
        for (int pidx = 0; pidx < param_.poly_deg; ++pidx) {
            ::std::cout << " " << param_.poly[pidx];
        }
        ::std::cout << ::std::endl;
        ::std::cout << "    inv poly";
        for (int pidx = 0; pidx < param_.invp_deg; ++pidx) {
            ::std::cout << " " << param_.invp[pidx];
        }
        ::std::cout << ::std::endl;
        ::std::cout << "***********************************************" << ::std::endl;
    }
    virtual ::std::string Name() const final {
        return ::std::string("CameraFisheye");
    }

protected:
    CameraParam param_;
};

} // namesapce factor
} // namespace hitnlls