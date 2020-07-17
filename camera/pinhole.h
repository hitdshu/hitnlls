#pragma once

#include "camera.h"
#include "distortion.h"

namespace nlls {

template <class D>
class Pinhole : public Camera<4, D> {
public:
    using ThisType = Pinhole;
    using BaseType = Camera<4, D>;
    static constexpr int CDOF = 4;
    static constexpr int DDOF = D::DOF;
    static constexpr int DOF = BaseType::DOF;
    using CameraParamVector = typename BaseType::CameraParamVector;
    using TotalParamVector = typename BaseType::TotalParamVector;

    explicit Pinhole() : BaseType() {}
    explicit Pinhole(int width, int height) : BaseType(width, height) {}
    explicit Pinhole(int width, int height, const TotalParamVector &tp) : BaseType(width, height, tp) {}

    float &Fx() { return BaseType::CameraData()[0]; }
    const float &Fx() const { return BaseType::CameraData()[0]; }
    float &Fy() { return BaseType::CameraData()[1]; }
    const float &Fy() const { return BaseType::CameraData()[1]; }
    float &Cx() { return BaseType::CameraData()[2]; }
    const float &Cx() const { return BaseType::CameraData()[2]; }
    float &Cy() { return BaseType::CameraData()[3]; }
    const float &Cy() const { return BaseType::CameraData()[3]; }

    virtual Vector2f Project(const Vector3f &pt3d, Matrix23f *dp = nullptr, Matrix<float, 2, DOF> *dc = nullptr) const override final {
        Vector2f result;
        if (dp && dc) {
            Matrix33f dddp;
            Matrix<float, 3, DDOF> dddd;
            Vector3f tmp = D::Distort(pt3d, BaseType::DistortionData(), &dddp, &dddd);
            result << Fx()*tmp[0]/tmp[2]+Cx(), Fy()*tmp[1]/tmp[2]+Cy();
            Matrix23f dpdd;
            dpdd << Fx()/tmp[2], 0, -Fx()*tmp[0]/nlls::pow(tmp[2], 2), 
                0, Fy()/tmp[2], -Fy()*tmp[1]/nlls::pow(tmp[2], 2);
            Matrix<float, 2, CDOF> dpdc;
            dpdc << tmp[0]/tmp[2], 0, 1, 0, 
                0, tmp[1]/tmp[2], 0, 1;
            *dp = dpdd * dddp;
            dc->Block(0, CDOF, 2, DDOF) = dpdd * dddd;
            dc->Block(0, 0, 2, CDOF) = dpdc;
        } else if (dp) {
            Matrix33f dddp;
            Vector3f tmp = D::Distort(pt3d, BaseType::DistortionData(), &dddp, nullptr);
            result << Fx()*tmp[0]/tmp[2]+Cx(), Fy()*tmp[1]/tmp[2]+Cy();
            Matrix23f dpdd;
            dpdd << Fx()/tmp[2], 0, -Fx()*tmp[0]/nlls::pow(tmp[2], 2), 
                0, Fy()/tmp[2], -Fy()*tmp[1]/nlls::pow(tmp[2], 2);
            *dp = dpdd * dddp;
        } else if (dc) {
            Matrix<float, 3, DDOF> dddd;
            Vector3f tmp = D::Distort(pt3d, BaseType::DistortionData(), nullptr, &dddd);
            result << Fx()*tmp[0]/tmp[2]+Cx(), Fy()*tmp[1]/tmp[2]+Cy();
            Matrix23f dpdd;
            dpdd << Fx()/tmp[2], 0, -Fx()*tmp[0]/nlls::pow(tmp[2], 2), 
                0, Fy()/tmp[2], -Fy()*tmp[1]/nlls::pow(tmp[2], 2);
            dc->Block(0, CDOF, 2, DDOF) = dpdd * dddd;
            Matrix<float, 2, CDOF> dpdc;
            dpdc << tmp[0]/tmp[2], 0, 1, 0, 
                0, tmp[1]/tmp[2], 0, 1;
            dc->Block(0, 0, 2, CDOF) = dpdc;
        } else {
            Vector3f tmp = D::Distort(pt3d, BaseType::DistortionData());
            result << Fx()*tmp[0]/tmp[2]+Cx(), Fy()*tmp[1]/tmp[2]+Cy();
        }
        return result;
    }
    virtual Vector3f UnProject(const Vector2f &pt2d) const override final {
        Vector3f ppt;
        ppt << (pt2d[0]-Cx())/Fx(), (pt2d[1]-Cy())/Fy(), 1;
        return D::Undistort(ppt, BaseType::DistortionData());
    }
};

using PinholeIdeal = Pinhole<Ideal>;
using PinholeRadtan = Pinhole<RadTan>;
using PinholeFisheye = Pinhole<Equidist>;

} // namespace nlls